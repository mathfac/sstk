'use strict';

var AssetLoader = require('assets/AssetLoader');
var Constants = require('Constants');
var PubSub = require('PubSub');
var Autocomplete = require('ui/Autocomplete');
var SolrQueryParser = require('search/SolrQueryParser');
var UIUtil = require('ui/UIUtil');
var _ = require('util');
require('jquery-pagination');

/**
 * Panel that shows search results in a grid layout (uses jquery-pagination for pagination)
 * @param params Panel configuration
 * @param [params.entriesPerRow=3] {int} Number of results to show per row
 * @param [params.nRows=33] {int} Limit on number of rows to query/display
 * @param [params.padding=50] {int} Amount of padding (in pixels) to have between result entries
 * @constructor
 * @memberOf ui
 */
function SearchPanel(params) {
  PubSub.call(this);

  this.searchModule = params.searchModule;
  this.source = params.source;

  this.entriesPerRow = 3;
  this.nRows = 33;
  this.sortOrder = '';
  this.additionalSortOrder = '';
  this.padding = 50;

  this.tooltipIncludeAll = true;
  this.tooltipIncludeFields =  [];
  this.showAnimatedOnHover = true;
  this.showSearchOptions = true;
  this.showSearchSortOption = true;
  this.showSearchSourceOption = true;
  this.loadImagesLazy = false;
  this.previewImageIndex = undefined;

  if (params) {
    // Application callback when search has succeeded
    // Parameters: resultList (list of results)
    this.searchSucceededCallback = params.searchSucceededCallback;
    // Application callback when search has failed
    // Parameters: error message
    this.searchFailedCallback = params.searchFailedCallback;
    // Application callback to get url of preview image
    // Parameters: result.id
    this.getImagePreviewUrlCallback = params.getImagePreviewUrlCallback;
    // Application callback when search result if clicked on
    // Parameters: source, result.id, result
    this.onClickResultCallback = params.onClickResultCallback;
    this.onMousedownResultCallback = params.onMousedownResultCallback;
    this.onHoverResultCallback = params.onHoverResultCallback;
    // Application callback to add more stuff to result element
    this.appendResultElemCallback = params.appendResultElemCallback;
    // Application callback for when source changed
    this.sourceChangedCallback = params.sourceChangedCallback;
    // Application callback for expanding grouped nodes
    this.expandGroupCallback = params.expandGroupCallback;
    if (!this.expandGroupCallback && params.allowGroupExpansion) {
      this.expandGroupCallback = function (group) {
        if (group.fields && group.fields.length === 1) {
          var field = group.fields[0];
          this.expandCategory(field.value);
        }
      }.bind(this);
    }

    // Search panel (overall search panel)
    this.container = params.searchPanel;

    // Number of results to show per row (default is 3 columns)
    if (params.entriesPerRow !== undefined) this.entriesPerRow = params.entriesPerRow;
    // Limit on number of rows to query/display
    if (params.nRows !== undefined) this.nRows = params.nRows;
    // Sources to allow
    if (params.sources !== undefined) this.sources = params.sources;
    // If include all returned fields in tooltip or just name
    if (params.tooltipIncludeAll !== undefined) this.tooltipIncludeAll = params.tooltipIncludeAll;
    if (params.tooltipIncludeFields !== undefined) this.tooltipIncludeFields = params.tooltipIncludeFields;
    if (params.tooltipIncludeExtraFields !== undefined) {
      this.tooltipIncludeFields = this.tooltipIncludeFields.concat(params.tooltipIncludeExtraFields);
    }
    // Show search options?
    if (params.showSearchOptions !== undefined) this.showSearchOptions = params.showSearchOptions;
    if (params.showSearchSortOption !== undefined) this.showSearchSortOption = params.showSearchSortOption;
    if (params.showSearchSourceOption !== undefined) this.showSearchSourceOption = params.showSearchSourceOption;
    if (params.showSearchBySize !== undefined) this.showSearchBySize = params.showSearchBySize;
    // Show animated gif onhover?
    if (params.showAnimatedOnHover !== undefined) this.showAnimatedOnHover = params.showAnimatedOnHover;

    // Load/Save support?
    if (params.showLoadFile !== undefined) this.showLoadFile = params.showLoadFile;
    if (params.allowSave !== undefined) this.allowSave = params.allowSave;

    if (params.loadImagesLazy !== undefined) this.loadImagesLazy = params.loadImagesLazy;

    if (params.sortOrder !== undefined) this.sortOrder = params.sortOrder;
    if (params.additionalSortOrder !== undefined) this.additionalSortOrder = params.additionalSortOrder;

    if (params.previewImageIndex !== undefined) this.previewImageIndex = params.previewImageIndex;

    if (params.padding !== undefined) this.padding = params.padding;

    // Predefined elements
    this.searchButton = params.searchButton;
    this.searchTextElem = params.searchTextElem;
  }

  this.limit = this.nRows * this.entriesPerRow;
  if (this.sources.indexOf(this.source) < 0) {
    this.source = this.sources[0];
  }

  if (this.container && this.container.length) {
    var scope = this;
    if (this.showSearchOptions) {
      // Binds elements
      this.searchOptionsElem = $('<div></div>').attr('class', 'searchOptions');
      // Text field from which to get the query term
      if (!this.searchTextElem || $(this.searchTextElem).length === 0) {
        this.searchTextElem = $('<input/>').attr('type', 'text').attr('class', 'search').attr('size','32');
        this.searchOptionsElem.append(this.searchTextElem);
      }

      // Allows people to select from different sources to search

      if (this.showSearchSourceOption) {
        this.sourceElem = $('<select></select>').attr('class', 'searchSource');
        this.searchOptionsElem.append(this.sourceElem);
      }

      this.container.append(this.searchOptionsElem);

      // Add search button
      if (!this.searchButton || $(this.searchButton).length === 0) {
        this.searchButton = $('<input type="button" value="Search" class="searchButton"/>');
        this.searchOptionsElem.append(this.searchButton);
      }
      this.searchButton.click(function(event) { this.startSearch(); }.bind(this));

      // Add sort div
      if (this.showSearchSortOption) {
        var sortDiv = $('<div></div>');
        // Sort order
        this.sortOrderElem = $('<input/>').attr('type', 'text').attr('class', 'sortOrder').attr('size', '8');
        this.sortOrderDropdownElem = $('<select></select>').attr('class', 'sortOrder');
        this.sortOrderDropdownElem.attr('title', 'Specify sort order');
        this.defaultSortOrders = ['score', 'id', 'time', 'size', 'popularity', 'quality', 'random', 'custom'];
        for (var i = 0; i < this.defaultSortOrders.length; i++) {
          var s = this.defaultSortOrders[i];
          this.sortOrderDropdownElem.append('<option value="' + s + '">' + s + '</option>');
        }
        this.sortOrderDropdownElem.change(function () {
          scope.sortOrderDropdownElem.find('option:selected').each(function () {
            var sort = $(this).val();
            if (sort === 'score') {
              scope.sortOrderElem.val('');
            } else if (sort === 'id') {
              scope.sortOrderElem.val('id asc');
            } else if (sort === 'time') {
              scope.sortOrderElem.val('timestamp desc');
            } else if (sort === 'popularity') {
              scope.sortOrderElem.val('popularity desc');
            } else if (sort === 'quality') {
              scope.sortOrderElem.val('modelQuality desc');
            } else if (sort === 'size') {
              if (Constants.assetSources.scene.indexOf(scope.source) >= 0) {
                scope.sortOrderElem.val('nmodels asc');
              } else if (Constants.assetSources.model.indexOf(scope.source) >= 0 || Constants.assetSources.scan.indexOf(scope.source) >= 0) {
                scope.sortOrderElem.val('nfaces asc');
              } else {
                scope.sortOrderElem.val('');
              }
            } else if (sort === 'random') {
              scope.sortOrderElem.val(scope.searchModule.getRandomSortOrder());
            } else {
              // custom sort
              //scope.sortOrderElem.val("");
            }
            scope.search(scope.searchTextElem.val());
          });
        });
        this.sortOrderDropdownElem.val('score');
        sortDiv.append('Sort:');
        sortDiv.append(this.sortOrderDropdownElem);
        sortDiv.append(this.sortOrderElem);

        this.searchOptionsElem.append(sortDiv);
      }

      if (this.showSearchBySize) {
        this.searchBySizeElem = $('<input/>').attr('type', 'checkbox');
        var searchBySizeLabel = $('<label></label>').append(this.searchBySizeElem).append('Use size');
        this.searchOptionsElem.append(searchBySizeLabel);
      }
    }
    var buttonGroup = $('<div class="btn-group"></div>');
    this.searchOptionsButtonGroup = buttonGroup;
    if (this.allowSave) {
      var saveListButton = $('<input class="btn btn-default" type="button" value="Save Ids" class="saveButton"/>');
      saveListButton.attr('title', 'Save search results');
      saveListButton.click(function() {
        var url = scope.getLastQueryURL({ format: 'csv', fields: 'id' });
        if (url) {
          window.open(url, 'Search Results');
        } else {
          UIUtil.showAlert(scope.container, 'Please perform a search before saving ids', 'alert-warning');
        }
      });
      buttonGroup.append(saveListButton);
    }
    if (this.showLoadFile) {
      var loadFileInput = UIUtil.createFileInput({
        id: 'search-loadFile',
        label: 'Load Ids',
        hideFilename: true,
        inline: true,
        loadFn: this.loadIdsFromFile.bind(this)
      });
      loadFileInput.file.attr('title', 'File with list of ids (limited to <1000)');
      //console.log(loadFileInput);
      buttonGroup.append(loadFileInput.group);
    }
    if (this.searchOptionsElem) {
      this.searchOptionsElem.append(buttonGroup);
    }

    this.searchInfoElem = $('<div></div>').attr('class', 'searchInfo');
    // Display pagination information
    this.pageElem = $('<div></div>');
    this.searchInfoElem.append(this.pageElem);
    this.container.append(this.searchInfoElem);

    // Results are displayed in the resultsElem in a tabular form
    this.resultsElem = $('<div></div>').attr('class', 'searchResults');
    this.__setResultsElemSize();
    this.container.append(this.resultsElem);

    // Add element to go back from expanded group
    if (this.expandGroupCallback) {
      this.__unexpandGroupElem = $('<div class="glyphicon glyphicon-chevron-left">Back</div>');
      this.__unexpandGroupElem.hide();
      this.searchInfoElem.append(this.__unexpandGroupElem);
      this.__unexpandGroupElem.click(this.unexpandGroup.bind(this));
    }

    // Bind search action to when search text changes
    if (this.searchTextElem) {
      this.searchTextElem.keyup(function(event) {
        event.stopPropagation();
        if (event.keyCode === 13) {
          this.startSearch();
        }
      }.bind(this));
    }
    // Bind search source to when source source changes
    if (this.sourceElem) {
      this.sourceElem.empty();
      for (var i = 0; i < this.sources.length; i++) {
        var s = this.sources[i];
        this.sourceElem.append('<option value="' + s + '">' + s + '</option>');
      }
      this.sourceElem.change(function () {
        scope.sourceElem.find('option:selected').each(function () {
          scope.source = $(this).val();
          scope.search(scope.searchTextElem.val());
          //                    console.log("Source now " + scope.source);
        });
        if (scope.sourceChangedCallback) {
          scope.sourceChangedCallback(scope.source);
        }
      });
      this.sourceElem.val(this.source);
      this.source = this.sourceElem.val();
    }
  }
}

SearchPanel.prototype = Object.create(PubSub.prototype);
SearchPanel.prototype.constructor = SearchPanel;

Object.defineProperty(SearchPanel.prototype, 'isSearchBySize', {
  get: function () {return this.searchBySizeElem? this.searchBySizeElem.prop('checked') : false; },
  set: function (v) {
    if (this.searchBySizeElem) {
      this.searchBySizeElem.prop('checked', v);
    }
  }
});

SearchPanel.prototype.loadIdsFromFile = function(file) {
  var scope = this;
  var loader = new AssetLoader();
  loader.load(file, 'UTF-8', function(data) {
      var ids = data.split('\n');
      var options = {};//{ ordering: ids };
      scope.searchModule.queryIds(ids, function(err, res) {
        if (err) {
          scope.searchFailed(err);
        } else {
          scope.searchSucceeded(options, res);
        }
      });
    },
    undefined,
    function(err) {
      var filename = (typeof file === 'string')? file : file.name;
      UIUtil.showAlert(scope.searchPanel, 'Error loading ids from ' + filename + ': ' + err);
      console.error('Error loading ids from ' + filename, err);
    });
};

SearchPanel.prototype.setAutocomplete = function (suggester, options) {
  if (this.searchTextElem) {
    this.autocomplete = new Autocomplete({
      input: this.searchTextElem,
      suggester: suggester,
      options: options
    });
  }
};

SearchPanel.prototype.insertBeforeSearchBox = function(element) {
  this.searchTextElem.before(element);
};

SearchPanel.prototype.addToSearchOptions = function(element) {
  this.searchOptionsButtonGroup.before(element);
};

SearchPanel.prototype.setSearchText = function (s) {
  if (this.searchTextElem) {
    this.searchTextElem.val(s);
    return true;
  } else {
    this._searchText = s;
    return false;
  }
};

SearchPanel.prototype.getSearchText = function (defaultValue) {
  if (defaultValue == undefined) { defaultValue = '*:*'; }
  if (this.searchTextElem) {
    var searchTerm = this.searchTextElem.val().trim();
    if (searchTerm === '') {
      this.searchTextElem.val(defaultValue);
      searchTerm = this.searchTextElem.val();
    }
    return searchTerm;
  } else {
    if (this._searchText === '' || this._searchText == undefined) {
      this._searchText = defaultValue;
    }
    return this._searchText;
  }
};

SearchPanel.prototype.showResultList = function (start) {
  this.showSearchResults(this.resultList, start);
  this.updatePaging(start, this.resultList.length, this.showResultList.bind(this));
};

SearchPanel.prototype.setResultList = function (source, ids) {
  var showResultsFrom = this.showResultList.bind(this);
  this.source = source;
  this.resultList = [];
  for (var i = 0; i < ids.length; i++) {
    this.resultList.push({ id: ids[i], name: ids[i] });
  }
  showResultsFrom(0);
};

SearchPanel.prototype.setResultMessage = function (message) {
  // Empty panel and display message
  if (this.resultsElem) {
    this.resultsElem.empty();
    this.resultsElem.append(
      $('<span style="font-size:18pt;">' + message + '</span>')
    );
  }
};

SearchPanel.prototype.updatePreviewImages = function (previewImageIndex) {
  if (!this.resultsElem) return;
  var SearchPanel = this;
  var searchResults = this.resultsElem.find('.searchResult');
  this.previewImageIndex = previewImageIndex;
  searchResults.each(function (index, elem) {
    var resId = $(this).attr('id');
    var idfields = resId.split('_');
    var source = idfields[1];
    var id = idfields[2];
    var img = $(this).find('img.resultImg');
    var metadata = img.data('metadata');
    var imagePreviewUrl = SearchPanel.getImagePreviewUrlCallback(source, id, previewImageIndex, metadata);
    var fields = ['src', 'data-src'];
    for (var i = 0; i < fields.length; i++) {
      var field = fields[i];
      var src = img.attr(field);
      if (src && src.trim()) {
        img.attr(field, imagePreviewUrl);
      }
    }
    img.data('src', imagePreviewUrl);
  });
};

SearchPanel.prototype.unexpandGroup = function(event) {
  if (this.__grouped) {
    this.__searchSucceededWithDocs(this.__grouped, this.__groupedOptions);
  }
};

SearchPanel.prototype.createSearchResult = function (result, index, dimstr, isLazy) {
  var title = result.name;
  var group = result._group_;
  if (!title && group) {
    title = group.name;
  }
  if (this.tooltipIncludeAll) {
    title = JSON.stringify(result, this.tooltipIncludeFields, ' ');
  }
  var source = (result.source) ? result.source : this.source;
  var elem = $('<div></div>')
    .attr('class', 'searchResult')
    .attr('id', 'model_' + source + '_' + result.id)
    .attr('title', title)
    .data('result', result);
  var imgElem;
  if (this.getImagePreviewUrlCallback) {
    var imagePreviewUrl = this.getImagePreviewUrlCallback(source, result.id, this.previewImageIndex, result);
    imgElem = $('<img' + dimstr + '/>');
    //imgElem.addClass('freezeframe');
    if (isLazy) {
      //                imgElem.attr('src', '../resources/images/loading.gif');
      imgElem.attr('src', '');
      imgElem.attr('class', 'resultImg lazy');
      imgElem.attr('data-src', imagePreviewUrl);
    } else {
      imgElem.attr('src', imagePreviewUrl);
    }
    imgElem.attr('alt', result.id);
    imgElem.data('src', imagePreviewUrl);
    imgElem.data('metadata', result);
    elem.append(imgElem);
  } else {
    elem.append('<span>' + result.name + '</span>');
  }
  if (group && this.expandGroupCallback) {
//      var expandElem = $('<span class="glyphicon glyphicon-option-horizontal"></span>');
    var expandElem = $('<div class="glyphicon glyphicon-chevron-right"></div>');
    elem.append(expandElem);
    expandElem.click(function(event) {
      this.expandGroupCallback(group, index);
      this.Publish('GroupClicked', group, index);
      // So event doesn't propagate up to parent (elem)
      return false;
    }.bind(this));
    expandElem.mousedown(function (event) {
      // Prevent propagation up to parent (elem)
      event.preventDefault();
      return false;
    });
  }
  if (Constants.assetSources.model.indexOf(source) >= 0 || Constants.assetSources.scan.indexOf(source) >= 0) {
    // Models?
    if (!result.hasModel) {
      elem.addClass('searchResultNoModel');
    }
  }
  if (this.onClickResultCallback) {
    elem.click(function () {
      this.unselectOnPage(this.curSelectedIndex);
      elem.addClass('searchResultClicked');
      this.curSelectedIndex = index;
      this.onClickResultCallback(source, result.id, result, elem, index);
      this.Publish('ResultClicked', result, index);
    }.bind(this));
  }
  if (this.onHoverResultCallback) {
    elem.hover(function () {
      this.onHoverResultCallback(source, result.id, result, elem, index);
    }.bind(this));
  } else if (imgElem && this.showAnimatedOnHover) {
    var animatedGIFUrl = this.getImagePreviewUrlCallback(source, result.id, Constants.AssetGroup.ROTATING_IMAGE_INDEX, result);
    if (animatedGIFUrl) {
      elem.hover(
        function () { imgElem.attr('src', animatedGIFUrl); },
        function () { imgElem.attr('src', imgElem.data('src')); }
      );
    }
  }
  // Allows for dragging in
  if (this.onMousedownResultCallback) {
    elem.mousedown(function (event) {
      event.preventDefault();
      this.onMousedownResultCallback(source, result.id, result, elem, index);
    }.bind(this));
    this.Publish('ResultMousedown', result, index);
  }
  if (this.appendResultElemCallback) {
    this.appendResultElemCallback(source, result.id, result, elem);
  }
  return elem;
};

SearchPanel.prototype.__setResultsElemSize = function () {
  var searchOptionsElemHeight = (this.searchOptionsElem) ? this.searchOptionsElem.outerHeight() : 0;
  var restHeight = searchOptionsElemHeight + this.searchInfoElem.outerHeight();
  var resultsElemHeight = Math.max(10, this.container.height() - restHeight - 25);
  this.resultsElem.css('height', resultsElemHeight + 'px');
  //this.resultsElem.css('position', 'absolute');
  //this.resultsElem.css('top', restHeight + 'px');
};

SearchPanel.prototype.__getEntryWidth = function () {
  return Math.max(10, (this.resultsElem.width() - this.padding) / this.entriesPerRow - 10);
};

SearchPanel.prototype.showSearchResults = function (resultList, start) {
  if (!this.resultsElem) return;  // No place to put results!!!
  if (!start || start < 0) start = 0;

  this.__setResultsElemSize();
  this.resultsElem.empty();
  this.Publish('ClearResults');
  // If there were no search results, notify the user of this
  if (resultList.length === 0) {
    this.resultsElem.append('<span>No Results</span>');
    return;
  }

  var w = this.__getEntryWidth();
  var dimstr = ' width=' + w + ' ';
  var table = $('<table margin=2></table>');
  var row;
  var limit = Math.min(resultList.length, start + this.limit);
  for (var i = start; i < limit; i++) {
    var result = resultList[i];
    if (result) {
      var elem = this.createSearchResult(resultList[i], this.curStart + i, dimstr, this.loadImagesLazy);
      var tdElem = $('<td></td>').append(elem);
      if ((i % this.entriesPerRow) === 0) {
        row = $('<tr></tr>');
        table.append(row);
      }
      row.append(tdElem);
    } else {
      console.log('No result for index ' + i);
    }
  }
  this.resultsElem.append(table);
  if (this.loadImagesLazy) {
    table.find('img.lazy').lazy({
      bind: 'event',
      threshold: 50,
      visibleOnly: true,
      parent: this.resultsElem,
      appendScroll: this.resultsElem
    });
  } else {
  }
};

SearchPanel.prototype.setSearchCallback = function (callback) {
  this.searchTextElem.unbind('keyup');
  this.searchButton.unbind('click');
  this.searchTextElem.keyup(function(event) {
    event.stopPropagation();
    if (event.keyCode === 13) {
      callback(event);
    }
  }.bind(this));
  this.searchButton.click(callback);
};

SearchPanel.prototype.onResize = function (options) {
  if (!this.resultsElem) return;  // No place to put results!!!

  this.__setResultsElemSize();
  if (options && options.adjustEntriesPerRow && options.minEntryWidth && options.totalWidth) {
    var minWidth = Math.max(10, options.minEntryWidth + 10);
    this.entriesPerRow = Math.max(1, (options.totalWidth - this.padding) / minWidth);
    if (options.minEntriesPerRow && this.entriesPerRow < options.minEntriesPerRow) this.entriesPerRow = this.minEntriesPerRow;
    if (options.maxEntriesPerRow && this.entriesPerRow > options.maxEntriesPerRow) this.entriesPerRow = this.maxEntriesPerRow;
  }

  var imageElems = this.resultsElem.find('.searchResult > img');
  var w = this.__getEntryWidth();
  imageElems.each(function (index, elem) { $(this).attr('width', w); });

  if (this.loadImagesLazy) {
    this.resultsElem.find('img.lazy').lazy({
      bind: 'event',
      threshold: 50,
      visibleOnly: true,
      parent: this.resultsElem,
      appendScroll: this.resultsElem
    });
  }
  if (this.totalResults && this.pageElem) {
    var paginationDiv = this.pageElem.find('div.pagination');
    if (paginationDiv.length) {
      if (this.pageElem.width() !== paginationDiv.width()) {
        // TODO: This logic don't seem right, the callback may change!!!
        var showResultsFrom = this.showMoreSearchResults.bind(this);
        this.updatePaging(this.curStart, this.totalResults, showResultsFrom);
      }
    }
  }
};

SearchPanel.prototype.getTotalResults = function () {
  return this.totalResults;
};

SearchPanel.prototype.unselectOnPage = function (index) {
  if (index >= this.curStart && index < this.curEnd) {
    var i = index - this.curStart;
    var result = this.resultList[i];
    var source = (result.source) ? result.source : this.source;
    var elem = $('#' + 'model_' + source + '_' + result.id);
    elem.removeClass('searchResultClicked');
  }
};

SearchPanel.prototype.selectOnPage = function (index) {
  if (index >= this.curStart && index < this.curEnd) {
    // On this current page - select
    var i = index - this.curStart;
    var result = this.resultList[i];
    var source = (result.source) ? result.source : this.source;
    var elem = $('#' + 'model_' + source + '_' + result.id);
    if (elem.length > 0) {
      elem.click();
    } else if (this.onClickResultCallback) {
      this.unselectOnPage(this.curSelectedIndex);
      this.curSelectedIndex = index;
      this.onClickResultCallback(source, result.id, result, elem, index);
      this.Publish('ResultClicked', result, index);
    }
  }
};

SearchPanel.prototype.selectResult = function (index) {
  if (index >= 0 && index < this.totalResults) {
    if (index >= this.curStart && index < this.curEnd) {
      this.selectOnPage(index);
    } else {
      // Find correct start
      var start = this.curStart;
      if (index < start) {
        while (index < start) {
          start = Math.max(start - this.limit, 0);
        }
      } else {
        while (index >= start + this.limit) {
          start = start + this.limit;
        }
      }
      this.selectIndexOnSearchSucceeded = index;
      this.showMoreSearchResults(start);
    }
  }
};

SearchPanel.prototype.selectNext = function () {
  if (this.curSelectedIndex + 1 < this.totalResults) {
    this.selectResult(this.curSelectedIndex + 1);
  }
};

SearchPanel.prototype.selectPrev = function () {
  if (this.curSelectedIndex > 0) {
    this.selectResult(this.curSelectedIndex - 1);
  }
};

SearchPanel.prototype.__handlePaginationClick = function (callback, page, paginationContainer) {
  callback(page * this.limit);
  return false;
};

SearchPanel.prototype.updatePaging = function (start, numFound, callback) {
  this.totalResults = numFound;
  this.curStart = start;
  this.curEnd = Math.min(this.totalResults, this.limit + this.curStart);

  if (!this.pageElem) return;  // No page element to update
  this.pageElem.empty();
  //console.log("numFound=" + numFound + ", start=" + start);
  var last = Math.min(start + this.limit, numFound);
  var message = 'Displaying ' + (start + 1) + ' to ' + last + ' of ' + numFound;
  var textDiv = $('<div></div>').text(message);
  this.pageElem.append(textDiv);
  if (this.totalResults > this.limit) {
    var paginationElem = $('<div></div>').attr('class', 'pagination');
    paginationElem.width(this.pageElem.width());
    paginationElem.height('40px');
    var nDisplayEntries = Math.ceil(paginationElem.width() / 60);
    var curPage = Math.floor(this.curStart / this.limit);
    paginationElem.pagination(this.totalResults, {
      items_per_page: this.limit,
      current_page: curPage,
      load_first_page: false,
      num_edge_entries: 1,
      num_display_entries: nDisplayEntries,
      prev_text: '<',
      next_text: '>',
      //prev_show_always: false,
      //next_show_always: false,
      callback: this.__handlePaginationClick.bind(this, callback)
    });
    this.pageElem.append(paginationElem);
  } else {
    textDiv.css('height','40px');
  }
};

SearchPanel.prototype.showMoreSearchResults = function (start) {
  // Need more results...
  if (this.lastQuery) {
    var callback = this.lastQuery.searchDisplayOptions? this.lastQuery.searchDisplayOptions.callback : undefined;
    this.search(_.defaults(Object.create(null), { start: start }, this.lastQuery.query), callback, this.lastQuery.searchDisplayOptions);
  }
};

/**
 * Search has succeeded.  Display search results.
 * @param options {ui.SearchPanel.SearchDisplayOptions} Options on how the search results are to be displayed
 * @param data {{response: {docs: []}}|{grouped: {}} Object with list of documents or grouped documents
 */
SearchPanel.prototype.searchSucceeded = function (options, data) {
  this.Publish('SearchSucceededPreparePanel');
  if (data.response && data.response.docs) {
    this.__searchSucceededWithDocs(data.response, options);
  } else if (data.grouped) {
    // Pick the first of the grouped and display it
    // This handles group with group.query (not group with group.field)
    var docs = [];
    for (var p in data.grouped) {
      if (data.grouped.hasOwnProperty(p)) {
        var group = data.grouped[p];
        if (group.doclist && group.doclist.docs && group.doclist.docs.length > 0) {
          var doc = group.doclist.docs[0];
          doc['_group_'] = {
            name: p,
            numFound: group.doclist.numFound,
            fields: SolrQueryParser.simpleParse(p)
          };
          docs.push(doc);
        }
      }
    }
    if (options.ensureGroupedModelsUnique) {
      //console.log('filtering to ensure models are unique');
      //console.log('before filter');
      //console.log(docs);
      // filter the docs so scope if there are two with the same modelId
      var docsByModelId = _.groupBy(docs, 'id');
      for (var modelId in docsByModelId) {
        if (docsByModelId.hasOwnProperty(modelId)) {
          var mdocs = docsByModelId[modelId];
          if (mdocs.length > 0) {
            var repdoc = mdocs[0];
            for (var i = 1; i < mdocs.length; i++) {
              var doc = mdocs[i];
              if (doc['_group_'].numFound > repdoc['_group_'].numFound) {
                repdoc = doc;
              }
            }
            repdoc['_group_'].keep = true;
          }
        }
      }
      docs = docs.filter( function(x) { return x['_group_'].keep; });
      //console.log('after filter');
      //console.log(docs);
    }
    var doclist = {
      start: 0,
      numFound: docs.length,
      docs: docs
    };
    //console.log(doclist);
    this.__grouped = doclist;
    this.__groupedOptions = options;
    this.__searchSucceededWithDocs(doclist, options);
  } else {
    console.error('Cannot handle search response');
    console.log(data);
  }
};

SearchPanel.prototype.__searchSucceededWithDocs = function (doclist, options) {
  var resultList = doclist.docs;
  var ordering = (options) ? options.ordering : undefined;
  var callback = (options) ? options.callback : undefined;
  if (this.__unexpandGroupElem) {
    if (options.isExpandedGroup) {
      this.__unexpandGroupElem.show();
    } else {
      this.__unexpandGroupElem.hide();
    }
  }
  if (ordering !== undefined && ordering !== null) {
    // Assuming ordering is a array of { "modelId", "score" }
    // Convert ordering into a ordering lookup: a map with rank and score
    var orderingLookup = {};
    var idField = options.idField || 'fullId';
    for (var i = 0; i < ordering.length; i++) {
      var item = ordering[i];
      item['rank'] = i;
      orderingLookup[item[idField]] = item;
    }
    // Reorder resultList...
    var hasOrder = resultList.filter(function (element, index, array) { return orderingLookup[element.fullId]; });
    var rest = resultList.filter(function (element, index, array) { return !orderingLookup[element.fullId]; });
    var sorted = [];
    hasOrder.forEach(function (element, index, array) {
        sorted[orderingLookup[element.fullId]['rank']] = element;
        element['score'] = orderingLookup[element.fullId]['score'];
      }
    );
    // Filter out items scope weren't found
    sorted = sorted.filter(function (element, index, array) { return element; });
    resultList = sorted.concat(rest);
  }

  // Save resultList...
  this.resultList = resultList;
  if (!doclist.start && doclist.numFound === doclist.docs.length) {
    // Everything is already here
    //console.log("Show results with all results");
    this.updatePaging(doclist.start, doclist.numFound);  // update paging element before display of result list so things are laid out nicely
    this.showResultList(0);
  } else {
    var showResultsFrom = this.showMoreSearchResults.bind(this);
    this.updatePaging(doclist.start, doclist.numFound, showResultsFrom);
    this.showSearchResults(resultList);
  }
  if (typeof (callback) === 'function') {
    callback(this.source, resultList);
  } else if (this.selectIndexOnSearchSucceeded !== undefined) {
    this.selectOnPage(this.selectIndexOnSearchSucceeded);
    this.selectIndexOnSearchSucceeded = undefined;
  } else if (this.searchSucceededCallback) {
    this.searchSucceededCallback(this.source, resultList);
  } else {
    // console.log("Search succeeded!");
  }
};

SearchPanel.prototype.searchFailed = function (err) {
  if (this.searchFailedCallback) {
    this.searchFailedCallback(err);
  } else {
    this.showSearchFailedMessage(err);
  }
};

SearchPanel.prototype.showSearchFailedMessage = function (message) {
  if (this.resultsElem) {
    this.resultsElem.empty();
    this.Publish('ClearResults');
    this.resultsElem.append(
      $('<span>There was an error when processing the search request</span>' +
        '<br/>' +
        '<span>Error: ' + message + '</span>')
    );
  } else {
    console.error('Search failed: ' + message);
  }
};

SearchPanel.prototype.getQuerySortOrder = function () {
  if (this.sortOrder) {
    // already has sortOrder
    // Append additional sort orders to it
    return this.sortOrder + ',' + this.additionalSortOrder;
  } else {
    // No sort order
    return this.additionalSortOrder;
  }
};

SearchPanel.prototype.startSearch = function (callback) {
  this.sortOrder = this.sortOrderElem ? this.sortOrderElem.val() : this.sortOrder || '';
  var searchTerm = this.getSearchText();
  this.Publish('startSearch', searchTerm);
  this.search({ searchText: searchTerm, start: 0}, callback);
};

SearchPanel.prototype.refreshSearch = function () {
  this.sortOrder = this.sortOrderElem ? this.sortOrderElem.val() : this.sortOrder || '';
  var searchTerm = this.getSearchText();
  this.search({ searchText: searchTerm, start: this.curStart});
};

/**
 * Performs search!
 * @param query {string|{}} Search text or query options
 * @param [callback] {function(source,results)}
 * @param [searchDisplayOptions] {ui.SearchPanel.SearchDisplayOptions} Additional options for search display
 */
SearchPanel.prototype.search = function (query, callback, searchDisplayOptions) {
  if (!query) return;
  var queryOpts = query;
  if (_.isString(query)) {
    queryOpts = { searchText: query };
  }
  var sortOrder = this.getQuerySortOrder();
  queryOpts = _.defaults(Object.create(null), queryOpts, { sort: sortOrder, source: this.source, start: 0, limit: this.limit });

  var updatedSearchDisplayOptions = _.defaults(Object.create(null), { callback: callback }, searchDisplayOptions || {});
  this.lastQuery = { query: queryOpts, searchDisplayOptions: searchDisplayOptions };

  var scope = this;
  console.log('queryOpts', queryOpts);
  this.searchModule.query(queryOpts, function(err, res) {
    if (err) {
      scope.searchFailed(err);
    } else {
      scope.searchSucceeded(updatedSearchDisplayOptions, res);
    }
  });
};

SearchPanel.prototype.selectSource = function (source) {
  this.source = source;
  if (this.sourceElem) { this.sourceElem.val(this.source); }
};

SearchPanel.prototype.hasSource = function (source) {
  return this.sources.indexOf(source) >= 0;
};

SearchPanel.prototype.registerSource = function(source) {
  if (this.sources.indexOf(source) < 0) {
    this.sources.push(source);
    if (this.sourceElem) {
      this.sourceElem.append('<option value="' + source + '">' + source + '</option>');
    }
  }
};

SearchPanel.prototype.getLastQueryURL = function (opts) {
  if (this.lastQuery) {
    return this.searchModule.getQueryUrl(_.defaults(Object.create(null), opts, this.lastQuery.query, { start: 0, limit: 1000000}));
  }
};

SearchPanel.prototype.searchCategories = function (categories, groupCategories, callback) {
  var searchOptions = {};
  if (categories && Array.isArray(categories) && categories.length) {
    searchOptions.limit = categories.length;
    if (groupCategories) {
      searchOptions['group'] = true;
      searchOptions['group.query'] = categories.map(function (x) { return 'category:' + x; });
      //customOptions['group.limit'] = 2; // number of entries to return per group
    }
  }

  searchOptions.searchText = this.searchModule.getQuery('category', categories);
  this.search(searchOptions, callback, {isExpandedGroup: false, ensureGroupedModelsUnique: true });
};

SearchPanel.prototype.expandCategory = function (category, callback) {
  var queryString = this.searchModule.getQuery('category', category);
  this.search({ searchText: queryString }, callback, { isExpandedGroup: true });
};

/**
 * Event indicating search panel is being cleared
 * @event ui.SearchPanel#ClearResults
 */

/**
 * Event indicating group is being clicked
 * @event ui.SearchPanel#GroupClicked
 * @param group
 * @param index
 */

/**
 * Event indicating result is being clicked
 * @event ui.SearchPanel#ResultClicked
 * @param result
 * @param index
 */

/**
 * Event indicating mousedown happening on result
 * @event ui.SearchPanel#ResultMousedown
 * @param result
 * @param index
 */

/**
 * Event indicating search started
 * @event ui.SearchPanel#startSearch
 * @param searchTerm
 */

/**
 * Event indicating search has succeeded and search panel is being prepared to show results
 * @event ui.SearchPanel#SearchSucceededPreparePanel
 */

/**
 * @typedef {Object} {SearchPanel.SearchDisplayOptions}
 * @memberOf ui
 * @property [ordering] {Array<{id:string}>} Explicit list of ordered assets
 * @property {idField} {string} What field to use for id
 * @property [isExpandedGroup=false] {boolean} Whether the group is expanded
 * @property [ensureGroupedModelsUnique=false] {boolean} Filter grouped models so they are unique
 */
module.exports = SearchPanel;
