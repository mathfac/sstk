// Create a grid that can be used for navigation
var BBox = require('geo/BBox');
var Constants = require('Constants');
var Colors = require('util/Colors');
var MeshHelpers = require('geo/MeshHelpers');
var Object3DUtil = require('geo/Object3DUtil');
var Graph = require('nav/Graph');
var PathFinder = require('nav/PathFinder');
var PubSub = require('PubSub');
var Timings = require('util/Timings');
var _ = require('util');

Graph.createFromJson = function(json, opts) {
    opts = _.clone(opts);
    opts.json = json;
    if (json.type === 'SceneGrid2D') {
        return new SceneGrid2D(opts);
    } else if (json.type === 'MultiLevelGrid') {
        return new MultiLevelGrid(opts);
    } else {
        console.warn('Invalid grid type: ' + json.type);
    }
};

Graph.ClearanceTileWeight = 1.0;
Graph.MaxTileWeight = 20000;

/**
 * 2D grid for a scene
 * @param opts
 * @param opts.cellSize {number} Size of a cell in virtual units
 * @param opts.bbox {geo.BBox} Bounding box of scene/level this grid will cover
 * @param opts.unit {number} Virtual unit (in m) used for the cellSize
 * @param opts.floorHeight {number} Floor height in virtual units
 * @param [opts.metadata] {Object} Additional metadata stored with the graph
 * @param [opts.useEdgeWeights=false] {boolean} Whether edge weights should be
 * used
 * @param [opts.allowDiagonalMoves=false] {boolean} Whether diagonal moves are
 * allowed
 * @param [opts.reverseEdgeOrder=false] {boolean} Whether edge order should be
 * reversed for nice stair casing
 * @extends nav.SquareGrid
 * @constructor
 * @memberOf nav
 */
function SceneGrid2D(opts) {
    var gridOpts = {
        metadata : opts.metadata,
        useEdgeWeights : opts.useEdgeWeights,
        dirs : opts.allowDiagonalMoves ? Graph.SquareGrid.DIRS8
                                       : Graph.SquareGrid.DIRS4,
        precomputeEdges : false,
        reverseEdgeOrder : opts.reverseEdgeOrder
    };
    // console.log('gridOpts', gridOpts);
    if (opts.json) {
        Graph.SquareGrid.call(this, 0, 0, gridOpts);
        // Construct from serialized json
        this.fromJson(opts.json, gridOpts);
    } else {
        this.cellSize = opts.cellSize;
        this.unit = opts.unit;

        var bbox = opts.bbox;
        var bbdims = bbox.dimensions();
        var centroid = bbox.centroid();
        var n1 = Math.ceil(bbdims.x / this.cellSize);
        var n2 = Math.ceil(bbdims.z / this.cellSize);

        var cellSize = this.cellSize;
        var size1 = n1 * cellSize;
        var size2 = n2 * cellSize;
        var min1 = centroid.x - size1 / 2;
        var min2 = centroid.z - size2 / 2;
        var minHeight = bbox.min.y;
        if (opts.floorHeight != undefined) {
            minHeight = opts.floorHeight;
        }
        this.min = [ min1, minHeight, min2 ];
        this.bbox = bbox; // Kept just as extra info
        Graph.SquareGrid.call(this, n1, n2, gridOpts);
    }
}

SceneGrid2D.prototype = Object.create(Graph.SquareGrid.prototype);
SceneGrid2D.constructor = SceneGrid2D;

Object.defineProperty(SceneGrid2D.prototype, 'floorHeight',
                      {get : function() { return this.min[1]; }});

SceneGrid2D.prototype.positionToCell = function(position, dir = null) {
    var i = Math.floor((position.x - this.min[0]) / this.cellSize);
    var j = Math.floor((position.z - this.min[2]) / this.cellSize);
    var isValid = (i >= 0 && i <= this.width) && (j >= 0 && j <= this.height);
    if (isValid && this.bbox) {
        isValid =
            position.y >= this.bbox.min.y && position.y <= this.bbox.max.y;
    }
    return {
        i : i,
        j : j,
        dir : dir,
        id : this.toId(i, j, dir),
        isValid : isValid
    };
};
SceneGrid2D.prototype.idToCell = function(id) {
    var ij = this.fromId(id);
    return {i : ij[0], j : ij[1], dir : ij[2], id : id};
};
SceneGrid2D.prototype.idToPosition = function(id) {
    var ij = this.fromId(id);
    var cellHeight = this.getCellAttribute('floorHeight', id);
    if (cellHeight == undefined || !isFinite(cellHeight)) {
        cellHeight = this.min[1];
    }
    return [
        this.min[0] + (ij[0] + 0.5) * this.cellSize, cellHeight,
        this.min[2] + (ij[1] + 0.5) * this.cellSize
    ];
};

// Represent this grid as a serializable json object
SceneGrid2D.prototype.toJson = function() {
    var json = Graph.SquareGrid.prototype.toJson.call(this);
    json['type'] = 'SceneGrid2D';
    json['cellSize'] = this.cellSize;
    json['min'] = this.min;
    json['unit'] = this.unit;
    json['bbox'] = this.bbox.toJSON();
    json['timings'] = this.timings;
    return json;
};
SceneGrid2D.prototype.fromJson = function(json, opts) {
    Graph.SquareGrid.prototype.fromJson.call(this, json, opts);
    this.cellSize = json.cellSize;
    this.min = json.min;
    this.timings = json.timings;
    this.unit = json.unit || 1; // Assume unit is 1 for saved out grids
    if (json.bbox) {
        this.bbox = new BBox(Object3DUtil.toVector3(json.bbox.min),
                             Object3DUtil.toVector3(json.bbox.max));
    }
};
SceneGrid2D.prototype.toVirtualUnits = function() {
    if (this.unit !== Constants.virtualUnitToMeters) {
        var scaleBy = this.unit / Constants.virtualUnitToMeters;
        this.min = _.map(this.min, function(x) { return x * scaleBy; });
        if (this.bbox) {
            this.bbox.min.multiplyScalar(scaleBy);
            this.bbox.max.multiplyScalar(scaleBy);
        }
        if (this._tileAttributes && this._tileAttributes['floorHeight']) {
            var floorHeights = this._tileAttributes['floorHeight'];
            for (var i = 0; i < floorHeights.length; i++) {
                floorHeights[i] = floorHeights[i] * scaleBy;
            }
        }
        this.cellSize = this.cellSize * scaleBy;
        this.unit = Constants.virtualUnitToMeters;
    }
};

/**
 * Refines grid based on options
 * @param opts.radius {number} Radius in virtual units
 * @param opts.clearance {number} Clearance in virtual units
 * @param opts.populateOccupancy {boolean}
 * @param opts.estimateFloorHeight {boolean}
 */
SceneGrid2D.prototype.refine = function(opts, levels) {
    // Refine grid based on options
    // Stuff current weights into occupancy attribute
    opts = opts || {};
    opts = _.defaults(opts, {radius : 0, clearance : 0, minFloorFrac : 0.7});
    var occupancy = this._tileAttributes['occupancy'];
    var floorHeights = this._tileAttributes['floorHeight'];
    var roomIndex = this._tileAttributes['roomIndex'];
    if (opts.adjustRoomIndex) {
        // TODO: Remove this logic
        console.log('Adjust stored room index');
        for (var i = 0; i < floorHeights.length; i++) {
            if (_.isFinite(floorHeights[i])) {
                roomIndex[i] += 1;
            }
        }
    }
    if (opts.populateOccupancy) {
        console.log('Populating occupancy');
        for (var i = 0; i < this._weights.length; i++) {
            occupancy[i] = _.isFinite(this._weights[i]) ? 0 : 1;
        }
    }
    if (opts.estimateFloorHeight) {
        console.log('estimate floor height');
        var position = new THREE.Vector3();
        var min1 = this.min[0];
        var min2 = this.min[2];
        for (var i = 0; i < this._weights.length; i++) {
            if (!_.isFinite(floorHeights[i])) {
                var cell = this.idToCell(i);
                position.set(min1 + (cell.i + 0.5) * this.cellSize, this.min[1],
                             min2 + (cell.j + 0.5) * this.cellSize);
                var data = opts.getCellData(position, opts.level || levels[0]);
                if (_.isFinite(data.floorHeight)) {
                    floorHeights[i] = data.floorHeight;
                    occupancy[i] = data.occupancy;
                    this._weights[i] = data.occupancy ? Infinity : 1;
                }
            }
        }
    }
    // Clearance and radius in number of tiles
    if (opts.updateWeights) {
        console.log('Updating Weights');
        var r = (opts.radius + 0.05) / this.cellSize;
        var r2 = r * r;
        var cr = (opts.radius + opts.clearance) / this.cellSize;
        var cr2 = cr * cr;
        var d = cr ? (Math.ceil(cr - 0.5)) : 0;
        // console.log('Got cell r and d', r, d);
        for (var i = 0; i < this._weights.length; i++) {
            // TODO: check if more than 70% of cells have valid floorHeights[i]
            if (d > 0) {
                var cell = this.idToCell(i);
                var x1 = Math.max(cell.i - d, 0);
                var y1 = Math.max(cell.j - d, 0);
                var x2 = Math.min(cell.i + d, this.width);
                var y2 = Math.min(cell.j + d, this.height);
                var nHasFloor = 0;
                var nCells = 0;
                this._weights[i] = 1;
                for (var x = x1; x < x2; x++) {
                    var dx = cell.i - x;
                    for (var y = y1; y < y2; y++) {
                        var dy = cell.j - y;
                        var dxy2 = (dx * dx + dy * dy);
                        if (dxy2 < r2) {
                            var id = this.toId(x, y);
                            nCells++;
                            if (_.isFinite(floorHeights[id])) {
                                nHasFloor++;
                            }
                            if (occupancy[id]) {
                                this._weights[i] = Infinity;
                                x = x2;
                                break;
                            }
                        } else if (dxy2 < cr2) {
                            var id = this.toId(x, y);
                            if (occupancy[id]) {
                                this._weights[i] = Graph.ClearanceTileWeight;
                            }
                        }
                    }
                }
                if (nHasFloor / nCells < opts.minFloorFrac) {
                    this._weights[i] = Infinity;
                }
            } else {
                this._weights[i] = occupancy[i] ? Infinity : 1;
            }
        }
    }

    // this.make_grid_finer();
};

SceneGrid2D.prototype.getLinePath = function(cell1, cell2) {
    var cellIds = [];
    var dx = cell2.i - cell1.i;
    var dy = cell2.j - cell1.j;
    var adx = Math.abs(dx);
    var ady = Math.abs(dy);
    if (adx === 0 && ady === 0) {
        cellIds.push(cell1.id);
    } else if (adx >= ady) {
        var c1 = cell1;
        var c2 = cell2;
        var swapped = false;
        if (dx < 0) {
            c1 = cell2;
            c2 = cell1;
            dx = -dx;
            dy = -dy;
            swapped = true;
        }
        var m = dy / dx;
        for (var i = 0; i <= dx; i++) {
            var y = Math.round(c1.j + i * m);
            var cid = this.toId(c1.i + i, y);
            cellIds.push(cid);
        }
        if (swapped) {
            cellIds.reverse();
        }
    } else {
        var c1 = cell1;
        var c2 = cell2;
        var swapped = false;
        if (dy < 0) {
            c1 = cell2;
            c2 = cell1;
            dx = -dx;
            dy = -dy;
        }
        var m = dx / dy;
        for (var j = 0; j <= dy; j++) {
            var x = Math.round(c1.i + j * m);
            var cid = this.toId(x, c1.j + j);
            cellIds.push(cid);
            swapped = true;
        }
        if (swapped) {
            cellIds.reverse();
        }
    }
    // if (cellIds.length > 1) {
    //   console.log('got path', cellIds, cell1, cell2);
    // }
    return cellIds;
};

// TODO: Enhance MultiLevelGrid so it has everything a graph needs to have
/**
 * Multi level grid for a scene
 * @param opts
 * @param opts.cellSize {number} Size of a cell in virtual units
 * @param opts.unit {number} Virtual unit (in m) used for the cellSize
 * @param [opts.metadata] {Object} Additional metadata stored with the graph
 * @constructor
 * @memberOf nav
 */
function MultiLevelGrid(opts) {
    PubSub.call(this, opts);
    opts = opts || {};
    // Multiple level grid
    this.metadata = opts.metadata;
    this.numNodes = 0;
    this.levelGrids = [];
    this.activeLevel = opts.activeLevel;
    if (opts.json) {
        // Construct from serialized json
        this.fromJson(opts.json, opts);
    } else {
        this.cellSize = opts.cellSize;
        this.unit = opts.unit;
    }
}

MultiLevelGrid.prototype = Object.create(PubSub.prototype);
MultiLevelGrid.constructor = MultiLevelGrid;

MultiLevelGrid.prototype.getLevel = function(i) { return this.levelGrids[i]; };
MultiLevelGrid.prototype.addLevel = function(levelGrid) {
    // NOTE: numNodes and nodeIdOffset need to be adjusted if individual level
    // grid changes
    levelGrid.nodeIdOffset = this.numNodes;
    levelGrid.level = this.levelGrids.length;
    this.levelGrids.push(levelGrid);
    this.numNodes += levelGrid.numNodes;
    // Forward all events from levelGrid
    var scope = this;
    levelGrid.SubscribeAll(this, function(event) {
        if (event.startsWith('tile')) {
            // Event for tile (update argument[1] which corresponds to the id
            arguments[1] = arguments[1] + levelGrid.nodeIdOffset;
        }
        scope.Publish.apply(scope, arguments);
    });
};

MultiLevelGrid.prototype.clear = function() {
    for (var i = 0; i < this.levelGrids.length; i++) {
        this.levelGrids[i].Unsubscribe(PubSub.ALL, this);
    }
    this.numNodes = 0;
    this.levelGrids = [];
    this.activeLevel = undefined;
};

MultiLevelGrid.prototype.fromJson = function(json, opts) {
    this.clear();
    var levelGrids = _.map(json.grids, function(grid) {
        return Graph.createFromJson(grid, opts);
    });
    for (var i = 0; i < levelGrids.length; i++) {
        this.addLevel(levelGrids[i]);
    }
    this.metadata = json.metadata;
    this.cellSize = json.cellSize;
    this.timings = json.timings;
    this.unit = json.unit;
};
MultiLevelGrid.prototype.toEncodedPixels = function(key, encodeFn) {
    var levels = _.map(this.levelGrids, function(levelGrid) {
        var pixels = levelGrid.toEncodedPixels(key, encodeFn);
        return pixels;
    });
    return levels;
};
MultiLevelGrid.prototype.toPixels = function(key) {
    var levels = _.map(this.levelGrids, function(levelGrid) {
        var pixels = levelGrid.toPixels(key);
        return pixels;
    });
    return levels;
};
MultiLevelGrid.prototype.toJson = function() {
    var gridsJson = _.map(this.levelGrids, function(levelGrid) {
        var gridJson = levelGrid.toJson();
        return gridJson;
    });
    return {
        type: 'MultiLevelGrid', metadata: this.metadata,
            cellSize: this.cellSize, unit: this.unit, grids: gridsJson,
            timings: this.timings
    }
};

MultiLevelGrid.prototype.idToLevelGrid = function(id) {
    for (var i = 0; i < this.levelGrids.length; i++) {
        var levelGrid = this.levelGrids[i];
        if (id >= levelGrid.nodeIdOffset &&
            id < (levelGrid.nodeIdOffset + levelGrid.numNodes)) {
            return levelGrid;
        }
    }
};
MultiLevelGrid.prototype.positionToCell = function(position) {
    // Assume there isn't that many levels, get cell for each and then figure
    // out which is best
    for (var i = 0; i < this.levelGrids.length; i++) {
        var levelGrid = this.levelGrids[i];
        var cell = levelGrid.positionToCell(position);
        cell.level = levelGrid.level;
        cell.id = cell.id + levelGrid.nodeIdOffset;
        // TODO: Better selection of cells
        if (cell.isValid) {
            return cell;
        }
    }
};
MultiLevelGrid.prototype.idToCell = function(id) {
    var levelGrid = this.idToLevelGrid(id);
    if (levelGrid) {
        var cell = levelGrid.idToCell(id - levelGrid.nodeIdOffset);
        cell.level = levelGrid.level;
        return cell;
    }
};
MultiLevelGrid.prototype.idToPosition = function(id) {
    var levelGrid = this.idToLevelGrid(id);
    return levelGrid ? levelGrid.idToPosition(id - levelGrid.nodeIdOffset)
                     : undefined;
};
MultiLevelGrid.prototype.toVirtualUnits = function() {
    for (var i = 0; i < this.levelGrids.length; i++) {
        var levelGrid = this.levelGrids[i];
        levelGrid.toVirtualUnits();
    }
};
MultiLevelGrid.prototype.tileWeight = function(id) {
    var levelGrid = this.idToLevelGrid(id);
    return levelGrid ? levelGrid.tileWeight(id - levelGrid.nodeIdOffset)
                     : Infinity;
};
MultiLevelGrid.prototype.setTileWeight = function(id, w) {
    var levelGrid = this.idToLevelGrid(id);
    if (levelGrid) {
        levelGrid.setTileWeight(id - levelGrid.nodeIdOffset, w)
    }
};
MultiLevelGrid.prototype.createCellAttributes = function(key, cellAttr, force) {
    for (var i = 0; i < this.levelGrids.length; i++) {
        var levelGrid = this.levelGrids[i];
        levelGrid.createCellAttributes(key, cellAttr, force);
    }
};
MultiLevelGrid.prototype.setCellAttribute = function(id, key, value) {
    var levelGrid = this.idToLevelGrid(id);
    if (levelGrid) {
        levelGrid.setCellAttribute(id - levelGrid.nodeIdOffset, key, value);
    }
};
MultiLevelGrid.prototype.getCellAttribute = function(id, key) {
    var levelGrid = this.idToLevelGrid(id);
    return levelGrid
               ? levelGrid.getCellAttribute(id - levelGrid.nodeIdOffset, key)
               : undefined;
};
MultiLevelGrid.prototype.getCellAttributeStatistics = function(key) {
    var stats;
    for (var i = 0; i < this.levelGrids.length; i++) {
        var levelGrid = this.levelGrids[i];
        var levelStats = levelGrid.getCellAttributeStatistics(key);
        // merge stats
        if (stats) {
            stats.min = Math.min(levelStats.min, stats.min);
            stats.max = Math.max(levelStats.max, stats.max);
            stats.sum = levelStats.sum + stats.sum;
            stats.size = levelStats.size + stats.size;
        } else {
            stats = levelStats;
        }
    }
    return stats;
};
MultiLevelGrid.prototype.setUserData = function(id, data) {
    var levelGrid = this.idToLevelGrid(id);
    if (levelGrid) {
        levelGrid.setUserData(id, data);
    }
};
MultiLevelGrid.prototype.getUserData = function(id) {
    var levelGrid = this.idToLevelGrid(id);
    return levelGrid ? levelGrid.getUserData(id - levelGrid.nodeIdOffset)
                     : undefined;
};
MultiLevelGrid.prototype.edgeWeight = function(id1, id2) {
    var levelGrid1 = this.idToLevelGrid(id1);
    var levelGrid2 = this.idToLevelGrid(id2);
    if (levelGrid1 !== levelGrid2) {
        return Infinity;
    } else {
        var levelGrid = levelGrid1;
        return levelGrid.edgeWeight(id1 - levelGrid.nodeIdOffset,
                                    id2 - levelGrid.nodeIdOffset);
    }
};
MultiLevelGrid.prototype.distance = function(id1, id2) {
    var levelGrid1 = this.idToLevelGrid(id1);
    var levelGrid2 = this.idToLevelGrid(id2);
    if (levelGrid1 !== levelGrid2) {
        return Infinity;
    } else {
        var levelGrid = levelGrid1;
        return levelGrid.distance(id1 - levelGrid.nodeIdOffset,
                                  id2 - levelGrid.nodeIdOffset);
    }
};
// Is there an edge from id1 to id2?
MultiLevelGrid.prototype.hasEdge = function(id1, id2) {
    var levelGrid1 = this.idToLevelGrid(id1);
    var levelGrid2 = this.idToLevelGrid(id2);
    if (levelGrid1 !== levelGrid2) {
        return Infinity;
    } else {
        var levelGrid = levelGrid1;
        return levelGrid.hasEdge(id1 - levelGrid.nodeIdOffset,
                                 id2 - levelGrid.nodeIdOffset);
    }
};
// Return index of edge between id1 and id2 (-1 if no edge)
MultiLevelGrid.prototype.edgeIndex = function(id1, id2) {
    var levelGrid1 = this.idToLevelGrid(id1);
    var levelGrid2 = this.idToLevelGrid(id2);
    if (levelGrid1 !== levelGrid2) {
        return Infinity;
    } else {
        var levelGrid = levelGrid1;
        return levelGrid.edgeIndex(id1 - levelGrid.nodeIdOffset,
                                   id2 - levelGrid.nodeIdOffset);
    }
};
// All edges from id
MultiLevelGrid.prototype.edgesFrom = function(id1) {
    var levelGrid = this.idToLevelGrid(id1);
    return levelGrid ? levelGrid.edgesFrom(id1 - levelGrid.nodeIdOffset) : [];
};
// All edges as a list of [id1, id2]
MultiLevelGrid.prototype.allEdges = function() {
    var all = [];
    for (var i = 0; i < this.levelGrids.length; i++) {
        var levelGrid = this.levelGrids[i];
        all.push.apply(all, levelGrid.allEdges());
    }
    return all;
};
MultiLevelGrid.prototype.fromId = function(id) {
    var c = this.idToCell(id);
    return [ c.i, c.j, c.level ];
};
MultiLevelGrid.prototype.getActive = function() {
    if (this.activeLevel != undefined) {
        return this.levelGrids[this.activeLevel];
    } else {
        // Full grid active!  So exciting!
        return this;
    }
};
MultiLevelGrid.prototype.refine = function(opts, levels) {
    // Refine grid based on options
    opts = opts || {};
    for (var i = 0; i < this.levelGrids.length; i++) {
        var level = levels[i];
        var levelOpts = _.defaults({level : level}, opts);
        var levelGrid = this.levelGrids[i];
        levelGrid.refine(levelOpts);
    }
};
MultiLevelGrid.prototype.getLinePath = function(cell1, cell2) {
    if (cell1.level === cell2.level) {
        var levelGrid = this.levelGrids[cell1.level];
        return levelGrid.getLinePath(cell1, cell2);
    }
};

/**
 * Navigation structure for a scene - includes support for visualization, graph
 * for path planning, and precomputed paths
 * @param opts Configuration
 * @param opts.sceneState {scene.SceneState} Scene state to build/load
 * navigation mesh for
 * @param opts.cellSize {number} Size of a navigation cell (in virtual units).
 * @param opts.getCellData {function(THREE.Vector3): {cost: number, roomIndex:
 * int, floorHeight: number}} Function returning the cell data at the given
 * positions
 * @param [opts.cellAttributes] {Object<string, {type: string, compute:
 * function(celldata): number}>} Array of dense attributes to keep in cells
 * @param [opts.isValid] {function} Function returning if the position is valid
 * (used only for visualization now)
 * @param [opts.useEdgeWeights=true] {boolean} Whether to use edge weights or
 * tile cost
 * @param [opts.allowDiagonalMoves=true] {boolean} Whether to allow movement to
 * an diagonal tile
 * @param [opts.reverseEdgeOrder=true] {boolean} Whether to reverse edge order
 * for stair casing
 * @param [opts.baseTileHeight=0] {number} Base tile height (for visualization
 * purposes)
 * @param [opts.tileOverlap=0] {number} How much should tiles overlap (for
 * visualization purposes).
 *   Should be between -1 to 1.  0 = no overlap, <0 will have a bit of a gap, >0
 * will have overlap.
 * @param [opts.tileOpacity=0.5] {number} Tile opacity (for visualization
 * purposes)
 * @param [opts.autoCreateGrid=true] {boolean} Whether to automatically create
 * the navgrid if not available from cache
 *   (this may take several minutes)
 * @param [opts.autoUpdate=true] {boolean} Whether to automatically update the
 * position of the agent and the shortest path
 * @param [opts.refineGrid] {{radius: number}} Options for how to refine loaded
 * grid
 * @param [opts.mapName=navmap] {string} Default name of navigation map to use
 * @constructor
 * @memberOf nav
 */
function NavScene(opts) {
    this.__timings = new Timings();
    this.metadata = opts.metadata; // Additional metadata to store away (not
                                   // used explicitly)
    this.sceneState = opts.sceneState;
    this.getCellData = opts.getCellData; // Returns the metadata for cell
                                         // position (e.g. cost at the given
                                         // positions }
    this.cellAttributes = opts.cellAttributes;
    this.isValid = opts.isValid; // Returns if the position is valid (used only
                                 // for visualization now)
    this.baseTileHeight = opts.baseTileHeight ||
                          0; // Base tile height (for visualization purposes)
    this.tileOverlap =
        opts.tileOverlap ||
        0; // How much should tiles overlap (for visualization purposes)
    this.tileOpacity =
        opts.tileOpacity || 0.5; // Tile opacity (for visualization purposes)
    this.cellSize = opts.cellSize;
    this.autoUpdate = (opts.autoUpdate !== undefined) ? opts.autoUpdate : true;
    this.autoCreateGrid =
        (opts.autoCreateGrid !== undefined) ? opts.autoCreateGrid : true;
    this.useEdgeWeights =
        (opts.useEdgeWeights !== undefined) ? opts.useEdgeWeights : true;
    this.allowDiagonalMoves = (opts.allowDiagonalMoves !== undefined)
                                  ? opts.allowDiagonalMoves
                                  : true;
    this.nrots = (opts.nrots !== undefined) ? opts.nrots : 40;
    this.refinment_factor =
        (opts.refinment_factor !== undefined) ? opts.refinment_factor : 6.0;

    this.reverseEdgeOrder =
        (opts.reverseEdgeOrder !== undefined) ? opts.reverseEdgeOrder : true;
    this.refineGrid = opts.refineGrid;
    this.mapName =
        opts.mapName || 'navmap'; // Name of field to use for navmap in asset
    this.pathFinder = new PathFinder();
    this.grid = null;
    this.mapModes = [ 'none', 'tileWeight', 'pathCost', 'path' ];
    if (this.cellAttributes) {
        var scope = this;
        _.each(this.cellAttributes,
               function(cellAttr, key) { scope.mapModes.push(key); });
    }
    this.__tmpPosition = new THREE.Vector3();
    this.__initGrid();
    this.clear();
}

NavScene.CellAttribute = Graph.CellAttribute;

// Initializes the grid
NavScene.prototype.__initGrid = function() {
    var gridInfo = this.sceneState.info[this.mapName];
    if (gridInfo && gridInfo.data) {
        // Want to load grid from data
        console.log('Parsing grid for ' + this.sceneState.info.fullId + ' ' +
                    this.mapName);
        var grid = this.__parseGrid(gridInfo.data);
        grid.toVirtualUnits();
        // Use grid cellSize and height...
        this.cellSize = grid.cellSize;
        var hasOccupancy = false;
        if (grid.metadata) {
            if (grid.metadata.cellAttributes) {
                hasOccupancy = grid.metadata.cellAttributes['occupancy'];
                this.cellAttributes = this.cellAttributes || {};
                var cellAttributes = this.cellAttributes;
                _.each(grid.metadata.cellAttributes, function(json) {
                    if (!cellAttributes[json.name]) {
                        cellAttributes[json.name] = new NavScene.CellAttribute(
                            json.name, json.type, json);
                    }
                });
            }
        }
        // TOOO: Consider keeping a copy of old grid
        if (this.refineGrid) {
            console.log('Refining grid');
            _.each(cellAttributes, function(cellAttr, key) {
                grid.createCellAttributes(key, cellAttr);
            });
            var refineOpts = _.defaults({}, this.refineGrid, {
                updateWeights : true,
                populateOccupancy : !hasOccupancy,
                getCellData : this.getCellData
            });
            grid.metadata.refineGrid = refineOpts;
            grid.refine(refineOpts, this.sceneState.getLevels());
        }
        this.grid = grid;
    } else if (this.autoCreateGrid) {
        // Create the grid
        console.log('Creating grid for ' + this.sceneState.info.fullId + ' ' +
                    this.mapName);
        this.grid = this.__createGrid();
        if (this.refineGrid) {
            console.log('Refining grid');
            var refineOpts = _.defaults(
                {}, this.refineGrid,
                {updateWeights : true, getCellData : this.getCellData});
            this.grid.metadata.refineGrid = refineOpts;
            this.grid.refine(this.refineGrid, this.sceneState.getLevels());
        }
    } else {
        // Don't want to create the grid
        console.warn('No precomputed grid for ' + this.sceneState.getFullID() +
                     ': aborting navscene creation');
        this.grid = null;
    }
    this.baseGrid = this.grid;
};

NavScene.prototype.setActiveLevel = function(level) {
    console.log('set active level', level);
    if (_.isString(level)) {
        level = parseInt(level);
    }
    this.baseGrid.activeLevel = level;
    this.grid = this.baseGrid.getActive();
    // console.log('got active grid', this.grid);
    this.clear();
};

NavScene.prototype.isGridValid = function() { return !!this.grid; };

/**
 * Clears the map (of precomputed costs to goals), visualization, and paths
 */
NavScene.prototype.clear = function() {
    this.clearVisualization();
    this.map = [];
    this.start = null;
    this.goals = [];
    this.__tiles = {};
    this.__shortestPath = null;
    this.__lastSearchState = null;
    this.__currentPosition = null;
    this.__currentOrientation = null;
};

/**
 * Initializes the navigation structure with agent, new start and goals
 * @param agent {sim.Agent}
 * @param start
 * @param goals
 * @param mapState {nav.MapState} State of precomputed map for shortest path
 */
NavScene.prototype.reset = function(agent, start, goals, mapState) {
    this.clear();
    this.start = start;
    this.goals = goals;

    var reuseMap = false;
    // TODO: Make mapState a proper class
    if (mapState) {
        this.map = mapState.map;
        this.__lastSearchState = mapState.__lastSearchState;
        reuseMap = true;
    }
    // Compute shortest path with cell information saved inside start and goals
    this.__currentPosition = start.position.clone();
    this.update(agent);
};

NavScene.prototype.update = function(agent) {
    var position = agent.position;
    this.__currentOrientation = agent.getOrientation();
    this.__shortestPath = this.computeShortestPath({
        position : position,
        dir : [ -this.__currentOrientation.x, -this.__currentOrientation.z ],
        maxTileWeight : Graph.MaxTileWeight,
        prevPosition : this.__currentPosition
    },
                                                   this.goals,
                                                   {reuseMap : true});
    this.__currentPosition.copy(position);
    if (this.__shortestPath && this.__shortestPath.next) {
        this.__shortestPath.direction =
            this.__getDirection(this.grid, agent, this.__shortestPath.next.id);
    }
};

NavScene.prototype.__parseGrid = function(json) {
    var opts = {
        useEdgeWeights : this.useEdgeWeights,
        allowDiagonalMoves : this.allowDiagonalMoves,
        reverseEdgeOrder : this.reverseEdgeOrder,
        nrots : this.nrots,
        refinment_factor : this.refinment_factor
    };
    return Graph.createFromJson(json, opts);
};

NavScene.prototype.__createGrid2D = function(opts) {
    this.__timings.start('createGrid2D');
    var grid = new SceneGrid2D({
        bbox : opts.bbox,
        cellSize : this.cellSize,
        floorHeight : opts.floorHeight,
        unit : Constants.virtualUnitToMeters,
        metadata : opts.metadata,
        useEdgeWeights : this.useEdgeWeights,
        allowDiagonalMoves : this.allowDiagonalMoves,
        reverseEdgeOrder : this.reverseEdgeOrder
    });
    var doors = this.sceneState.findModelInstances(function(mi) {
        return mi.model.isDoor();
    });
    var position = new THREE.Vector3();
    var tileBBox = new BBox();
    var cellSize = this.cellSize;
    var min1 = grid.min[0];
    var min2 = grid.min[2];
    var dy = 0.02 * Constants.metersToVirtualUnit;
    console.log('navScene gridSize', grid.width, grid.height);
    var cellAttributes = this.cellAttributes; // Per cell attributes
    if (cellAttributes) {
        console.log('navScene cellAttributes', _.keys(cellAttributes));
        _.each(cellAttributes, function(cellAttr, key) {
            grid.createCellAttributes(
                key, cellAttr); // allocate memory for cell attributes
        });
    }
    var floorHeight = grid.floorHeight;
    for (var i = 0; i < grid.width; i++) {
        for (var j = 0; j < grid.height; j++) {
            var id = grid.toId(i, j);
            position.set(min1 + (i + 0.5) * cellSize, floorHeight,
                         min2 + (j + 0.5) * cellSize);
            var celldata = this.getCellData(position, opts.level);
            celldata.position = position;
            grid.setTileWeight(id, celldata.cost);
            if (cellAttributes) {
                _.each(cellAttributes, function(cellAttr, key) {
                    var value = cellAttr.compute(celldata);
                    grid.setCellAttribute(id, key, value); // set cell attribute
                });
            }
            // check if this tile intersect a door
            tileBBox.min.set(min1 + i * cellSize, floorHeight,
                             min2 + j * cellSize);
            tileBBox.max.set(min1 + (i + 1) * cellSize, floorHeight + dy,
                             min2 + (j + 1) * cellSize);
            var intersected = _.filter(doors, function(mi) {
                var bbox = Object3DUtil.getBoundingBox(mi.object3D);
                return bbox.intersects(tileBBox);
            });
            intersected = _.map(intersected,
                                function(x) { return x.object3D.userData.id; });
            if (intersected.length > 0) {
                grid.setUserData(id, {doors : intersected});
            }
        }
    }
    this.__timings.stop('createGrid2D');
    grid.timings = {createMillis : this.__timings.getDuration('createGrid2D')};
    return grid;
};

NavScene.prototype.__createGrid = function() {
    var levels = this.sceneState.getLevels();
    var metadata = _.defaults({
        sceneId : this.sceneState.getFullID(),
        cellAttributes : _.map(this.cellAttributes,
                               function(cellAttr, key) {
                                   return {
                                       name: cellAttr.name,
                                           type: cellAttr.type.name,
                                           dataType: cellAttr.dataType
                                   }
                               })
    },
                              this.metadata || {});
    if (levels && levels.length > 1) {
        console.log('creating MultiLevelGrid with ' + levels.length);
        this.__timings.start('createMultiLevelGrid');
        var multilevelGrid = new MultiLevelGrid({
            metadata : metadata,
            cellSize : this.cellSize,
            unit : Constants.virtualUnitToMeters
        });
        for (var i = 0; i < levels.length; i++) {
            var level = levels[i];
            var levelGrid = this.__createGrid2D({
                bbox : Object3DUtil.getBoundingBox(level),
                floorHeight : this.sceneState.getFloorHeight(level),
                level : level
            });
            multilevelGrid.addLevel(levelGrid);
        }
        this.__timings.stop('createMultiLevelGrid');
        multilevelGrid.timings = {
            createMillis : this.__timings.getDuration('createMultiLevelGrid')
        };
        return multilevelGrid;
    } else {
        var defaultFloorHeight = this.sceneState.getFloorHeight();
        return this.__createGrid2D({
            metadata : metadata,
            bbox : this.sceneState.getBBox(),
            floorHeight : this.sceneState.getFloorHeight(this.sceneState.scene,
                                                         defaultFloorHeight),
            level : this.sceneState.scene
        });
    }
};

NavScene.prototype.getCellId = function(p, saveCell) {
    var cell = this.grid.positionToCell(p.position, p.dir);
    if (saveCell) {
        p.cell = cell;
    }
    return cell.id;
};

NavScene.prototype.__getLinePathCellIds = function(start, end) {
    var startCell = this.grid.positionToCell(start);
    var endCell = this.grid.positionToCell(end);
    var cellIds = this.grid.getLinePath(startCell, endCell);
    return cellIds;
};

/**
 * Make sure all cells between start and end are unoccupied
 * @param start {THREE.Vector3} Start position (world coordinates)
 * @param end {THREE.Vector3} End position (world coordinates)
 * @param radius
 * @param filter
 */
NavScene.prototype.checkLinePathUnoccupied = function(start, end, radius,
                                                      filter) {
    // console.log('checkLinePathUnoccupied', start, end);
    // console.log('checkLinePath', startCell, endCell, cellIds);
    var cellIds = this.__getLinePathCellIds(start, end);
    if (!cellIds) {
        return false;
    }
    for (var i = 0; i < cellIds.length; i++) {
        var cellId = cellIds[i];
        var tileWeight = this.grid.tileWeight(cellId);
        // console.log('got tileWeight', cellId, tileWeight);
        // TODO: use filter function instead of hardcoded weight comparison
        var okay = !this.grid.getCellAttribute(cellId, 'occupancy');
        if (!okay) {
            return false;
        }
    }
    return true;
};

/**
 * Make sure position is unoccupied
 * @param pos
 * @param radius
 * @param filter
 */
NavScene.prototype.checkPositionUnoccupied = function(pos, radius, filter) {
    // TODO: take radius into account by checking all cells in radius
    var cell = this.grid.positionToCell(pos);
    if (cell) {
        var cellId = cell.id;
        radius /= this.grid.cellSize;
        var r2 = radius * radius;
        for (var dy = -radius; dy <= radius; dy++) {
            for (var dx = -radius; dx <= radius; dx++) {
                if (dy * dy + dx + dx > r2)
                    continue;

                var ij = this.grid.fromId(cellId);
                ij[0] += dx;
                ij[1] += dy;

                var testId = this.grid.toId(ij[0], ij[1]);
                if (this.grid.getCellAttribute(testId, 'occupancy')) {
                    return false;
                }
            }
        }
        return !this.grid.getCellAttribute(cellId, 'occupancy');
    } else {
        return false;
    }

};

NavScene.prototype.getCellAttribute = function(pos, key) {
    var cellId = this.getCellId({position : pos});
    var value = this.grid.getCellAttribute(cellId, key);
    return value;
};

NavScene.prototype.getRoom = function(pos) {
    var roomIndex = this.getCellAttribute(pos, 'roomIndex');
    return this.sceneState.getRoomByIndex1(roomIndex);
};

NavScene.prototype.__getDirection = (function() {
    var nextPos = new THREE.Vector3();
    var localPos = new THREE.Vector3();
    return function(graph, agent, cellId) {
        if (agent) {
            var pos2 = graph.idToPosition(cellId);
            nextPos.set(pos2[0],
                        pos2[1] +
                            agent.originHeight * Constants.metersToVirtualUnit,
                        pos2[2]);
            agent.worldToLocalPositionNoScaling(nextPos, localPos);
            localPos.normalize();
            return localPos.toArray();
        }
    };
})();

NavScene.prototype.__getPathDoors = function(path) {
    // Return doors we encountered in this path
    var graph = this.grid;
    var doors = [];
    var lastDoors = [];
    for (var i = 0; i < path.length; i++) {
        var data = graph.getUserData(path[i]);
        if (data && data.doors) {
            for (var j = 0; j < data.doors.length; j++) {
                if (lastDoors.indexOf(data.doors[j]) < 0) {
                    doors.push(data.doors[j]);
                }
            }
            lastDoors = data.doors;
        } else {
            lastDoors = [];
        }
    }
    return doors;
};

NavScene.prototype.__getPathRoomIndices = function(path) {
    // Return rooms we encountered in this path
    var graph = this.grid;
    var roomIndices = [];
    for (var i = 0; i < path.length; i++) {
        var index = graph.getCellAttribute(path[i], 'roomIndex');
        if (index > 0) {
            if (roomIndices.indexOf(index) < 0) {
                roomIndices.push(index);
            }
        }
    }
    return roomIndices;
};

NavScene.prototype.__getPath = function(cellId) {
    var scope = this;
    var map = this.map;
    var path = this.pathFinder.getPath(map, cellId);
    var realPath =
        _.map(path, function(id) { return scope.grid.idToPosition(id); });

    var result = {path : path, realPath : realPath};
    var firstCell = map[cellId];
    if (firstCell) {
        result.start = firstCell;
        result.cost = firstCell.cost_so_far;
        result.distance = firstCell.distance * this.grid.cellSize;
        result.doors = this.__getPathDoors(path);
        result.rooms = this.__getPathRoomIndices(path);
    }
    if (firstCell && firstCell.parent) {
        var parentId = firstCell.parent;
        var parentCell = map[parentId];
        if (parentCell) {
            result.next = parentCell;
        }
    }
    result.isValid = !!firstCell && _.isFinite(result.cost);
    return result;
};

/**
 * Computes shorted path from current position to a set of goals
 * @param position
 * @param goals
 * @param opts
 * @returns shortest path
 */
NavScene.prototype.computeShortestPath = function(position, goals, opts) {
    opts = opts || {};
    // Figure out grid cell for start and grid cell for end
    // Convert from starts/goal to grid cells
    var scope = this;
    var mapState = opts.mapState ? opts.mapState : this;
    mapState.map = (mapState.map && opts.reuseMap) ? mapState.map : [];
    mapState.map = [];
    if (mapState.map.length === 0) {
        mapState.__lastSearchState = null;
    }
    var currCell = scope.getCellId(position, opts.saveCell);
    if (position.maxTileWeight != undefined) {
        // HACK! Make sure all tiles between previous position and position are
        // updated to be traversable!!!
        var cellIds;
        if (opts.prevPosition) {
            cellIds = this.__getLinePathCellIds(opts.prevPosition, position);
        } else {
            cellIds = [ currCell ];
        }
        for (var i = 0; i < cellIds.length; i++) {
            var cellId = cellIds[i];
            if (this.grid.tileWeight(cellId) > position.maxTileWeight) {
                this.grid.setTileWeight(cellId, position.maxTileWeight);
                mapState.__lastSearchState =
                    null; // Force recompute (slower but more correct)
                mapState.map = [];
            }
        }
    }

    if (!mapState.map[currCell]) {
        var goalCells = _.map(goals, function(x) {
            var cid = scope.getCellId(x, opts.saveCell);
            return cid;

            var radius = 1.0 / scope.grid.cellSize;
            var r2 = radius * radius;

            var res = [];
            for (var dz = -radius; dz <= radius; dz++) {
                for (var dx = -radius; dx <= radius; dx++) {
                    if ((dz * dz + dx * dx) > r2)
                        continue;
                    var xy = scope.grid.fromId(cid);
                    xy[0] += dx;
                    xy[1] += dz;
                    var new_cid = scope.grid.toId(xy[0], xy[1]);
                    if (scope.grid.isTraversable(xy[0], xy[1]) &&
                        (scope.grid.getCellAttribute(cid, 'roomIndex') ===
                         scope.grid.getCellAttribute(new_cid, 'roomIndex'))) {
                        res.push(new_cid);
                    }
                }
            }
            return _.chain(res)
                .sortBy(function(id) { return scope.grid.distance(cid, id); })
                .first()
                .value();
        });
        // console.time('NavScene.getShortestPath');
        // console.log('currCell and goals', currCell, goalCells);
        var heuristic = this.allowDiagonalMoves
                            ? PathFinder.heuristics.octile_grid2d
                            : PathFinder.heuristics.manhattan_grid2d;
        mapState.__lastSearchState = this.pathFinder.theta_star_search(
            goalCells, currCell, this.grid, mapState.map,
            mapState.__lastSearchState, heuristic);
        // console.timeEnd('NavScene.getShortestPath');
    }
    var path = this.__getPath(currCell);
    return path;
};

NavScene.prototype.getShortestPath = function() { return this.__shortestPath; };

NavScene.prototype.__checkPath = function(path, opts) {
    var okay = path && path.isValid;
    if (okay && opts) {
        if (opts.minSteps != undefined) {
            if (path.start.steps < opts.minSteps) {
                return false;
            }
        }
        if (opts.maxSteps != undefined) {
            if (path.start.steps > opts.maxSteps) {
                return false;
            }
        }
        if (opts.minDist != undefined) {
            if (path.distance < opts.minDist) {
                return false;
            }
        }
        if (opts.maxDist != undefined) {
            if (path.distance > opts.maxDist) {
                return false;
            }
        }
        if (opts.minCost != undefined) {
            if (path.cost < opts.minCost) {
                return false;
            }
        }
        if (opts.maxCost != undefined) {
            if (path.cost > opts.maxCost) {
                return false;
            }
        }
        if (false && opts.minRooms != undefined) {
            if ((path.rooms.length - 1) < opts.minRooms) {
                return false;
            }
        }
        if (opts.maxRooms != undefined) {
            if ((path.rooms.length - 1) > opts.maxRooms) {
                return false;
            }
        }
    }
    return okay;
};

NavScene.prototype.hasValidPath = function(opts) {
    if (this.__shortestPath) {
        this.__shortestPath.valid = this.__checkPath(this.__shortestPath, opts);
        return this.__shortestPath.valid;
    } else {
        return false;
    }
};

NavScene.defaultTileColors = {
    obstacle : new THREE.Color('gray'),
    traversable : new THREE.Color('blue'),
    // Default colors (used in top down rendering)
    'start' : new THREE.Color('purple'),
    'goal' : new THREE.Color('green'),
    'shortestPath' : new THREE.Color('blue'),
    'trace' : new THREE.Color('red'),
    'position' : new THREE.Color('darkred')
};

NavScene.defaultTileColorFn = function(w) {
    return (w === Infinity) ? NavScene.defaultTileColors.obstacle
                            : NavScene.defaultTileColors.traversable;
};

NavScene.prototype.__getTileMaterial = function(color, opacity) {
    if (!this.__materialCache) {
        this.__materialCache = {};
    }
    color = Object3DUtil.getColor(color);
    var key = color.getHexString() + ':' + opacity;
    if (this.__materialCache[key]) {
        return this.__materialCache[key];
    } else {
        var material = Object3DUtil.getBasicMaterial(color, opacity);
        this.__materialCache[key] = material;
        return material;
    }
};

NavScene.prototype.__createTiles = function(opts) {
    opts = opts || {};
    var baseHeight = this.baseTileHeight || 0;
    var cellSize = this.cellSize;
    var dy = 0.02 * Constants.metersToVirtualUnit;
    var d = cellSize * (0.5 + this.tileOverlap);
    var min = new THREE.Vector3();
    var max = new THREE.Vector3();
    var grid = this.grid;
    var tiles = [];
    var tilesGroup = new THREE.Group();
    tilesGroup.name = opts.name || 'NavSceneTiles';
    var getColor = opts.getColor || NavScene.defaultTileColorFn;
    var tmpVector = new THREE.Vector3();
    for (var i = 0; i < grid.numNodes; i++) {
        var position = grid.idToPosition(i);
        if (this.isValid) {
            tmpVector.set(position[0], position[1], position[2]);
            var isValid = this.isValid(tmpVector);
            if (!isValid)
                continue; // Skip non valid tiles
        }
        var tileWeight =
            opts.getWeight ? opts.getWeight(i) : grid.tileWeight(i);
        var tileHeight = baseHeight + position[1];
        min.set(position[0] - d, tileHeight, position[2] - d);
        max.set(position[0] + d, tileHeight + dy, position[2] + d);
        var color = getColor(tileWeight, i);
        var mat = this.__getTileMaterial(color, this.tileOpacity);
        // console.log('tile at ' + i, min, max, color);
        var tile = new MeshHelpers.BoxMinMax(min, max, mat);
        tile.name = 'Tile' + i;
        tilesGroup.add(tile);
        tiles[i] = tile;
    }
    console.log('tile bbox', Object3DUtil.getBoundingBox(tilesGroup));
    return {group : tilesGroup, tiles : tiles};
};

NavScene.prototype.__recolorTiles = function(opts) {
    opts = opts || {};
    var getColor = opts.getColor || NavScene.defaultTileColorFn;
    var grid = this.grid;
    var tiles = opts.tiles;
    for (var i = 0; i < grid.numNodes; i++) {
        var tile = tiles[i];
        if (tile) {
            var tileWeight =
                opts.getWeight ? opts.getWeight(i) : grid.tileWeight(i);
            var color = getColor(tileWeight, i);
            var newmat = this.__getTileMaterial(color, tile.material.opacity);
            if (newmat !== this.material) {
                tile.material = newmat;
            }
            // tile.material.color.set(color);
        }
    }
};

NavScene.prototype.clearVisualization = function() {
    if (this.__tiles) {
        var node = this.sceneState.debugNode;
        _.each(this.__tiles,
               function(tiles, key) { node.remove(tiles.group); });
    }
};

NavScene.prototype.__visualizeTiles = function(opts) {
    opts = _.defaults(opts, {type : 'nav'});
    var tiles = this.__tiles[opts.type];
    if (tiles) {
        opts.tiles = tiles.tiles;
        this.__recolorTiles(opts);
    } else {
        tiles = this.__createTiles(opts);
        this.__tiles[opts.type] = tiles;
    }
    this.sceneState.debugNode.add(tiles.group);
    // Make sure debugNode is visible
    Object3DUtil.setVisible(this.sceneState.debugNode, true);
};

NavScene.prototype.visualizeTileAttribute = function(attrname) {
    // console.log('visualizing navigation tiles');
    var grid = this.grid;
    var attrinfo = this.cellAttributes[attrname];
    var isCategorical = attrinfo.dataType === 'categorical';
    // TODO: get counts for categorical data
    var stats =
        isCategorical ? undefined : grid.getCellAttributeStatistics(attrname);
    this.__visualizeTiles({
        getWeight : function(id) {
            var attribute = grid.getCellAttribute(id, attrname);
            return attribute;
        },
        getColor :
            isCategorical
                ? function(w) { return Colors.createColor(w); }
                : Colors.getColorFunction({
                      type : 'interpolate',
                      space : 'hsl',
                      infinity : NavScene.defaultTileColors.obstacle,
                      colors : [
                          new THREE.Color('blue'), new THREE.Color('white'),
                          new THREE.Color('red')
                      ],
                      minWeight : (attrinfo.min != null) ? attrinfo.min
                                                         : stats.min,
                      maxWeight : (attrinfo.max != null) ? attrinfo.max
                                                         : stats.max
                  })
    });
};

NavScene.prototype.visualizeTileWeight = function() {
    // console.log('visualizing navigation tiles');
    this.__visualizeTiles({
        getColor : Colors.getColorFunction({
            type : 'interpolate',
            space : 'hsl',
            infinity : NavScene.defaultTileColors.obstacle,
            colors : [
                new THREE.Color('blue'), new THREE.Color('orange'),
                new THREE.Color('red')
            ],
            minWeight : 1,
            maxWeight : Graph.MaxTileWeight
        })
    });
};

NavScene.prototype.visualizeTraversable = function(color) {
    // console.log('visualizing navigation tiles');
    this.__visualizeTiles({
        getColor : function(w) {
            if (w === Infinity) {
                return NavScene.defaultTileColors.obstacle;
            } else {
                return color;
            }
        }
    });
};

NavScene.prototype.visualizePathCost = function(opts) {
    opts = opts || {};
    var showPathOnly = opts.showPathOnly;
    // console.log('visualizing navigation map');
    var colorFn = Colors.getColorFunction({
        type : 'interpolate',
        space : 'hsl',
        infinity : NavScene.defaultTileColors.obstacle,
        colors :
            opts.colors || [ new THREE.Color('blue'), new THREE.Color('red') ],
        maxWeight : Math.ceil(Math.sqrt(this.grid.numNodes))
    });
    var startId = this.start.cell ? this.start.cell.id : -1;
    var goalIds = this.goals.map(function(g) { return g ? g.cell.id : -1; });
    var pathIds = this.__shortestPath ? this.__shortestPath.path : [];
    var map = this.map;
    this.__visualizeTiles({
        name : 'NavScenePathCosts',
        type : 'map',
        getWeight : function(id) {
            var cell = map[id];
            return cell ? cell.steps : Infinity;
        },
        getColor : function(weight, id) {
            if (id === startId) {
                return Colors.lighten('yellow');
            } else if (goalIds.indexOf(id) >= 0) {
                return Colors.lighten('green');
            } else {
                var onpath = pathIds.indexOf(id) >= 0;
                if (showPathOnly) {
                    return onpath ? colorFn(weight)
                                  : NavScene.defaultTileColors.obstacle;
                } else {
                    var color = colorFn(weight);
                    return onpath ? Colors.lighten(color) : color;
                }
            }
        }
    });
};

/**
 * Show map as overlay in scene
 * @param mode {string}
 *  tileWeight - weight of tile
 *  pathCost - cost of tile as distance to goal
 *  path - shortest path
 */
NavScene.prototype.showMap = function(mode) {
    var mapMode = _.isString(mode) ? mode : this.mapModes[mode];
    // console.log('Show map for mode ' + mapMode);
    this.clearVisualization();
    if (mapMode === 'tileWeight') {
        this.visualizeTileWeight();
    } else if (mapMode === 'pathCost') {
        this.visualizePathCost({showPathOnly : false});
    } else if (mapMode === 'path') {
        this.visualizePathCost({showPathOnly : true});
    } else if (this.cellAttributes) {
        var cellAttr = this.cellAttributes[mapMode];
        if (cellAttr) {
            this.visualizeTileAttribute(mapMode);
        }
    }
};

NavScene.prototype.__setPixelColors = function(pixelLayers, cellIds, color) {
    for (var i = 0; i < cellIds.length; i++) {
        var cellId = cellIds[i];
        if (cellId >= 0) {
            var c = this.grid.idToCell(cellId);
            if (Array.isArray(pixelLayers)) {
                var layer = pixelLayers[c.level];
                if (layer) {
                    if (c.i >= 0 && c.i < layer.width && c.j >= 0 &&
                        c.j <= layer.height) {
                        var pi = 4 * (c.i + c.j * layer.width);
                        layer.data[pi] = Math.floor(color.r * 255);
                        layer.data[pi + 1] = Math.floor(color.g * 255);
                        layer.data[pi + 2] = Math.floor(color.b * 255);
                        layer.data[pi + 3] = 255;
                    }
                }
            } else {
                var layer = pixelLayers;
                if (c.i >= 0 && c.i < layer.width && c.j >= 0 &&
                    c.j <= layer.height) {
                    var pi = 4 * (c.i + c.j * layer.width);
                    layer.data[pi] = Math.floor(color.r * 255);
                    layer.data[pi + 1] = Math.floor(color.g * 255);
                    layer.data[pi + 2] = Math.floor(color.b * 255);
                    layer.data[pi + 3] = 255;
                }
            }
        }
    }
};

/**
 * Fetch map as image
 * @param opts
 */
NavScene.prototype.__getEncodedMap = function(opts) {
    // console.log('get encoded map', opts);
    // Then we have some overlays ('shortestPath', 'trace', 'start', 'goal',
    // 'position')
    var pixelType = opts.baseMap;
    var encodingFn = this.getEncodeDataFn(pixelType); // Update
    var pixels = this.grid.toEncodedPixels(
        (pixelType === 'tileWeight' || pixelType === 'tileTraversibility')
            ? null
            : pixelType,
        encodingFn);
    // Overlays
    if (opts.overlays.indexOf('shortestPath') >= 0) {
        // Add shortest path
        var pathIds = this.__shortestPath ? this.__shortestPath.path : [];
        this.__setPixelColors(pixels, pathIds, opts.colors['shortestPath']);
    }
    if (opts.overlays.indexOf('trace') >= 0) {
        // Add agent trace of positions
        // this.__setPixelColors(pixels, pathIds, opts.colors['trace']);
    }
    if (opts.overlays.indexOf('start') >= 0) {
        // Add start
        var startId = this.start.cell ? this.start.cell.id : -1;
        this.__setPixelColors(pixels, [ startId ], opts.colors['start']);
    }
    if (opts.overlays.indexOf('goal') >= 0) {
        // Add goals
        var goalIds =
            this.goals.map(function(g) { return g ? g.cell.id : -1; });
        this.__setPixelColors(pixels, goalIds, opts.colors['goal']);
    }
    if (opts.overlays.indexOf('position') >= 0) {
        // Add current position
        if (this.__currentPosition) {
            // var currCell = this.grid.positionToCell(this.__currentPosition);
            // this.__setPixelColors(pixels, [currCell.id],
            // opts.colors['position']);
            this.__tmpPosition.copy(this.__currentPosition)
                .addScaledVector(this.__currentOrientation, 5 * this.cellSize);
            var cellIds = this.__getLinePathCellIds(this.__currentPosition,
                                                    this.__tmpPosition);
            this.__setPixelColors(pixels, cellIds, opts.colors['position']);
        }
    }
    return pixels;
};

NavScene.prototype.getEncodedMap = function(opts) {
    opts = _.defaults(Object.create(null), opts || {}, {
        baseMap : 'tileTraversibility',
        overlays : [ 'shortestPath', 'trace', 'start', 'goal', 'position' ],
        colors : NavScene.defaultTileColors
    });
    return this.__getEncodedMap(opts);
};

NavScene.prototype.getEncodeDataFn = function(attrname) {
    if (attrname === 'tileWeight') {
        return function(d) {
            return _.isFinite(d) ? [ 0, 0, d, 255 ] : [ 0, 0, 0, 0 ];
        };
    } else if (attrname === 'tileTraversibility') {
        var stats = this.grid.getCellAttributeStatistics('tileWeight');
        var cf = Colors.getColorFunction({
            type : 'interpolate',
            space : 'hsl',
            infinity : NavScene.defaultTileColors.obstacle,
            colors : [ new THREE.Color('white'), new THREE.Color('yellow') ],
            minWeight : stats.min,
            maxWeight : stats.max
        });
        return function(d) {
            if (_.isFinite(d)) {
                var c = cf(d);
                return [
                    Math.floor(c.r * 255), Math.floor(c.g * 255),
                    Math.floor(c.b * 255), 255
                ];
            } else {
                return [ 127, 127, 127, 255 ];
            }
        };
    } else {
        var grid = this.grid;
        var attrinfo = this.cellAttributes[attrname];
        var isCategorical = attrinfo.dataType === 'categorical';
        if (isCategorical) {
            return function(d) {
                var c = Colors.createColor(d, Constants.defaultPalette);
                if (_.isFinite(d)) {
                    return [
                        Math.floor(c.r * 255), Math.floor(c.g * 255),
                        Math.floor(c.b * 255), 255
                    ];
                } else {
                    return [ 0, 0, 0, 0 ];
                }
            };
        } else {
            var stats = grid.getCellAttributeStatistics(attrname);
            var cf = Colors.getColorFunction({
                type : 'interpolate',
                space : 'hsl',
                infinity : NavScene.defaultTileColors.obstacle,
                colors : [ new THREE.Color('yellow'), new THREE.Color('red') ],
                minWeight : (attrinfo.min != null) ? attrinfo.min : stats.min,
                maxWeight : (attrinfo.max != null) ? attrinfo.max : stats.max
            });
            return function(d) {
                if (_.isFinite(d)) {
                    var c = cf(d);
                    return [
                        Math.floor(c.r * 255), Math.floor(c.g * 255),
                        Math.floor(c.b * 255), 255
                    ];
                } else {
                    return [ 0, 0, 0, 0 ];
                }
            };
        }
    }
};

NavScene.prototype.__pathScorer = function(opts) {
    var scope = this;
    var mapState = opts.mapState || {map : []};
    var position = new THREE.Vector3();
    return function(i) {
        var ij = scope.grid.fromId(i);
        if (scope.grid.isTraversable(ij[0], ij[1])) {
            var p = scope.grid.idToPosition(i);
            position.set(p[0], p[1], p[2]);
            var path = scope.computeShortestPath(
                {position : position}, opts.goals,
                {mapState : mapState, reuseMap : true, saveCell : true});
            var okay = scope.__checkPath(path, opts);
            return okay ? 1 : 0;
        } else {
            return 0;
        }
    }
};

NavScene.prototype.__simpleScorer = function(opts) {
    var scope = this;
    return function(i) {
        var cost = scope.grid.tileWeight(i);
        if (_.isFinite(cost)) {
            return 1;
        } else {
            return 0;
        }
    };
};

function convertSample(grid, s) {
    // sampleWeight = s.weight;
    if (s.weight > 0) {
        var c = grid.idToCell(s.value);
        c.position = Object3DUtil.toVector3(grid.idToPosition(s.value));
        c.roomIndex = grid.getCellAttribute(s.value, 'roomIndex');
        return c;
    }
}

NavScene.prototype.sample = function(opts) {
    var sampler = opts.sampler;
    var sampled = sampler.sample({
        elements : _.getIterator(0, this.grid.numNodes, 1),
        scorer : opts.scorer === 'path' ? this.__pathScorer(opts)
                                        : this.__simpleScorer(opts),
        nsamples : opts.nsamples
    });
    var grid = this.grid;
    if (_.isArray(sampled)) {
        sampled = _.filter(_.map(sampled, function(s) {
            return convertSample(grid, s);
        }), function(s) { return s; });
    } else if (sampled) {
        sampled = convertSample(grid, sampled);
    }
    // console.log(sampled);
    return sampled;
};

module.exports = NavScene;
