// Adapted from grid.js at http://www.redblobgames.com/pathfinding/
// Copyright 2014 Red Blob Games
// License: Apache v2
var PubSub = require('PubSub');
var TypeUtils = require('data/TypeUtils');
var _ = require('util');

function CellAttribute(name, type, opts) {
    opts = opts || {};
    if (_.isString(type)) {
        type = TypeUtils.nameToTypedArray(type) || type;
    }
    this.name = name;
    this.type = type || Float32Array;
    this.dataType = opts.dataType;
    this.compute =
        opts.compute || function(cellData) { return cellData[name]; };
}

/**
 * Graph class
 * Weights are assigned to *nodes* not *edges*, but the pathfinder
 * will need edge weights, so we treat an edge A->B as having the
 * weight of tile B. If a tile weight is Infinity we don't expose
 * edges to it.
 * @param numNodes {int} Number of nodes in the graph
 * @param [opts.useEdgeWeights=false] {boolean} Whether to use edge weights or
 * not
 * @param [opts.precomputeEdges] {boolean} Whether edges are precomputed
 * @param [opts.metadata] {Object} Additional metadata stored with the graph
 * @constructor
 * @memberOf nav
 */
function Graph(numNodes, opts) {
    PubSub.call(this);
    opts = opts || {};
    this.metadata = opts.metadata; // Other metadata we want to store
    this.useEdgeWeights = opts.useEdgeWeights;
    this.numNodes = numNodes;
    this._edges = [];       // node id to list of node ids
    this._edgeWeights = []; // node id to edge weight for that node
    this._weights = [];     // node id to number (could be Infinity)
    this._tileAttributes =
        {}; // dense tile attributes (map of attribute name to dense vector)
    this._userData = {}; // sparse per tile user data (arbitrary objects)
    for (var id = 0; id < numNodes; id++) {
        this._weights[id] = 1;
        if (opts.precomputeEdges !== false) {
            this._edges[id] = [];
            if (this.useEdgeWeights) {
                this._edgeWeights[id] = [];
            }
        }
    }
}

Graph.prototype = Object.create(PubSub.prototype);
Graph.prototype.constructor = Graph;

Graph.prototype.createCellAttributes = function(key, cellAttr, force) {
    if (force || !this._tileAttributes[key]) {
        var type = cellAttr.type || Float32Array;
        this._tileAttributes[key] = new type(this.numNodes);
    }
};
Graph.prototype.setCellAttribute = function(id, key, value) {
    this._tileAttributes[key][id] = value;
};
Graph.prototype.getCellAttribute = function(id, key) {
    return this._tileAttributes[key] ? this._tileAttributes[key][id]
                                     : undefined;
};
Graph.prototype.getCellAttributeStatistics = function(key) {
    var values = this._tileAttributes[key];
    if (key === 'tileWeight') {
        values = this._weights;
    }
    var min = Infinity, max = -Infinity, sum = 0, n = 0;
    for (var i = 0; i < values.length; i++) {
        var v = values[i];
        if (_.isFinite(v)) {
            min = Math.min(v, min);
            max = Math.max(v, max);
            sum += v;
            n++;
        }
    }
    return {min : min, max : max, sum : sum, size : n};
};
// Weights are given to tiles, not edges, but the search interface
// will only ask about edges. Weight of edge id1->id2 is of tile id2.
Graph.prototype.getOrderedWeights = function() {
    var weightCounts = [];
    this._weights.forEach(function(w) {
        if (weightCounts[w])
            weightCounts[w] = weightCounts[w] + 1;
        else
            weightCounts[w] = 1;
    });
    var orderedWeights =
        Object.keys(this._weights).map(function(x) { return parseInt(x); });
    return orderedWeights;
};

Graph.prototype.setUserData = function(id, data) {
    var oldData = this._userData[id];
    this._userData[id] = data;
    this.Publish('tileUserDataUpdated', id, data, oldData);
};
Graph.prototype.getUserData = function(id) { return this._userData[id]; };
Graph.prototype.tileWeight = function(id) { return this._weights[id]; };
Graph.prototype.setTileWeight = function(id, w) {
    if (this._weights[id] !== w) {
        var oldWeight = this._weights[id];
        this._weights[id] = w;
        this.Publish('tileWeightUpdated', id, w, oldWeight);
    }
};
Graph.prototype.setTileWeights = function(ids, w) {
    var _this = this;
    ids.forEach(function(id) {
        if (_this._weights[id] !== w) {
            var oldWeight = _this._weights[id];
            _this._weights[id] = w;
            this.Publish('tileWeightUpdated', id, w, oldWeight);
        }
    });
};
Graph.prototype.tilesWithGivenWeight = function(w) {
    if (w === void 0) {
        w = Infinity;
    }
    var tiles = [];
    for (var i = 0; i < this.numNodes; i++) {
        if (this._weights[i] === w) {
            tiles.push(i);
        }
    }
    return tiles;
};
Graph.prototype.__addEdge = function(id1, id2, weight) {
    // console.log('addEdge', id1, id2, weight);
    var edgeIndex = this.edgeIndex(id1, id2);
    if (edgeIndex >= 0) {
        if (this.useEdgeWeights) {
            this._edgeWeights[id1][edgeIndex] = weight;
        }
    } else {
        this._edges[id1].push(id2);
        if (this.useEdgeWeights) {
            this._edgeWeights[id2].push(weight);
        }
    }
};

Graph.prototype.__removeEdge = function(id1, id2) {
    var edgeIndex = this.edgeIndex(id1, id2);
    if (edgeIndex >= 0) {
        // console.log('removeEdge', id1, id2);
        this._edges = this._edges.splice(edgeIndex, 1);
        if (this.useEdgeWeights) {
            this._edgeWeights = this._edgeWeights.splice(edgeIndex, 1);
        }
    }
};

Graph.prototype.edgeWeight = function(id1, id2) {
    if (this.useEdgeWeights) {
        var edgeIndex = this.edgeIndex(id1, id2);
        return edgeIndex >= 0 ? this._edgeWeights[id1][edgeIndex] : Infinity;
    } else {
        if (!this.hasEdge(id1, id2)) {
            return Infinity;
        }
        if (this._weights[id2] === undefined) {
            return 1;
        }
        return this._weights[id2];
    }
};
Graph.prototype.distance = function(id1, id2) {
    // Assume edge weight is same as distance (not true!!!)
    return this.edgeWeight(id1, id2);
};
// Is there an edge from id1 to id2?
Graph.prototype.hasEdge = function(id1, id2) {
    // NOTE: Be careful of type of id2
    return this._edges[id1] && this._edges[id1].indexOf(id2) >= 0;
};
// Return index of edge between id1 and id2 (-1 if no edge)
Graph.prototype.edgeIndex = function(id1, id2) {
    // NOTE: Be careful of type of id2
    return this._edges[id1] ? this._edges[id1].indexOf(id2) : -1;
};
// All edges from id
Graph.prototype.edgesFrom = function(id1) {
    var _this = this;
    var edges = this._edges[id1].filter(function(id2) {
        return _this.tileWeight(id2) !== Infinity;
    });
    return edges;
};
// All edges as a list of [id1, id2]
Graph.prototype.allEdges = function() {
    var all = [];
    for (var id1 = 0; id1 < this.numNodes; id1++) {
        this._edges[id1].forEach(function(id2) {
            return all.push([ id1, id2 ]);
        });
    }
    return all;
};
// Represent this graph as a serializable json object
Graph.prototype.toJson = function() {
    var encoding = 'rle';
    var json = {
        type : 'Graph',
        metadata : this.metadata,
        edges : this._edges,
        weights : TypeUtils.arrayToJson(
            this._weights, {encoding : encoding}), // this._weights,
        userData : this._userData
    };
    if (this._tileAttributes) {
        json.tileAttributes = {};
        _.each(this._tileAttributes, function(array, key) {
            json.tileAttributes[key] =
                TypeUtils.arrayToJson(array, {encoding : encoding});
        })
    }
    return json;
};
Graph.prototype.fromJson = function(json, opts) {
    opts = opts || {};
    this.useEdgeWeights = opts.useEdgeWeights;
    this._edges = json.edges || [];
    this._weights = TypeUtils.jsonToArray(json.weights, Infinity);
    this._weights = this._weights.map(function(x) {
        if (x === null)
            return Infinity;
        else
            return x;
    });
    this._userData = json.userData;
    if (json.tileAttributes) {
        var tileAttributes = {};
        _.each(json.tileAttributes, function(array, key) {
            tileAttributes[key] = TypeUtils.jsonToArray(array, NaN);
        });
        this._tileAttributes = tileAttributes;
    }
    this.metadata = json.metadata;
    this.numNodes = this._weights.length;
};

// Returns active portion of this graph (by default, the whole graph!)
Graph.prototype.getActive = function() { return this; };
// Make a proxy graph object, to share some things but override
// some methods for comparison diagrams
Graph.prototype.makeProxy = function() {
    var proxy = {};
    for (var field in this) {
        proxy[field] = this[field];
    }
    return proxy;
};

/**
 * A grid of squares, to be used as a graph.
 * The class creates the structure of the grid; the client can
 * directly set the weights on nodes.
 * @param W {int} Width of grid
 * @param H {int} Height of grid
 * @param opts.dirs What directions each cell has neighbors
 * @param opts.precomputeEdges {boolean} Whether edges should be precomputed
 * (this is deprecated and to be removed)
 * @param opts.reverseEdgeOrder {boolean} Whether edge order should be reversed
 * for nice alternativing stair case pattern
 * @extends nav.Graph
 * @constructor
 * @memberOf nav
 */
function SquareGrid(W, H, opts) {
    opts = opts || {};
    Graph.call(this, W * H, opts);
    this.width = W;
    this.height = H;
    this._dirs = opts.dirs || SquareGrid.DIRS;
    this._precomputeEdges = false;
    this._reverseEdgeOrder = opts.reverseEdgeOrder;
    this.__populateEdges();
    if (this._precomputeEdges) {
        this.Subscribe('tileWeightUpdated', this,
                       this.__onTileWeightUpdated.bind(this));
    }

    this.lvl_stride = 1000000;
}

SquareGrid.prototype = Object.create(Graph.prototype);
SquareGrid.constructor = SquareGrid;

SquareGrid.prototype.__populateEdges = function() {
    // console.time('SquareGrid.__populateEdges');
    var scope = this;
    this._dirIds =
        this._dirs.map(function(d) { return scope.toId(d[0], d[1]); });
    this.dir = this._dirs;
    if (!this._precomputeEdges)
        return;

    var W = this.width;
    var H = this.height;
    if (!this._edges) {
        this._edges = [];       // node id to list of node ids
        this._edgeWeights = []; // node id to list of edge weights
        for (var id = 0; id < this.numNodes; id++) {
            this._edges[id] = [];
            if (this.useEdgeWeights) {
                this._edgeWeights[id] = [];
            }
        }
    }
    for (var x = 0; x < W; x++) {
        for (var y = 0; y < H; y++) {
            var id = this.toId(x, y);
            scope._dirs.forEach(function(dir) {
                var x2 = x + dir[0], y2 = y + dir[1];
                if (scope.isTraversable(x2, y2)) {
                    var id2 = scope.toId(x2, y2);
                    scope._edges[id].push(id2);
                    if (scope.useEdgeWeights) {
                        scope._edgeWeights[id].push(scope.__computeEdgeWeight(
                            dir, scope._weights[id2]));
                    }
                }
            });
        }
    }
    // console.timeEnd('SquareGrid.__populateEdges');
};

SquareGrid.prototype.__computeEdgeWeight = function(dir, targetTileWeight) {
    return Math.sqrt(dir[0] * dir[0] + dir[1] * dir[1]) * targetTileWeight;
};

SquareGrid.prototype.__onTileWeightUpdated = function(id, w, oldWeight) {
    var xy = this.fromId(id);
    // console.log('update tile weight', id, w, xy);
    if (!this.isValidCell(xy[0], xy[1])) {
        // console.log('Skipping invalid cell');
        return;
    }
    var scope = this;
    var cellIsOkay = this.isTraversable(xy[0], xy[1]);
    this._dirs.forEach(function(dir) {
        var x2 = xy[0] + dir[0], y2 = xy[1] + dir[1];
        // Update edges from (x2,y2) to (x1,y1)
        var id2 = scope.toId(x2, y2);
        if (scope.isValidCell(x2, y2)) {
            if (cellIsOkay) {
                scope.__addEdge(id2, id,
                                scope.__computeEdgeWeight([ -dir[0], -dir[1] ],
                                                          scope._weights[id2]));
            } else {
                scope.__removeEdge(id2, id);
            }
        }
    });
};

SquareGrid.prototype.__edgesFrom = function(id1) {
    var edges = [];
    var xy = this.fromId(id1);
    var scope = this;
    this._dirs.forEach(function(dir) {
        var x2 = xy[0] + dir[0], y2 = xy[1] + dir[1];
        if (scope.isTraversable(x2, y2)) {
            edges.push(scope.toId(x2, y2));
        }
    });
    return edges;
};

SquareGrid.prototype.edgesFrom = function(id1) {
    var edges = this.__edgesFrom(id1);
    var xy = this.fromId(id1);
    if (this._reverseEdgeOrder) {
        if ((xy[0] + xy[1]) % 2 === 0) {
            // This is purely for aesthetic purposes on grids -- using a
            // checkerboard pattern, flip every other tile's edges so
            // that paths along diagonal lines end up stair stepping
            // instead of doing all east/west movement first and then
            // all north/south.
            edges.reverse();
        }
    }
    return edges;
};
// Encode/decode grid locations (x,y) to integers (id)
SquareGrid.prototype.isValidCell = function(x, y) {
    return 0 <= x && x < this.width && 0 <= y && y < this.height;
};

SquareGrid.prototype.isTraversable = function(x, y) {
    return this.isValidCell(x, y) &&
           this.tileWeight(this.toId(x, y)) !== Infinity;
};

SquareGrid.prototype.toId = function(x, y, dir = null) {
    return x + y * this.width;
};
SquareGrid.prototype.fromId = function(id) {
    return [ id % this.width, Math.floor(id / this.width), [ 0, 0 ] ];
};
// Represent this grid as a serializable json object
SquareGrid.prototype.toJson = function() {
    var json = Graph.prototype.toJson.call(this);
    json['type'] = 'SquareGrid';
    json['width'] = this.width;
    json['height'] = this.height;
    delete json['edges']; // Don't need to have edges (we can recompute them
                          // since they are computed from the grid structure)
    return json;
};
SquareGrid.prototype.fromJson = function(json, opts) {
    opts = opts || {};
    Graph.prototype.fromJson.call(this, json, opts);
    this.width = json.width;
    this.height = json.height;
    this._dirs = opts.dirs || SquareGrid.DIRS;
    this.__populateEdges();
};
// Encode this grid as a set of pixels
SquareGrid.prototype.toPixels = function(key) {
    if (key) {
        if (this._tileAttributes[key]) {
            return {
                data: this._tileAttributes[key], width: this.width,
                    height: this.height
            }
        }
    } else {
        return {data : this._weights, width : this.width, height : this.height};
    }
};
// Read pixels and use it to initialize the grid
SquareGrid.prototype.fromPixels = function(pixels, opts) {
    this.width = pixels.width;
    this.height = pixels.height;
    this._weights = pixels.data;
    this._dirs = opts.dirs || SquareGrid.DIRS;
    this.__populateEdges();
};
// To encoded
SquareGrid.prototype.toEncodedPixels = function(key, encodeFn) {
    var pixels = this.toPixels(key);
    var data = new Uint8Array(pixels.width * pixels.height * 4);
    for (var i = 0; i < pixels.data.length; i++) {
        var d = pixels.data[i];
        var v = encodeFn(d);
        var j = i * 4;
        data[j] = v[0];
        data[j + 1] = v[1];
        data[j + 2] = v[2];
        data[j + 3] = v[3];
    }
    pixels.data = data;
    pixels.shape = [ pixels.width, pixels.height, 4 ];
    return pixels;
};

// Interface of SquareGrid that don't require precomputation of edges
SquareGrid.prototype.distance = function(id1, id2) {
    var xy1 = this.fromId(id1);
    var xy2 = this.fromId(id2);
    var dx = xy2[0] - xy1[0];
    var dy = xy2[1] - xy1[1];
    var dir = [ dx, dy ];
    return Math.sqrt(dir[0] * dir[0] + dir[1] * dir[1]);
};

SquareGrid.prototype.__line_tile_weight = function(xy0, xy1) {
    var x1 = xy0[0], y1 = xy0[1];
    var x2 = xy1[0], y2 = xy1[1];

    var xstep = 1, ystep = 1;
    var error, errorprev;
    var y = y1, x = x1;
    var ddy, ddx;
    var dx = x2 - x1;
    var dy = y2 - y1;

    if (Math.abs(dy) <= 1 && Math.abs(dx) <= 1) {
        return 0.0;
    }

    if (dy < 0) {
        dy = -dy;
        ystep = -1;
    }

    if (dx < 0) {
        dx = -dx;
        xstep = -1;
    }

    ddy = 2 * dy;
    ddx = 2 * dx;

    var max_tile_weight = 0.0;
    if (ddx >= ddy) {
        errorprev = error = dx;
        for (var i = 0; i < dx; i++) {
            x += xstep;
            error += ddy;
            if (error > ddx) {
                y += ystep;
                error -= ddx;

                if (error + errorprev < ddx) {
                    max_tile_weight =
                        Math.max(max_tile_weight,
                                 this._weights[this.toId(x, y - ystep)]);
                } else if (error + errorprev > ddx) {
                    max_tile_weight =
                        Math.max(max_tile_weight,
                                 this._weights[this.toId(x - xstep, y)]);
                } else {
                    max_tile_weight =
                        Math.max(max_tile_weight,
                                 this._weights[this.toId(x, y - ystep)]);
                    max_tile_weight =
                        Math.max(max_tile_weight,
                                 this._weights[this.toId(x - xstep, y)]);
                }
            }

            max_tile_weight =
                Math.max(max_tile_weight, this._weights[this.toId(x, y)]);

            errorprev = error;
        }
    } else {
        errorprev = error = dy;
        for (var i = 0; i < dy; i++) {
            y += ystep;
            error += ddx;
            if (error > ddy) {
                x += xstep;
                error -= ddy;

                if (error + errorprev < ddy) {
                    max_tile_weight =
                        Math.max(max_tile_weight,
                                 this._weights[this.toId(x - xstep, y)]);
                } else if (error + errorprev > ddy) {
                    max_tile_weight =
                        Math.max(max_tile_weight,
                                 this._weights[this.toId(x, y - ystep)]);
                } else {
                    max_tile_weight =
                        Math.max(max_tile_weight,
                                 this._weights[this.toId(x - xstep, y)]);
                    max_tile_weight =
                        Math.max(max_tile_weight,
                                 this._weights[this.toId(x, y - ystep)]);
                }
            }

            max_tile_weight =
                Math.max(max_tile_weight, this._weights[this.toId(x, y)]);

            errorprev = error;
        }
    }

    return max_tile_weight;
};

SquareGrid.prototype.edgeWeight = function(id1, id2) {
    var xy1 = this.fromId(id1);
    var xy2 = this.fromId(id2);
    var dx = xy2[0] - xy1[0];
    var dy = xy2[1] - xy1[1];
    return this.__computeEdgeWeight(
        [ dx, dy ],
        Math.max(this._weights[id2], this.__line_tile_weight(xy1, xy2)));
};
SquareGrid.prototype.__hasEdgeXY = function(xy1, xy2) {
    var dx = xy2[0] - xy1[0];
    var dy = xy2[1] - xy1[1];
    return this._dirIds.indexOf(this.toId(dx, dy)) >= 0 &&
           this.isTraversable(xy2[0], xy2[1]);
};
// Is there an edge from id1 to id2?
function __ensureInt(x) {
    if (typeof x === 'string') {
        return parseInt(x);
    }
}
SquareGrid.prototype.hasEdge = function(id1, id2) {
    return (this._precomputeEdges)
               ? Graph.prototype.hasEdge.call(this, __ensureInt(id1),
                                              __ensureInt(id2))
               : this.__hasEdgeXY(this.fromId(id1), this.fromId(id2));
};
// Index of edge
SquareGrid.prototype.edgeIndex = function(id1, id2) {
    return Graph.prototype.edgeIndex.call(this, __ensureInt(id1),
                                          __ensureInt(id2))
};
// All edges as a list of [id1, id2]
SquareGrid.prototype.allEdges = function() {
    if (this._precomputeEdges) {
        return Graph.prototype.allEdges.call(this);
    } else {
        var all = [];
        var W = this.width;
        var H = this.height;
        var scope = this;
        for (var x = 0; x < W; x++) {
            for (var y = 0; y < H; y++) {
                var id = this.toId(x, y);
                this._dirs.forEach(function(dir) {
                    var x2 = x + dir[0], y2 = y + dir[1];
                    if (scope.isTraversable(x2, y2)) {
                        all.push([ id, scope.toId(x2, y2) ]);
                    }
                });
            }
        }
        return all;
    }
};

/**
 * A grid of squares, to be used as a graph. Distnace is in action space!
 * The class creates the structure of the grid; the client can
 * directly set the weights on nodes.
 * @param W {int} Width of grid
 * @param H {int} Height of grid
 * @param opts.dirs What directions each cell has neighbors
 * @param opts.precomputeEdges {boolean} Whether edges should be precomputed
 * (this is deprecated and to be removed)
 * @param opts.reverseEdgeOrder {boolean} Whether edge order should be reversed
 * for nice alternativing stair case pattern
 * @extends nav.Graph
 * @constructor
 * @memberOf nav
 */
function ActionGrid(W, H, opts) {
    opts = opts || {};
    SquareGrid.call(this, W, H, opts);
    this.dirs = this._dirs;
    this._n_dirs = this._dirs.length;
    this._rots = (opts.nrots !== undefined) ? opts.nrots : 40;

    this.refinement_factor =
        (opts.refinement_factor !== undefined) ? opts.refinement_factor : 6.0;

    this._cells_per_step = this.refinement_factor;

    this.lvl_stride = this.numNodes;
    this._lvl_stride = this.numNodes;
}

ActionGrid.prototype = Object.create(SquareGrid.prototype);
ActionGrid.constructor = ActionGrid;

ActionGrid.prototype.make_grid_finer = function() {
    new_weights = [];
    var W = this.width;
    var H = this.height;

    for (var y = 0; y < H * this.refinement_factor; y++) {
        var old_y = Math.floor(y / this.refinement_factor);
        for (var x = 0; x < W * this.refinement_factor; x++) {
            var old_x = Math.floor(x / this.refinement_factor);
            var old_id = old_x + old_y * W;
            new_weights.push(this._weights[old_id]);
        }
    }

    this._weights = new_weights;
    this.numNodes = this._weights.length;
    this.width *= this.refinement_factor;
    this.height *= this.refinement_factor;

    this.lvl_stride = this.numNodes;
    this._lvl_stride = this.numNodes;

    this.cellSize /= this.refinement_factor;
    console.log('CELL SIZE: ' + this.cellSize);

    var dirs = [];
    for (var r = 0; r < 2.0 * Math.PI; r += 2.0 * Math.PI / this._rots) {
        dirs.push([
            Math.round(Math.cos(r) * this._cells_per_step),
            Math.round(Math.sin(r) * this._cells_per_step)
        ]);
    }
    this.dirs = dirs;
    this._dirs = dirs;
    this._n_dirs = dirs.length;

    console.log('NUM DIRS: ' + this._n_dirs);
};

function _dot(d1, d2) {
    var v = 0.0;
    for (var i = 0; i < d1.length; i++) {
        v += d1[i] * d2[i];
    }
    return v;
}

function _mag(d) {
    var v = 0.0;
    for (var i = 0; i < d.length; i++) {
        v += d[i] * d[i];
    }
    return Math.sqrt(v);
}

var __eps = 1e-8;
var __1_minus_eps = 1.0 - __eps;
function _angular_distance(d1, d2) {
    var angle = Math.acos(
        Math.max(Math.min(_dot(d1, d2) / (_mag(d1) * _mag(d2)), __1_minus_eps),
                 -__1_minus_eps));
    return angle / Math.PI;
}

ActionGrid.prototype.__computeEdgeWeight = function(id1, id2,
                                                    targetTileWeight) {
    return 1.0 * targetTileWeight;

    var xyd1 = this.fromId(id1);
    var xyd2 = this.fromId(id2);
    var trans = [ xyd2[0] - xyd1[0], xyd2[1] - xyd1[1] ];
    var angle_dist = _angular_distance(xyd1[2], xyd2[2]);
    console.assert(angle_dist >= 0.0 && angle_dist <= 1.0);

    return (_mag(trans) + angle_dist * this._rots / 2.0) * targetTileWeight;
};

ActionGrid.prototype.__edgesFrom = function(id1) {
    var edges = [];
    var lvl = Math.floor(id1 / this._lvl_stride);
    var xy_id = id1 % this._lvl_stride;

    var lvl1 = lvl - 1;
    if (lvl1 < 0) {
        lvl1 = this._n_dirs + lvl1;
    }
    edges.push(lvl1 * this._lvl_stride + xy_id);

    var lvl2 = lvl + 1;
    if (lvl2 >= this._n_dirs) {
        lvl2 = lvl2 - this._n_dirs;
    }
    edges.push(lvl2 * this._lvl_stride + xy_id);

    var xy = this.fromId(id1);
    var dir = xy[2];
    var x2 = xy[0] + dir[0], y2 = xy[1] + dir[1];
    if (this.isTraversable(x2, y2)) {
        edges.push(this.toId(x2, y2, dir));
    }

    /* edges = [];
    var scope = this;
    this._dirs.forEach(function(dir) {
        var x2 = xy[0] + dir[0], y2 = xy[1] + dir[1];
        if (scope.isTraversable(x2, y2)) {
            edges.push(scope.toId(x2, y2, dir));
        }
    }); */

    return edges;
};

ActionGrid.prototype.edgesFrom = function(id1) {
    var edges = this.__edgesFrom(id1);
    return edges;
};

ActionGrid.prototype.toId = function(x, y, dir = null) {
    var lvl = 0;
    if (dir == null) {
        lvl = 0;
    } else {
        var dists =
            _.map(this.dirs, function(d) { return _angular_distance(d, dir); });

        var best_dist = 1.0;
        for (var i = 0; i < dists.length; i++) {
            var d = dists[i];
            if (d < best_dist) {
                lvl = i;
                best_dist = d;
            }
        }
    }

    return lvl * this._lvl_stride + SquareGrid.prototype.toId.call(this, x, y);
};

ActionGrid.prototype.fromId = function(id) {
    var lvl = Math.floor(id / this._lvl_stride);

    console.assert(lvl >= 0 && lvl < this._n_dirs, lvl, this._n_dirs);

    id = id % this._lvl_stride;
    var xy = SquareGrid.prototype.fromId.call(this, id);
    xy[2] = this.dirs[lvl];
    return xy;
};

ActionGrid.prototype.fromJson = function(json, opts) {
    opts = opts || {};
    SquareGrid.prototype.fromJson.call(this, json, opts);

    this.dirs = this._dirs;
    this._n_dirs = this.dirs.length;

    this.lvl_stride = this.numNodes;
    this._lvl_stride = this.numNodes;
};
ActionGrid.prototype.tileWeight = function(id) {
    id = id % this._lvl_stride;
    return this._weights[id];
};
ActionGrid.prototype.setTileWeight = function(id, w) {
    id = id % this._lvl_stride;
    this._weights[id] = w;
};

// Interface of SquareGrid that don't require precomputation of edges
ActionGrid.prototype.distance = function(id1, id2) {
    var xy1 = this.fromId(id1);
    var xy2 = this.fromId(id2);
    var dx = xy2[0] - xy1[0];
    var dy = xy2[1] - xy1[1];
    var dir = [ dx, dy ];
    return Math.sqrt(dir[0] * dir[0] + dir[1] * dir[1]);
};

ActionGrid.prototype.edgeWeight = function(id1, id2) {
    return this.__computeEdgeWeight(id1, id2, this.tileWeight(id2));
};

SquareGrid.DIRS4 = [ [ 1, 0 ], [ 0, 1 ], [ -1, 0 ], [ 0, -1 ] ];
SquareGrid.DIRS8 = [
    [ 1, 0 ], [ 1, 1 ], [ 0, 1 ], [ -1, 1 ], [ -1, 0 ], [ -1, -1 ], [ 0, -1 ],
    [ 1, -1 ]
];
SquareGrid.DIRS = SquareGrid.DIRS4;

module.exports = {
    Graph : Graph,
    SquareGrid : SquareGrid,
    ActionGrid : ActionGrid,
    CellAttribute : CellAttribute
};
