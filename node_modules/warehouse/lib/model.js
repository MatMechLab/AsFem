'use strict';

const { EventEmitter } = require('events');
const cloneDeep = require('rfdc')();
const Promise = require('bluebird');
const { parseArgs, getProp, setGetter, shuffle } = require('./util');
const Document = require('./document');
const Query = require('./query');
const Schema = require('./schema');
const Types = require('./types');
const WarehouseError = require('./error');
const PopulationError = require('./error/population');
const Mutex = require('./mutex');

class Model extends EventEmitter {

  /**
   * Model constructor.
   *
   * @param {string} name Model name
   * @param {Schema|object} [schema] Schema
   */
  constructor(name, schema_) {
    super();

    let schema;

    // Define schema
    if (schema_ instanceof Schema) {
      schema = schema_;
    } else if (typeof schema_ === 'object') {
      schema = new Schema(schema_);
    } else {
      schema = new Schema();
    }

    // Set `_id` path for schema
    if (!schema.path('_id')) {
      schema.path('_id', {type: Types.CUID, required: true});
    }

    this.name = name;
    this.data = {};
    this._mutex = new Mutex();
    this.schema = schema;
    this.length = 0;

    class _Document extends Document {
      constructor(data) {
        super(data);

        // Apply getters
        schema._applyGetters(this);
      }
    }

    this.Document = _Document;

    _Document.prototype._model = this;
    _Document.prototype._schema = schema;

    class _Query extends Query {}

    this.Query = _Query;

    _Query.prototype._model = this;
    _Query.prototype._schema = schema;

    // Apply static methods
    Object.assign(this, schema.statics);

    // Apply instance methods
    Object.assign(_Document.prototype, schema.methods);
  }

  /**
   * Creates a new document.
   *
   * @param {object} data
   * @return {Document}
   */
  new(data) {
    return new this.Document(data);
  }

  /**
   * Finds a document by its identifier.
   *
   * @param {*} id
   * @param {object} options
   *   @param {boolean} [options.lean=false] Returns a plain JavaScript object
   * @return {Document|object}
   */
  findById(id, options_) {
    const raw = this.data[id];
    if (!raw) return;

    const options = Object.assign({
      lean: false
    }, options_);

    const data = cloneDeep(raw);
    return options.lean ? data : this.new(data);
  }

  /**
   * Checks if the model contains a document with the specified id.
   *
   * @param {*} id
   * @return {boolean}
   */
  has(id) {
    return Boolean(this.data[id]);
  }

  /**
   * Acquires write lock.
   *
   * @param {*} id
   * @return {Promise}
   * @private
   */
  _acquireWriteLock(id) {
    const mutex = this._mutex;

    return new Promise((resolve, reject) => {
      mutex.lock(resolve);
    }).disposer(() => {
      mutex.unlock();
    });
  }

  /**
   * Inserts a document.
   *
   * @param {Document|object} data
   * @return {Promise}
   * @private
   */
  _insertOne(data_) {
    const schema = this.schema;

    // Apply getters
    const data = data_ instanceof this.Document ? data_ : this.new(data_);
    const id = data._id;

    // Check ID
    if (!id) {
      return Promise.reject(new WarehouseError('ID is not defined', WarehouseError.ID_UNDEFINED));
    }

    if (this.has(id)) {
      return Promise.reject(new WarehouseError('ID `' + id + '` has been used', WarehouseError.ID_EXIST));
    }

    // Apply setters
    const result = data.toObject();
    schema._applySetters(result);

    // Pre-hooks
    return execHooks(schema, 'pre', 'save', data).then(data => {
      // Insert data
      this.data[id] = result;
      this.length++;

      this.emit('insert', data);
      return execHooks(schema, 'post', 'save', data);
    });
  }

  /**
   * Inserts a document.
   *
   * @param {object} data
   * @param {function} [callback]
   * @return {Promise}
   */
  insertOne(data, callback) {
    return Promise.using(this._acquireWriteLock(), () => this._insertOne(data)).asCallback(callback);
  }

  /**
   * Inserts documents.
   *
   * @param {object|array} data
   * @param {function} [callback]
   * @return {Promise}
   */
  insert(data, callback) {
    if (Array.isArray(data)) {
      return Promise.mapSeries(data, item => this.insertOne(item)).asCallback(callback);
    }

    return this.insertOne(data, callback);
  }

  /**
   * Inserts the document if it does not exist; otherwise updates it.
   *
   * @param {object} data
   * @param {function} [callback]
   * @return {Promise}
   */
  save(data, callback) {
    const id = data._id;

    if (!id) return this.insertOne(data, callback);

    return Promise.using(this._acquireWriteLock(), () => {
      if (this.has(id)) {
        return this._replaceById(id, data);
      }

      return this._insertOne(data);
    }).asCallback(callback);
  }

  /**
   * Updates a document with a compiled stack.
   *
   * @param {*} id
   * @param {array} stack
   * @return {Promise}
   * @private
   */
  _updateWithStack(id, stack) {
    const schema = this.schema;

    const data = this.data[id];

    if (!data) {
      return Promise.reject(new WarehouseError('ID `' + id + '` does not exist', WarehouseError.ID_NOT_EXIST));
    }

    // Clone data
    let result = cloneDeep(data);

    // Update
    for (let i = 0, len = stack.length; i < len; i++) {
      stack[i](result);
    }

    // Apply getters
    const doc = this.new(result);

    // Apply setters
    result = doc.toObject();
    schema._applySetters(result);

    // Pre-hooks
    return execHooks(schema, 'pre', 'save', doc).then(data => {
      // Update data
      this.data[id] = result;

      this.emit('update', data);
      return execHooks(schema, 'post', 'save', data);
    });
  }

  /**
   * Finds a document by its identifier and update it.
   *
   * @param {*} id
   * @param {object} update
   * @param {function} [callback]
   * @return {Promise}
   */
  updateById(id, update, callback) {
    return Promise.using(this._acquireWriteLock(), () => {
      const stack = this.schema._parseUpdate(update);
      return this._updateWithStack(id, stack);
    }).asCallback(callback);
  }

  /**
   * Updates matching documents.
   *
   * @param {object} query
   * @param {object} data
   * @param {function} [callback]
   * @return {Promise}
   */
  update(query, data, callback) {
    return this.find(query).update(data, callback);
  }

  /**
   * Finds a document by its identifier and replace it.
   *
   * @param {*} id
   * @param  {object} data
   * @return {Promise}
   * @private
   */
  _replaceById(id, data_) {
    const schema = this.schema;

    if (!this.has(id)) {
      return Promise.reject(new WarehouseError('ID `' + id + '` does not exist', WarehouseError.ID_NOT_EXIST));
    }

    data_._id = id;

    // Apply getters
    const data = data_ instanceof this.Document ? data_ : this.new(data_);

    // Apply setters
    const result = data.toObject();
    schema._applySetters(result);

    // Pre-hooks
    return execHooks(schema, 'pre', 'save', data).then(data => {
      // Replace data
      this.data[id] = result;

      this.emit('update', data);
      return execHooks(schema, 'post', 'save', data);
    });
  }

  /**
   * Finds a document by its identifier and replace it.
   *
   * @param {*} id
   * @param {object} data
   * @param {function} [callback]
   * @return {Promise}
   */
  replaceById(id, data, callback) {
    return Promise.using(this._acquireWriteLock(), () => this._replaceById(id, data)).asCallback(callback);
  }

  /**
   * Replaces matching documents.
   *
   * @param {object} query
   * @param {object} data
   * @param {function} [callback]
   * @return {Promise}
   */
  replace(query, data, callback) {
    return this.find(query).replace(data, callback);
  }

  /**
   * Finds a document by its identifier and remove it.
   *
   * @param {*} id
   * @param {function} [callback]
   * @return {Promise}
   * @private
   */
  _removeById(id) {
    const schema = this.schema;

    const data = this.data[id];

    if (!data) {
      return Promise.reject(new WarehouseError('ID `' + id + '` does not exist', WarehouseError.ID_NOT_EXIST));
    }

    // Pre-hooks
    return execHooks(schema, 'pre', 'remove', data).then(data => {
      // Remove data
      this.data[id] = null;
      this.length--;

      this.emit('remove', data);
      return execHooks(schema, 'post', 'remove', data);
    });
  }

  /**
   * Finds a document by its identifier and remove it.
   *
   * @param {*} id
   * @param {function} [callback]
   * @return {Promise}
   */
  removeById(id, callback) {
    return Promise.using(this._acquireWriteLock(), () => this._removeById(id)).asCallback(callback);
  }

  /**
   * Removes matching documents.
   *
   * @param {object} query
   * @param {object} [callback]
   * @return {Promise}
   */
  remove(query, callback) {
    return this.find(query).remove(callback);
  }

  /**
   * Deletes a model.
   */
  destroy() {
    this._database._models[this.name] = null;
  }

  /**
   * Returns the number of elements.
   *
   * @return {number}
   */
  count() {
    return this.length;
  }

  /**
   * Iterates over all documents.
   *
   * @param {function} iterator
   * @param {object} [options] See {@link Model#findById}.
   */
  forEach(iterator, options) {
    const keys = Object.keys(this.data);
    let num = 0;

    for (let i = 0, len = keys.length; i < len; i++) {
      const data = this.findById(keys[i], options);
      if (data) iterator(data, num++);
    }
  }

  /**
   * Returns an array containing all documents.
   *
   * @param {Object} [options] See {@link Model#findById}.
   * @return {Array}
   */
  toArray(options) {
    const result = new Array(this.length);

    this.forEach((item, i) => {
      result[i] = item;
    }, options);

    return result;
  }

  /**
   * Finds matching documents.
   *
   * @param {Object} query
   * @param {Object} [options]
   *   @param {Number} [options.limit=0] Limits the number of documents returned.
   *   @param {Number} [options.skip=0] Skips the first elements.
   *   @param {Boolean} [options.lean=false] Returns a plain JavaScript object.
   * @return {Query|Array}
   */
  find(query, options_) {
    const options = options_ || {};
    const filter = this.schema._execQuery(query);
    const keys = Object.keys(this.data);
    const len = keys.length;
    let limit = options.limit || this.length;
    let skip = options.skip;
    const data = this.data;
    const arr = [];

    for (let i = 0; limit && i < len; i++) {
      const key = keys[i];
      const item = data[key];

      if (item && filter(item)) {
        if (skip) {
          skip--;
        } else {
          arr.push(this.findById(key, options));
          limit--;
        }
      }
    }

    return options.lean ? arr : new this.Query(arr);
  }

  /**
   * Finds the first matching documents.
   *
   * @param {Object} query
   * @param {Object} [options]
   *   @param {Number} [options.skip=0] Skips the first elements.
   *   @param {Boolean} [options.lean=false] Returns a plain JavaScript object.
   * @return {Document|Object}
   */
  findOne(query, options_) {
    const options = options_ || {};
    options.limit = 1;

    const result = this.find(query, options);
    return options.lean ? result[0] : result.data[0];
  }

  /**
   * Sorts documents. See {@link Query#sort}.
   *
   * @param {String|Object} orderby
   * @param {String|Number} [order]
   * @return {Query}
   */
  sort(orderby, order) {
    const sort = parseArgs(orderby, order);
    const fn = this.schema._execSort(sort);

    return new this.Query(this.toArray().sort(fn));
  }

  /**
   * Returns the document at the specified index. `num` can be a positive or
   * negative number.
   *
   * @param {Number} i
   * @param {Object} [options] See {@link Model#findById}.
   * @return {Document|Object}
   */
  eq(i_, options) {
    let index = i_ < 0 ? this.length + i_ : i_;
    const data = this.data;
    const keys = Object.keys(data);

    for (let i = 0, len = keys.length; i < len; i++) {
      const key = keys[i];
      const item = data[key];

      if (!item) continue;

      if (index) {
        index--;
      } else {
        return this.findById(key, options);
      }
    }
  }

  /**
   * Returns the first document.
   *
   * @param {Object} [options] See {@link Model#findById}.
   * @return {Document|Object}
   */
  first(options) {
    return this.eq(0, options);
  }

  /**
   * Returns the last document.
   *
   * @param {Object} [options] See {@link Model#findById}.
   * @return {Document|Object}
   */
  last(options) {
    return this.eq(-1, options);
  }

  /**
   * Returns the specified range of documents.
   *
   * @param {Number} start
   * @param {Number} [end]
   * @return {Query}
   */
  slice(start_, end_) {
    const total = this.length;

    let start = start_ | 0;
    if (start < 0) start += total;
    if (start > total - 1) return new this.Query([]);

    let end = end_ | 0 || total;
    if (end < 0) end += total;

    let len = start > end ? 0 : end - start;
    if (len > total) len = total - start;
    if (!len) return new this.Query([]);

    const arr = new Array(len);
    const keys = Object.keys(this.data);
    const keysLen = keys.length;
    let num = 0;

    for (let i = 0; num < len && i < keysLen; i++) {
      const data = this.findById(keys[i]);
      if (!data) continue;

      if (start) {
        start--;
      } else {
        arr[num++] = data;
      }
    }

    return new this.Query(arr);
  }

  /**
   * Limits the number of documents returned.
   *
   * @param {Number} i
   * @return {Query}
   */
  limit(i) {
    return this.slice(0, i);
  }

  /**
   * Specifies the number of items to skip.
   *
   * @param {Number} i
   * @return {Query}
   */
  skip(i) {
    return this.slice(i);
  }

  /**
   * Returns documents in a reversed order.
   *
   * @return {Query}
   */
  reverse() {
    return new this.Query(this.toArray().reverse());
  }

  /**
   * Returns documents in random order.
   *
   * @return {Query}
   */
  shuffle() {
    return new this.Query(shuffle(this.toArray()));
  }

  /**
   * Creates an array of values by iterating each element in the collection.
   *
   * @param {Function} iterator
   * @param {Object} [options]
   * @return {Array}
   */
  map(iterator, options) {
    const result = new Array(this.length);
    const keys = Object.keys(this.data);
    const len = keys.length;

    for (let i = 0, num = 0; i < len; i++) {
      const data = this.findById(keys[i], options);
      if (data) {
        result[num] = iterator(data, num);
        num++;
      }
    }

    return result;
  }

  /**
   * Reduces a collection to a value which is the accumulated result of iterating
   * each element in the collection.
   *
   * @param {Function} iterator
   * @param {*} [initial] By default, the initial value is the first document.
   * @return {*}
   */
  reduce(iterator, initial) {
    const arr = this.toArray();
    const len = this.length;
    let i, result;

    if (initial === undefined) {
      i = 1;
      result = arr[0];
    } else {
      i = 0;
      result = initial;
    }

    for (; i < len; i++) {
      result = iterator(result, arr[i], i);
    }

    return result;
  }

  /**
   * Reduces a collection to a value which is the accumulated result of iterating
   * each element in the collection from right to left.
   *
   * @param {Function} iterator
   * @param {*} [initial] By default, the initial value is the last document.
   * @return {*}
   */
  reduceRight(iterator, initial) {
    const arr = this.toArray();
    const len = this.length;
    let i, result;

    if (initial === undefined) {
      i = len - 2;
      result = arr[len - 1];
    } else {
      i = len - 1;
      result = initial;
    }

    for (; i >= 0; i--) {
      result = iterator(result, arr[i], i);
    }

    return result;
  }

  /**
   * Creates a new array with all documents that pass the test implemented by the
   * provided function.
   *
   * @param {Function} iterator
   * @param {Object} [options]
   * @return {Query}
   */
  filter(iterator, options) {
    const arr = [];

    this.forEach((item, i) => {
      if (iterator(item, i)) arr.push(item);
    }, options);

    return new this.Query(arr);
  }

  /**
   * Tests whether all documents pass the test implemented by the provided
   * function.
   *
   * @param {Function} iterator
   * @return {Boolean}
   */
  every(iterator) {
    const keys = Object.keys(this.data);
    const len = keys.length;
    let num = 0;

    if (!len) return true;

    for (let i = 0; i < len; i++) {
      const data = this.findById(keys[i]);

      if (data) {
        if (!iterator(data, num++)) return false;
      }
    }

    return true;
  }

  /**
   * Tests whether some documents pass the test implemented by the provided
   * function.
   *
   * @param {Function} iterator
   * @return {Boolean}
   */
  some(iterator) {
    const keys = Object.keys(this.data);
    const len = keys.length;
    let num = 0;

    if (!len) return false;

    for (let i = 0; i < len; i++) {
      const data = this.findById(keys[i]);

      if (data) {
        if (iterator(data, num++)) return true;
      }
    }

    return false;
  }

  /**
   * Returns a getter function for normal population.
   *
   * @param {Object} data
   * @param {Model} model
   * @param {Object} options
   * @return {Function}
   * @private
   */
  _populateGetter(data, model, options) {
    let hasCache = false;
    let cache;

    return () => {
      if (!hasCache) {
        cache = model.findById(data);
        hasCache = true;
      }

      return cache;
    };
  }

  /**
   * Returns a getter function for array population.
   *
   * @param {Object} data
   * @param {Model} model
   * @param {Object} options
   * @return {Function}
   * @private
   */
  _populateGetterArray(data, model, options) {
    const Query = model.Query;
    let hasCache = false;
    let cache;

    return () => {
      if (!hasCache) {
        let arr = [];

        for (let i = 0, len = data.length; i < len; i++) {
          arr.push(model.findById(data[i]));
        }

        if (options.match) {
          cache = new Query(arr).find(options.match, options);
        } else if (options.skip) {
          if (options.limit) {
            arr = arr.slice(options.skip, options.skip + options.limit);
          } else {
            arr = arr.slice(options.skip);
          }

          cache = new Query(arr);
        } else if (options.limit) {
          cache = new Query(arr.slice(0, options.limit));
        } else {
          cache = new Query(arr);
        }

        if (options.sort) {
          cache = cache.sort(options.sort);
        }

        hasCache = true;
      }

      return cache;
    };
  }

  /**
   * Populates document references with a compiled stack.
   *
   * @param {Object} data
   * @param {Array} stack
   * @return {Object}
   * @private
   */
  _populate(data, stack) {
    const models = this._database._models;

    for (let i = 0, len = stack.length; i < len; i++) {
      const item = stack[i];
      const model = models[item.model];

      if (!model) {
        throw new PopulationError('Model `' + item.model + '` does not exist');
      }

      const path = item.path;
      const prop = getProp(data, path);

      if (Array.isArray(prop)) {
        setGetter(data, path, this._populateGetterArray(prop, model, item));
      } else {
        setGetter(data, path, this._populateGetter(prop, model, item));
      }
    }

    return data;
  }

  /**
   * Populates document references.
   *
   * @param {String|Object} path
   * @return {Query}
   */
  populate(path) {
    if (!path) throw new TypeError('path is required');

    const stack = this.schema._parsePopulate(path);
    const arr = new Array(this.length);

    this.forEach((item, i) => {
      arr[i] = this._populate(item, stack);
    });

    return new Query(arr);
  }

  /**
   * Imports data.
   *
   * @param {Array} arr
   * @private
   */
  _import(arr) {
    const len = arr.length;
    const data = this.data;
    const schema = this.schema;

    for (let i = 0; i < len; i++) {
      const item = arr[i];
      data[item._id] = schema._parseDatabase(item);
    }

    this.length = len;
  }

  /**
   * Exports data.
   *
   * @return {String}
   * @private
   */
  _export() {
    return JSON.stringify(this.toJSON());
  }

  toJSON() {
    const result = new Array(this.length);
    const { data, schema } = this;
    const keys = Object.keys(data);
    const { length } = keys;

    for (let i = 0, num = 0; i < length; i++) {
      const raw = data[keys[i]];
      if (raw) {
        result[num++] = schema._exportDatabase(cloneDeep(raw));
      }
    }
    return result;
  }
}

Model.prototype.get = Model.prototype.findById;

function execHooks(schema, type, event, data) {
  const hooks = schema.hooks[type][event];
  if (!hooks.length) return Promise.resolve(data);

  return Promise.each(hooks, hook => hook(data)).thenReturn(data);
}

Model.prototype.size = Model.prototype.count;

Model.prototype.each = Model.prototype.forEach;

Model.prototype.random = Model.prototype.shuffle;

module.exports = Model;
