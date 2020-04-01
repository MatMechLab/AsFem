'use strict';

const cloneDeep = require('rfdc')();

class Document {

  /**
   * Document constructor.
   *
   * @param {object} data
   */
  constructor(data) {
    if (data) {
      Object.assign(this, data);
    }
  }

  /**
   * Saves the document.
   *
   * @param {function} [callback]
   * @return {Promise}
   */
  save(callback) {
    return this._model.save(this, callback);
  }

  /**
   * Updates the document.
   *
   * @param {object} data
   * @param {function} [callback]
   * @return {Promise}
   */
  update(data, callback) {
    return this._model.updateById(this._id, data, callback);
  }

  /**
   * Replaces the document.
   *
   * @param {object} data
   * @param {function} [callback]
   * @return {Promise}
   */
  replace(data, callback) {
    return this._model.replaceById(this._id, data, callback);
  }

  /**
   * Removes the document.
   *
   * @param {function} [callback]
   * @return {Promise}
   */
  remove(callback) {
    return this._model.removeById(this._id, callback);
  }

  /**
   * Returns a plain JavaScript object.
   *
   * @return {object}
   */
  toObject() {
    const keys = Object.keys(this);
    const obj = {};

    for (let i = 0, len = keys.length; i < len; i++) {
      const key = keys[i];
      // Don't deep clone getters in order to avoid "Maximum call stack size
      // exceeded" error
      obj[key] = isGetter(this, key) ? this[key] : cloneDeep(this[key]);
    }

    return obj;
  }

  /**
   * Returns a string representing the document.
   *
   * @return {String}
   */
  toString() {
    return JSON.stringify(this);
  }

  /**
   * Populates document references.
   *
   * @param {String|Object} expr
   * @return {Document}
   */
  populate(expr) {
    const stack = this._schema._parsePopulate(expr);
    return this._model._populate(this, stack);
  }
}

function isGetter(obj, key) {
  return Object.getOwnPropertyDescriptor(obj, key).get;
}

module.exports = Document;
