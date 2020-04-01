'use strict';

const SchemaType = require('../schematype');
const { setGetter } = require('../util');

/**
 * Virtual schema type.
 */
class SchemaTypeVirtual extends SchemaType {

  /**
   * Add a getter.
   *
   * @param {Function} fn
   * @chainable
   */
  get(fn) {
    if (typeof fn !== 'function') {
      throw new TypeError('Getter must be a function!');
    }

    this.getter = fn;

    return this;
  }

  /**
   * Add a setter.
   *
   * @param {Function} fn
   * @chainable
   */
  set(fn) {
    if (typeof fn !== 'function') {
      throw new TypeError('Setter must be a function!');
    }

    this.setter = fn;

    return this;
  }

  /**
   * Applies getters.
   *
   * @param {*} value
   * @param {Object} data
   * @return {*}
   */
  cast(value, data) {
    if (typeof this.getter !== 'function') return;

    const getter = this.getter;
    let hasCache = false;
    let cache;

    setGetter(data, this.name, () => {
      if (!hasCache) {
        cache = getter.call(data);
        hasCache = true;
      }

      return cache;
    });
  }

  /**
   * Applies setters.
   *
   * @param {*} value
   * @param {Object} data
   */
  validate(value, data) {
    if (typeof this.setter === 'function') {
      this.setter.call(data, value);
    }
  }
}

module.exports = SchemaTypeVirtual;
