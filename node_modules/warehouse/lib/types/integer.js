'use strict';

const SchemaTypeNumber = require('./number');
const ValidationError = require('../error/validation');

/**
 * Integer schema type.
 */
class SchemaTypeInteger extends SchemaTypeNumber {

  /**
   * Casts a integer.
   *
   * @param {*} value
   * @param {Object} data
   * @return {Number}
   */
  cast(value_, data) {
    const value = super.cast(value_, data);

    return parseInt(value, 10);
  }

  /**
   * Validates an integer.
   *
   * @param {*} value
   * @param {Object} data
   * @return {Number|Error}
   */
  validate(value_, data) {
    const value = super.validate(value_, data);

    if (value % 1 !== 0) {
      throw new ValidationError(`\`${value}\` is not an integer!`);
    }

    return value;
  }
}

module.exports = SchemaTypeInteger;
