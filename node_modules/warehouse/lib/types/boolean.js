'use strict';

const SchemaType = require('../schematype');
const ValidationError = require('../error/validation');

/**
 * Boolean schema type.
 */
class SchemaTypeBoolean extends SchemaType {

  /**
   * Casts a boolean.
   *
   * @param {*} value
   * @param {Object} data
   * @return {Boolean}
   */
  cast(value_, data) {
    const value = super.cast(value_, data);

    if (value === 'false' || value === '0') return false;

    return Boolean(value);
  }

  /**
   * Validates a boolean.
   *
   * @param {*} value
   * @param {Object} data
   * @return {Boolean|Error}
   */
  validate(value_, data) {
    const value = super.validate(value_, data);

    if (value != null && typeof value !== 'boolean') {
      throw new ValidationError(`\`${value}\` is not a boolean!`);
    }

    return value;
  }

  /**
   * Parses data and transform them into boolean values.
   *
   * @param {*} value
   * @param {Object} data
   * @return {Boolean}
   */
  parse(value, data) {
    return Boolean(value);
  }

  /**
   * Transforms data into number to compress the size of database files.
   *
   * @param {Boolean} value
   * @param {Object} data
   * @return {Number}
   */
  value(value, data) {
    return +value;
  }
}

module.exports = SchemaTypeBoolean;
