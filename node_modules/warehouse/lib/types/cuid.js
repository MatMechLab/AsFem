'use strict';

const SchemaType = require('../schematype');
const cuid = require('cuid');
const ValidationError = require('../error/validation');

/**
 * [CUID](https://github.com/ericelliott/cuid) schema type.
 */
class SchemaTypeCUID extends SchemaType {

  /**
   * Casts data. Returns a new CUID only if value is null and the field is
   * required.
   *
   * @param {String} value
   * @param {Object} data
   * @return {String}
   */
  cast(value, data) {
    if (value == null && this.options.required) {
      return cuid();
    }

    return value;
  }

  /**
   * Validates data. A valid CUID must be started with `c` and 25 in length.
   *
   * @param {*} value
   * @param {Object} data
   * @return {String|Error}
   */
  validate(value, data) {
    if (value && (value[0] !== 'c' || value.length !== 25)) {
      throw new ValidationError(`\`${value}\` is not a valid CUID`);
    }

    return value;
  }
}

module.exports = SchemaTypeCUID;
