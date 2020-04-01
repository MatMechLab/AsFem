'use strict';

const WarehouseError = require('../error');

class ValidationError extends WarehouseError {}

ValidationError.prototype.name = 'ValidationError';

module.exports = ValidationError;
