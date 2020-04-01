'use strict';

const WarehouseError = require('../error');

class PopulationError extends WarehouseError {}

PopulationError.prototype.name = 'PopulationError';

module.exports = PopulationError;
