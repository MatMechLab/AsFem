'use strict';

const compress = require('compression');

module.exports = function(app) {
  const config = this.config.server || {};
  if (!config.compress) return;

  app.use(compress());
};
