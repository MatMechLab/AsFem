'use strict';

const morgan = require('morgan');

module.exports = function(app) {
  const { config } = this;
  const { args = {} } = this.env;
  let logger = args.l || args.log || config.server.log;

  if (!logger && !args.debug) return;
  if (typeof logger !== 'string') logger = 'dev';

  app.use(morgan(logger));
};
