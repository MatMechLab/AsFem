'use strict';

module.exports = function(app) {
  const config = this.config.server || {};
  if (!config.header) return;

  app.use((req, res, next) => {
    res.setHeader('X-Powered-By', 'Hexo');
    next();
  });
};
