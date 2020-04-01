'use strict';

module.exports = function(app) {
  const { root } = this.config;
  if (root === '/') return;

  // If root url is not `/`, redirect to the correct root url
  app.use((req, res, next) => {
    if (req.method !== 'GET' || req.url !== '/') return next();

    res.statusCode = 302;
    res.setHeader('Location', root);
    res.end('Redirecting');
  });
};
