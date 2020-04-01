'use strict';

const pathFn = require('path');
const mime = require('mime');

module.exports = function(app) {
  const { config, route } = this;
  const { args = {} } = this.env;
  const { root } = config;

  if (args.s || args.static) return;

  app.use(root, (req, res, next) => {
    const { method } = req;
    if (method !== 'GET' && method !== 'HEAD') return next();

    let url = route.format(decodeURIComponent(req.url));
    const data = route.get(url);
    const extname = pathFn.extname(url);

    // When the URL is `foo/index.html` but users access `foo`, redirect to `foo/`.
    if (!data) {
      if (extname) return next();

      url = encodeURI(url);
      res.statusCode = 302;
      res.setHeader('Location', `${root + url}/`);
      res.end('Redirecting');
      return;
    }

    res.setHeader('Content-Type', extname ? mime.getType(extname) : 'application/octet-stream');

    if (method === 'GET') {
      data.pipe(res).on('error', next);
    } else {
      res.end();
    }
  });
};
