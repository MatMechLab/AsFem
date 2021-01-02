'use strict';

const pathFn = require('path');
const mime = require('mime');

module.exports = function(app) {
  const { config, route } = this;
  const { args = {} } = this.env;
  const { root } = config;

  if (args.s || args.static) return;

  const { pretty_urls } = config;
  const { trailing_index, trailing_html } = pretty_urls ? pretty_urls : {};

  app.use(root, (req, res, next) => {
    const { method, url: requestUrl } = req;
    if (method !== 'GET' && method !== 'HEAD') return next();

    let url = route.format(decodeURIComponent(requestUrl));
    let data = route.get(url);
    let extname = pathFn.extname(url);

    if (!data) {
      if (route.get(url + '.html')) {
        // location `foo/bar.html`; request `foo/bar`; proxy to the location
        extname = '.html';
        data = route.get(url + extname);
        res.setHeader('Content-Type', 'text/html');
        req.url = encodeURI('/' + url + extname);
        data.pipe(res).on('error', next);
        return;
      } else if (route.get(url + '/index.html')) {
        // location `foo/index.html`; request `foo`; redirect to `foo/`
        url = encodeURI(url);
        res.statusCode = 301;
        res.setHeader('Location', `${root + url}/`);
        res.end('Redirecting');
        return;
      } else if (route.get(url.replace(/\/index\.html$/i, '') + '.html')) {
        // location `foo/bar.html`; request `foo/bar/`; redirect to `foo`
        // request with trailing slash is appended with index.html by route.format()
        url = encodeURI(url.replace(/\/index\.html$/i, ''));
        res.statusCode = 301;
        res.setHeader('Location', root + url);
        res.end('Redirecting');
        return;
      } return next();
    }
    if (trailing_html === false && !requestUrl.endsWith('/index.html') && requestUrl.endsWith('.html')) {
      // location `foo/bar.html`; request `foo/bar.html`; redirect to `foo/bar`
      url = encodeURI(url.replace(/\.html$/i, ''));
      res.statusCode = 301;
      res.setHeader('Location', root + url);
      res.end('Redirecting');
      return;
    } else if (trailing_index === false && requestUrl.endsWith('/index.html')) {
      // location `foo/index.html`; request `foo/index.html`; redirect to `foo/`
      url = encodeURI(url.replace(/index\.html$/i, ''));
      res.statusCode = 301;
      res.setHeader('Location', root + url);
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
