/* global hexo */

'use strict';

const renderer = require('./lib/renderer');

hexo.config.marked = Object.assign({
  gfm: true,
  pedantic: false,
  breaks: true,
  smartLists: true,
  smartypants: true,
  modifyAnchors: 0,
  autolink: true,
  mangle: true,
  sanitizeUrl: false,
  headerIds: true,
  anchorAlias: false,
  lazyload: false,
  // TODO: enable prependRoot by default in v4
  prependRoot: false,
  postAsset: false,
  external_link: {
    enable: false,
    exclude: [],
    nofollow: false
  }
}, hexo.config.marked);

renderer.disableNunjucks = Boolean(hexo.config.marked.disableNunjucks);

hexo.extend.renderer.register('md', 'html', renderer, true);
hexo.extend.renderer.register('markdown', 'html', renderer, true);
hexo.extend.renderer.register('mkd', 'html', renderer, true);
hexo.extend.renderer.register('mkdn', 'html', renderer, true);
hexo.extend.renderer.register('mdwn', 'html', renderer, true);
hexo.extend.renderer.register('mdtxt', 'html', renderer, true);
hexo.extend.renderer.register('mdtext', 'html', renderer, true);
