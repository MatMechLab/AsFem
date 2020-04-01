/* global hexo */

'use strict';

const renderer = require('./lib/renderer');

hexo.extend.renderer.register('ejs', 'html', renderer, true);
