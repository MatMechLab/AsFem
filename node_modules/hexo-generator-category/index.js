/* global hexo */

'use strict';

hexo.config.category_generator = Object.assign({
  per_page: typeof hexo.config.per_page === 'undefined' ? 10 : hexo.config.per_page
}, hexo.config.category_generator);

hexo.extend.generator.register('category', require('./lib/generator'));
