/* global hexo */
'use strict';

hexo.config.tag_generator = Object.assign({
  per_page: hexo.config.per_page == null ? 10 : hexo.config.per_page
}, hexo.config.tag_generator);

hexo.extend.generator.register('tag', require('./lib/generator'));
