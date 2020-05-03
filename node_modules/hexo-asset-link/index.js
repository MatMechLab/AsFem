'use strict';

const url = require('url');
const chalk = require('chalk');

hexo.extend.filter.register('before_post_render', function (data) {
  // Need post asset folder option enabled and asset_dir attribute available
  if (!hexo.config.post_asset_folder || !data.asset_dir) return;
  // Make sure path delimiter is slash rather than backslash.
  let asset_dir = data.asset_dir.replace(/\\/g, '/');
  hexo.log.d('Post asset folder full path:', chalk.magenta(asset_dir));
  // Split hierarchy, filter empty string, last one is asset folder's name.
  let asset_dir_name = asset_dir.split('/').filter(i => i).pop();
  hexo.log.d('Post asset folder name:', chalk.magenta(asset_dir_name));
  // Start with './' or not, end with '/', this is how user write asset links in markdown.
  let path_markdown = RegExp('(\.\/)?' + asset_dir_name + '\/', 'g');
  if (!path_markdown.test(data.content)) return; // no asset link found, do nothing
  // Permalink's pathname, supposed to start with '/'
  let pathname = url.parse(data.permalink).pathname;
  hexo.log.d('Post html path name:', chalk.magenta(pathname));
  // Strip any suffix if exists, supposed to start and end with '/', this is where assets would be in html.
  let path_html = pathname.replace(/\.[^/.]+$/, '/');
  data.content = data.content.replace(path_markdown, path_html);
  hexo.log.i('Path converted:', chalk.yellow(path_markdown.toString()), '→', chalk.green(path_html));
});
