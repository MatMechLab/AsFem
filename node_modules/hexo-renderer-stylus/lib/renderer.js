'use strict';

const stylus = require('stylus');

function getProperty(obj, name) {
  name = name.replace(/\[(\w+)\]/g, '.$1').replace(/^\./, '');

  const split = name.split('.');
  let key = split.shift();

  if (!Object.prototype.hasOwnProperty.call(obj, key)) return '';

  let result = obj[key];
  const len = split.length;

  if (!len) return result || '';
  if (typeof result !== 'object') return '';

  for (let i = 0; i < len; i++) {
    key = split[i];
    if (!Object.prototype.hasOwnProperty.call(result, key)) return '';

    result = result[split[i]];
    if (typeof result !== 'object') return result;
  }

  return result;
}

function applyPlugins(stylusConfig, plugins) {
  plugins.forEach(plugin => {
    const factoryFn = require(plugin.trim());
    stylusConfig.use(factoryFn());
  });
}

function stylusFn(data, options, callback) {
  const config = this.config.stylus || {};
  const self = this;
  const plugins = ['nib'].concat(config.plugins || []);

  function defineConfig(style) {
    style.define('hexo-config', data => {
      return getProperty(self.theme.config, data.val);
    });
  }

  const stylusConfig = stylus(data.text);

  applyPlugins(stylusConfig, plugins);

  stylusConfig
    .use(defineConfig)
    .use(style => this.execFilterSync('stylus:renderer', style, {context: this}))
    .set('filename', data.path)
    .set('sourcemap', config.sourcemaps)
    .set('compress', config.compress)
    .set('include css', true)
    .render(callback);
}

stylusFn.disableNunjucks = true;

module.exports = stylusFn;
