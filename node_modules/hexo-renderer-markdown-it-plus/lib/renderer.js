/**
 * @Author: CHC
 * @Date:   2017-07-24T14:15:39+08:00
 * @Email:  chenhuachaoxyz@gmail.com
 * @Filename: renderer.js
 * @Last modified by:   CHC
 * @Last modified time: 2018-03-18T11:22:07+08:00
 * @License: MIT
 * @Copyright: 2017
 */

'use strict';

var hljs = require('highlight.js');

var def_pugs_lst = [
    'markdown-it-emoji',
    'markdown-it-sub',
    'markdown-it-sup',
    'markdown-it-deflist',
    'markdown-it-abbr',
    'markdown-it-footnote',
    'markdown-it-ins',
    'markdown-it-mark',
    // 'markdown-it-katex',
    '@iktakahiro/markdown-it-katex',
    'markdown-it-toc-and-anchor'
];

var utils = require('markdown-it/lib/common/utils');

/**
 * General Default markdown-it config.
 * @param  {Object} config configuration
 * @param  {Object} res    configuration
 * @param  {String} key    the key of config
 * @param  {String} trueVal    the val of config[key] == true
 * @param  {String} falseVal    the val of config[key] == false
 * @return {null}
 */
function checkValue(config, res, key, trueVal, falseVal) {
    res[key] = (config[key] == true || config[key] == undefined || config[key] == null) ? trueVal: falseVal;
}

/**
 * markdown-it default config
 * @param  {Object} config configuration
 * @return {Object}        valid configuration
 */
function checkConfig(config, md) {
    var _res = {};
    var hljs_class = config['hljs_class'];
    var pre_class = config['pre_class'];
    if(!pre_class) pre_class = 'highlight';
    checkValue(config, _res, 'highlight', function(str, lang) {
        if (lang && hljs.getLanguage(lang)) {
            try {
                return '<pre class="' + pre_class + '"><code class="' + lang + '">' + hljs.highlight(lang, str, true).value + '</code></pre>';
            } catch (__) {}
        }
        return '<pre class="' + pre_class + '"><code class="' + lang + '">' + utils.escapeHtml(str) + '</code></pre>';
    }, function(str, lang) {
        return '<pre class="' + pre_class + '"><code class="' + lang + '">' + utils.escapeHtml(str) + '</code></pre>';
    });
    checkValue(config, _res, 'html', true, false);
    checkValue(config, _res, 'xhtmlOut', true, false);
    checkValue(config, _res, 'breaks', true, false);
    checkValue(config, _res, 'linkify', true, false);
    checkValue(config, _res, 'typographer', true, false);
    _res['langPrefix'] = config['langPrefix'] ? config['langPrefix'] : '';
    _res['quotes'] = config['quotes'] ? config['quotes'] : '“”‘’';
    return _res;
}

/**
 * General default plugin config
 * @param  {List} pugs plugin List.
 * @return {List}        plugin List.
 */
function checkPlugins(pugs) {
    var def_pugs_obj = {};
    for(var i = 0;i < def_pugs_lst.length;i++)
        def_pugs_obj[def_pugs_lst[i]] = {'name': def_pugs_lst[i], 'enable': true};
    var _t = [];
    for(var i = 0;i < pugs.length;i++) {
        if(!(pugs[i] instanceof Object) || !(pugs[i].plugin instanceof Object)) continue;
        var pug_name = pugs[i].plugin.name;
        if(!pug_name) continue;
        if(pugs[i].plugin.enable == null || pugs[i].plugin.enable == undefined || pugs[i].plugin.enable != true)
            pugs[i].plugin.enable = false;
        if(def_pugs_obj[pug_name]) {
            def_pugs_obj[pug_name] = pugs[i].plugin;
        }
        else _t.push(pugs[i].plugin);
    }

    for(var i = def_pugs_lst.length - 1;i >= 0;i--) {
        _t.unshift(def_pugs_obj[def_pugs_lst[i]]);
    }
    return _t;
}

module.exports = function(data, options) {
    var config = this.config.markdown_it_plus;
    if(!(config instanceof Object)) {
        config = {};
    }
    var parseConfig = checkConfig(config)
    var md = require('markdown-it')(parseConfig);
    if(config.plugins == undefined || config.plugins == null) {
        config.plugins = [];
    }

    // config.plugins =
    var plugins = checkPlugins(config.plugins);

    md = plugins.reduce(function (md, pugs) {
        if(pugs.enable) {
            if(pugs.name == 'markdown-it-toc-and-anchor') {
                if(pugs.options == null) pugs.options = {}
                if(!pugs.options.anchorLinkSymbol) pugs.options.anchorLinkSymbol = '';
                if(!pugs.options.tocFirstLevel) pugs.options.tocFirstLevel = 2;
                return md.use(require('./markdown-it-toc-and-anchor/index.js').default, pugs.options);
            }
            else {
                let plugin = require(pugs.name);
                if(typeof plugin !== 'function' && typeof plugin.default === 'function')
                    plugin = plugin.default;
                if(pugs.options) return md.use(plugin, pugs.options);
                else return md.use(plugin);
            }
        }
        else return md;
    }, md);

    return md.render(data.text);
}
