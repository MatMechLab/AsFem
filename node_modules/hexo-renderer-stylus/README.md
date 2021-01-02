# hexo-renderer-stylus

[![Build Status](https://travis-ci.org/hexojs/hexo-renderer-stylus.svg?branch=master)](https://travis-ci.org/hexojs/hexo-renderer-stylus)
[![NPM version](https://badge.fury.io/js/hexo-renderer-stylus.svg)](https://www.npmjs.com/package/hexo-renderer-stylus)
[![Coverage Status](https://img.shields.io/coveralls/hexojs/hexo-renderer-stylus.svg)](https://coveralls.io/r/hexojs/hexo-renderer-stylus?branch=master)

Add support for [Stylus] with [nib] and other plugins.

## Install

Prerequisites:
- Hexo 3: >= 0.2
- Hexo 2: 0.1.x

``` bash
$ npm install hexo-renderer-stylus --save
```

## Options

You can configure this plugin in `_config.yml`.

``` yaml
stylus:
  compress: false
  sourcemaps:
    comment: true
    inline: true
    sourceRoot: ''
    basePath: .
  plugins: 'nib'
```

- **compress** - Compress generated CSS (default: `false`)
- **sourcemaps**
  - **comment** - Adds a comment with the `sourceMappingURL` to the generated CSS (default: `true`)
  - **inline** - Inlines the sourcemap with full source text in base64 format (default: `false`)
  - **sourceRoot** - `sourceRoot` property of the generated sourcemap
  - **basePath** - Base path from which sourcemap and all sources are relative (default: `.`)
- **plugins** - Stylus plugin(s) (default: `nib`)

## Setting Stylus variables

It is possible to set variables that can be used in Stylus.
The purpose of setting variable is to avoid direct modification of the Sylus code,
and thus to make themes more generic

For example, instead of hardcoding:
```stylus
div
 color #FFCC44
```

You can refer to a variable:
```stylus
div
 color convert(hexo-config("moody_red"))
```

And in your **theme's** configuration, you can define this variable:
```yml
moody_red: "#8B0001"
```

(The "convert" function above is here to convert the string into an actual stylus color)

You can also use the theme_config variable in the main `_config.yml`:
```yml
theme_config:
  moody_red: "#8B0001"
```

[Stylus]: http://stylus-lang.com/
[nib]: http://stylus.github.io/nib/

## Extensibility

This plugin provide a filter `stylus:renderer` to allows you extend it. When there’s something you cannot do in Stylus, define it in JavaScript!

For example, to define some global variable:

```js
hexo.extend.filter.register('stylus:renderer', function(style) {
  style
    // we may define a global variable by passing a `Node`
    .define('has-canvas', require('stylus').nodes.false);
    // stylus also casts JavaScript values to their Stylus equivalents when possible
    .define('families', ['Helvetica Neue', 'Helvetica', 'sans-serif'])
    // also allows you to provide a JavaScript-defined function to Stylus
    .define('get-list', function(){
      return ['foo', 'bar', 'baz'];
    });
})
```

Save the file in "scripts/" folder and run Hexo as usual.

Notice: for more JavaScript api, refer to stylus's [documentation](http://stylus-lang.com/docs/js.html).
