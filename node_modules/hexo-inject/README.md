# hexo-inject [![Build Status](https://travis-ci.org/akfish/hexo-inject.svg?branch=master)](https://travis-ci.org/akfish/hexo-inject)

Dynamic script & style (and more) injection for Hexo

## Usage

This plugin is **for plugin/theme developers** to inject custom code into rendered HTML.

### 0. Overview

Injection is called once per complete HTML page (ones that have both `<head>` and `<body>` section).

There are 4 injection points:

Name | API
---- | ---------
`head_begin`  | `inject.headBegin`
`head_end`    | `inject.headEnd`
`body_begin`  | `inject.bodyBegin`
`body_end`    | `inject.bodyEnd`

```html
<!DOCTYPE html>
<html>
  <head>
    <!-- head_begin -->
    <!-- ... -->
    <!-- head_end -->
  </head>
  <body>
    <!-- body_begin -->
    <!-- ... -->
    <!-- body_end -->
  </body>
</html>
```

### 1. Install

Ask your user to run `npm install --save hexo-inject`.

Or add `postinstall` script to your plugin's `package.json`:

```json
{
  "scripts": {
    "postinstall": "npm install --save hexo-inject"
  }
}
```

### 2. Register to `inject_ready` filter

hexo-inject will execute `inject_ready` filter to pool all installed plugins for injection configuration once hexo's `after_init` is fired.

In your plugin:

```js
hexo.extend.filter.register('inject_ready', (inject) => {
  // Configure injections here
  // Inject raw html at head_begin
  inject.raw('head_begin', 'injected content')
  // Or short hand
  inject.headBegin.raw('injected content')
})
```

### 3. Simple injection helpers

hexo-inject provides a few helpers for simple HTML content injection:

* `tag (injectionPoint, name, attrs, content, endTag, opts)`
* `script (injectionPoint, attrs, content, opts)`
* `style (injectionPoint, attrs, content, opts)`
* `link (injectionPoint, attrs, opts)`

Examples:
```js
inject.link('head_begin', { href: '/foo/bar.css', rel: 'stylesheet' })
inject.headBegin.script({}, 'var foo = 1;', { shouldInject: (src) => determinedBy(src) })
```

Notes:
* `injectionPoint` is omitted if the helper is called from short-hand form (e.g `inject.headBegin`)
* All values in `attrs` and `content` can be a `string`, a `Promise` that returns a `string`, or a `function` that returns a `string` or a `Promise`
* `opts.shouldInject` can be a `boolean` value or a `function` that takes current page's HTML source and returns a `boolean` value. If `shouldInject` returns `false`, the content will not be injected to that page.

### 4. Inject files

hexo-inject also provides `require (injectionPoint, module, opts)` helper for file injection.

The workflow is:
* File path is specified by `module` and is resolved relative to the callsite script's folder
* Once resolved, the file will be renderer by hexo's renderer (determined by file's extension). The extension will be changed to output format's. (i.e. `.swig` -> `.html`)
* The rendered content will be processed by loader (again, determined by file's extension) and the result will be injected
* If `opts.inline == false`, hexo-inject will serve the file and reference it accordingly (i.e. via `<script src="/path/to/served.js"></script>`). Otherwise the content will be injected directly.

Valid `opts` fileds are:
* `inline` - a boolean value
* `src` - custom path for serving the file. Default to `/injected/${module.fileName}${module.ext}`
* `data` - passed to renderer
* `shouldInject`

#### Custom loader

hexo-inject provides loader for `.js` and `.css` by default. If you need to handle other formats, you should implement your own loader:

```js
inject.loader.register('.foo', (content, opts) => {
  return opts.inline
    ? `<Foo src=${opts.src}></Foo>`
    : `<Foo>${content}</Foo>`
})
```

Note that you might need to handle `opts.inline` accordingly and know that `content` will be an empty string if `inline == false`.
