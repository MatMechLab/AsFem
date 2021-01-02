# hexo-renderer-marked

[![Build Status](https://travis-ci.org/hexojs/hexo-renderer-marked.svg?branch=master)](https://travis-ci.org/hexojs/hexo-renderer-marked)
[![NPM version](https://badge.fury.io/js/hexo-renderer-marked.svg)](https://www.npmjs.com/package/hexo-renderer-marked)
[![Coverage Status](https://img.shields.io/coveralls/hexojs/hexo-renderer-marked.svg)](https://coveralls.io/r/hexojs/hexo-renderer-marked?branch=master)
[![NPM Dependencies](https://img.shields.io/librariesio/release/npm/hexo-renderer-marked.svg)](https://libraries.io/npm/hexo-renderer-marked)

Add support for [Markdown]. This plugin uses [marked] as its render engine.

## Installation

``` bash
$ npm install hexo-renderer-marked --save
```

- Hexo 4: >= 2.0
- Hexo 3: >= 0.2
- Hexo 2: 0.1.x

## Options

You can configure this plugin in `_config.yml`.

``` yaml
marked:
  gfm: true
  pedantic: false
  breaks: true
  smartLists: true
  smartypants: true
  quotes: '“”‘’'
  modifyAnchors: 0
  anchorAlias: false
  autolink: true
  mangle: true
  sanitizeUrl: false
  headerIds: true
  lazyload: false
  prependRoot: false
  postAsset: false
  external_link:
    enable: false
    exclude: []
    nofollow: false
  disableNunjucks: false
```

- **gfm** - Enables [GitHub flavored markdown](https://help.github.com/articles/github-flavored-markdown)
- **pedantic** - Conform to obscure parts of `markdown.pl` as much as possible. Don't fix any of the original markdown bugs or poor behavior.
- **breaks** - Enable GFM [line breaks](https://help.github.com/articles/github-flavored-markdown#newlines). This option requires the `gfm` option to be true.
- **smartLists** - Use smarter list behavior than the original markdown.
- **smartypants** - Use "smart" typograhic punctuation for things like quotes and dashes.
- **quotes** - Defines the double and single quotes used for substituting regular quotes if **smartypants** is enabled.
  * Example: '«»“”'
    * "double" will be turned into «single»
    * 'single' will be turned into “single”
  * Both double and single quotes substitution must be specified, otherwise it will be silently ignored.
- **modifyAnchors** - Transform the anchorIds into lower case (`1`) or upper case (`2`).
- **autolink** - Enable autolink for URLs. E.g. `https://hexo.io` will become `<a href="https://hexo.io">https://hexo.io</a>`.
- **mangle** - Escape autolinked email address with HTML character references.
  * This is to obscure email address from _basic_ crawler used by spam bot, while still readable to web browsers.
- **sanitizeUrl** - Remove URLs that start with `javascript:`, `vbscript:` and `data:`.
- **headerIds** - Insert header id, e.g. `<h1 id="value">text</h1>`. Useful for inserting anchor link to each paragraph with a heading.
- **anchorAlias** - Enables custom header id
  * Example: `## [foo](#bar)`, id will be set as "bar".
  * Requires **headerIds** to be enabled.
- **lazyload** - Lazy loading images via `loading="lazy"` attribute.
- **prependRoot** - Prepend root value to (internal) image path.
  * Example `_config.yml`:
  ``` yml
  root: /blog/
  ```
  * `![text](/path/to/image.jpg)` becomes `<img src="/blog/path/to/image.jpg" alt="text">`
- **postAsset** - Resolve post asset's image path to relative path and prepend root value when [`post_asset_folder`](https://hexo.io/docs/asset-folders) is enabled.
  * "image.jpg" is located at "/2020/01/02/foo/image.jpg", which is a post asset of "/2020/01/02/foo/".
  * `![](image.jpg)` becomes `<img src="/2020/01/02/foo/image.jpg">`
  * Requires **prependRoot** to be enabled.
- **external_link**
  * **enable** - Open external links in a new tab.
  * **exclude** - Exclude hostname. Specify subdomain when applicable, including `www`.
    - Example: `[foo](http://bar.com)` becomes `<a href="http://bar.com" target="_blank" rel="noopener">foo</a>`
  * **nofollow** - Add `rel="noopener external nofollow noreferrer"` to all external links for security, privacy and SEO. [Read more](https://developer.mozilla.org/en-US/docs/Web/HTML/Link_types). _This can be enabled regardless of `external_link.enable`_
    - Example: `[foo](http://bar.com)` becomes `<a href="http://bar.com" rel="noopener external nofollow noreferrer">foo</a>`
- **disableNunjucks**: If true, Nunjucks tags `{{ }}` or `{% %}` (usually used by [tag plugins](https://hexo.io/docs/tag-plugins)) will not be rendered.

For more options, see [Marked](https://marked.js.org/using_advanced#options). Due to the customizations implemented by this plugin, some of the Marked's options may not work as expected. Feel free to raise an [issue](https://github.com/hexojs/hexo-renderer-marked/issues) to us for clarification.

## Extras

### Definition/Description Lists

`hexo-renderer-marked` also implements description/definition lists using the same syntax as [PHP Markdown Extra][PHP Markdown Extra].

This Markdown:

```markdown
Definition Term
:    This is the definition for the term
```

will generate this HTML:

```html
<dl>
  <dt>Definition Term</dt>
  <dd>This is the definition for the term</dd>
</dl>
```

Note: There is currently a limitation in this implementation. If multiple definitions are provided, the rendered HTML will be incorrect.

For example, this Markdown:

```markdown
Definition Term
:    Definition 1
:    Definition 2
```

will generate this HTML:

```html
<dl>
  <dt>Definition Term<br>: Definition 1</dt>
  <dd>Definition 2</dd>
</dl>
```

If you've got ideas on how to support multiple definitions, please provide a pull request. We'd love to support it.

### Extensibility

This plugin overrides some default behaviours of how [marked] plugin renders the markdown into html, to integrate with the Hexo ecosystem. It is possible to override this plugin too, without resorting to forking the whole thing.

For example, to override how heading like `# heading text` is rendered:

``` js
hexo.extend.filter.register('marked:renderer', function(renderer) {
  const { config } = this; // Skip this line if you don't need user config from _config.yml
  renderer.heading = function(text, level) {
    // Default behaviour
    // return `<h${level}>${text}</h${level}>`;
    // outputs <h1>heading text</h1>

    // If you want to insert custom class name
    return `<h${level} class="headerlink">${text}</h${level}>`;
    // outputs <h1 class="headerlink">heading text</h1>
  }
})
```

Save the file in "scripts/" folder and run Hexo as usual.

Notice `renderer.heading = function (text, level) {` corresponds to [this line](https://github.com/hexojs/hexo-renderer-marked/blob/a93ebeb1e8cc11e754630c0a1506da9a1489b2b0/lib/renderer.js#L21). Refer to [renderer.js](https://github.com/hexojs/hexo-renderer-marked/blob/master/lib/renderer.js) on how this plugin overrides the default methods. For other methods not covered by this plugin, refer to marked's [documentation](https://marked.js.org/using_pro#renderer).

#### Tokenizer

It is also possible to customize the [tokenizer](https://marked.js.org/using_pro#tokenizer).

``` js
const { escape } = require('marked/src/helpers');

// https://github.com/markedjs/marked/blob/b6773fca412c339e0cedd56b63f9fa1583cfd372/src/Lexer.js#L8-L24
// Replace dashes only
const smartypants = (str) => {
  return str
    // em-dashes
    .replace(/---/g, '\u2014')
    // en-dashes
    .replace(/--/g, '\u2013')
};

hexo.extend.filter.register('marked:tokenizer', function(tokenizer) {
  const { smartypants: isSmarty } = this.config.marked;
  tokenizer.inlineText = function(src, inRawBlock) {
    const { rules } = this;

    // https://github.com/markedjs/marked/blob/b6773fca412c339e0cedd56b63f9fa1583cfd372/src/Tokenizer.js#L643-L658
    const cap = rules.inline.text.exec(src);
    if (cap) {
      let text;
      if (inRawBlock) {
        text = cap[0];
      } else {
        text = escape(isSmarty ? smartypants(cap[0], quotes) : cap[0]);
      }
      return {
        // `type` value is a corresponding renderer method
        // https://marked.js.org/using_pro#inline-level-renderer-methods
        type: 'text',
        raw: cap[0],
        text
      };
    }
  }
});
```

[Markdown]: https://daringfireball.net/projects/markdown/
[marked]: https://github.com/chjj/marked
[PHP Markdown Extra]: https://michelf.ca/projects/php-markdown/extra/#def-list
