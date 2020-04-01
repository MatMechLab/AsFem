# hexo-generator-index

[![Build Status](https://travis-ci.org/hexojs/hexo-generator-index.svg?branch=master)](https://travis-ci.org/hexojs/hexo-generator-index)
[![NPM version](https://badge.fury.io/js/hexo-generator-index.svg)](https://www.npmjs.com/package/hexo-generator-index)
[![Coverage Status](https://img.shields.io/coveralls/hexojs/hexo-generator-index.svg)](https://coveralls.io/r/hexojs/hexo-generator-index?branch=master)

Index generator for [Hexo].

## Installation

``` bash
$ npm install hexo-generator-index --save
```

## Options
Add or modify the following section to your root _config.yml file

``` yaml
index_generator:
  path: ''
  per_page: 10
  order_by: -date
```

- **path**: Root path for your blog's index page. 
  - default: ""
- **per_page**: Posts displayed per page.
  - default: [`config.per_page`](https://hexo.io/docs/configuration.html#Pagination) as specified in the official Hexo docs (if present), otherwise `10`
  - `0` disables pagination
- **order_by**: Posts order. 
  - default: date descending
- **pagination_dir**: URL format.
  - default: 'page'
  - `awesome-page` makes the URL ends with 'awesome-page/<page number>' for second page and beyond.

## License

MIT

[Hexo]: http://hexo.io/
