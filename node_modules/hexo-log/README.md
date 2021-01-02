# hexo-log

[![Build Status](https://travis-ci.org/hexojs/hexo-log.svg?branch=master)](https://travis-ci.org/hexojs/hexo-log)
[![NPM version](https://badge.fury.io/js/hexo-log.svg)](https://www.npmjs.com/package/hexo-log)
[![Coverage Status](https://coveralls.io/repos/hexojs/hexo-log/badge.svg?branch=master)](https://coveralls.io/r/hexojs/hexo-log?branch=master)

Logger for Hexo.

## Installation

``` bash
$ npm install hexo-log --save
```

## Usage

``` js
const log = require('hexo-log')({
  debug: false,
  silent: false
})

log.info('Hello world');
```

Option | Description | Default
--- | --- | ---
`debug` | Display debug message. | `false`
`silent` | Don't display any message in console. | `false`

## License

MIT
