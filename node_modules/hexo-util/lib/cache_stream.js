'use strict';

const { Transform } = require('stream');

class CacheStream extends Transform {
  constructor() {
    super();

    this._cache = [];
  }

  _transform(chunk, enc, callback) {
    const buf = chunk instanceof Buffer ? chunk : Buffer.from(chunk, enc);
    this._cache.push(buf);
    this.push(buf);
    callback();
  }

  getCache() {
    return Buffer.concat(this._cache);
  }
}

module.exports = CacheStream;
