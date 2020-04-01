# Warehouse

[![Build Status](https://travis-ci.com/hexojs/warehouse.svg?branch=master)](https://travis-ci.com/hexojs/warehouse)
[![NPM version](https://badge.fury.io/js/warehouse.svg)](https://www.npmjs.com/package/warehouse)
[![Coverage Status](https://coveralls.io/repos/github/hexojs/warehouse/badge.svg?branch=master)](https://coveralls.io/github/hexojs/warehouse?branch=master)

A JSON database with Models, Schemas, and a flexible querying interface. It powers the wildly successful static site generator [Hexo](https://hexo.io).

## Installation

``` bash
$ npm install warehouse
```

## 3.0 BREAKING CHANGE

In `warehouse@3`, the constructor has been changed from function declaration to `class` declaration or definition by `class` expression.
Derived classes of classes defined by `class` declarations and `class` expressions must also be defined in `class` declaration, `class` expression.
Anyone who created their own SchemaType will need to change.

``` js
const SchemaType = require('warehouse/schematype');

class MySchemaType extends SchemaType {
  constructor(name, options = {}) {
    super(name, Object.assign({ foo: 'foo' }, options));
  }
}
```

It changes to a class declaration or a class expression, but it does not need to deal with other than the definition of the constructor.

``` js
// It work!

MySchemaType.prototype.cast = function (value, data) {
  let result = SchemaType.prototype.cast.call(this, value, data);
  return result ? result : '';
}
```

## Usage

``` js
var Database = require('warehouse');
var db = new Database();

var Post = db.model('posts', {
  title: String,
  created: {type: Date, default: Date.now}
});

Post.insert({
  title: 'Hello world'
}).then(function(post){
  console.log(post);
});
```

+ [API documentation](https://hexojs.github.io/warehouse/)
+ [Examples of `Model`, `Schema`, and `SchemaType`](https://github.com/hexojs/hexo/tree/master/lib/models)
+ [More examples in `./test/scripts`](./test/scripts)

## Test

``` bash
$ npm test
```