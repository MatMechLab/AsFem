'use strict';

Object.defineProperty(exports, '__esModule', {
  value: true
});
exports['default'] = mixin;

var _privateUtils = require('./private/utils');

var defineProperty = Object.defineProperty;

function handleClass(target, mixins) {
  if (!mixins.length) {
    throw new SyntaxError('@mixin() class ' + target.name + ' requires at least one mixin as an argument');
  }

  for (var i = 0, l = mixins.length; i < l; i++) {
    var descs = (0, _privateUtils.getOwnPropertyDescriptors)(mixins[i]);

    for (var key in descs) {
      if (!(key in target.prototype)) {
        defineProperty(target.prototype, key, descs[key]);
      }
    }
  }
}

function mixin() {
  for (var _len = arguments.length, mixins = Array(_len), _key = 0; _key < _len; _key++) {
    mixins[_key] = arguments[_key];
  }

  if (typeof mixins[0] === 'function') {
    return handleClass(mixins[0], []);
  } else {
    return function (target) {
      return handleClass(target, mixins);
    };
  }
}

module.exports = exports['default'];