'use strict';

Object.defineProperty(exports, '__esModule', {
  value: true
});
exports['default'] = decorate;

function _toArray(arr) { return Array.isArray(arr) ? arr : Array.from(arr); }

function _toConsumableArray(arr) { if (Array.isArray(arr)) { for (var i = 0, arr2 = Array(arr.length); i < arr.length; i++) arr2[i] = arr[i]; return arr2; } else { return Array.from(arr); } }

var _privateUtils = require('./private/utils');

function handleDescriptor(target, key, _ref, _ref2) {
  var fn = _ref.value;

  var _ref22 = _toArray(_ref2);

  var decorator = _ref22[0];

  var args = _ref22.slice(1);

  return {
    configurable: true,
    enumerable: false,
    value: decorator.apply(undefined, [fn].concat(_toConsumableArray(args)))
  };
}

function decorate() {
  for (var _len = arguments.length, args = Array(_len), _key = 0; _key < _len; _key++) {
    args[_key] = arguments[_key];
  }

  return (0, _privateUtils.decorate)(handleDescriptor, args);
}

module.exports = exports['default'];