'use strict';

const { vsprintf } = require('sprintf-js');

function i18n(options = {}) {
  this.data = {};
  this.languages = options.languages || ['default'];

  if (!Array.isArray(this.languages)) {
    this.languages = [this.languages];
  }
}

i18n.prototype.get = function(languages) {
  const { data } = this;
  const result = {};

  if (languages) {
    if (!Array.isArray(languages)) {
      languages = [languages];
    }
  } else {
    languages = this.languages;
  }

  for (let i = 0, leni = languages.length; i < leni; i++) {
    const lang = languages[i];
    const langData = data[lang];
    if (!langData) continue;

    const keys = Object.keys(langData);

    for (let j = 0, lenj = keys.length; j < lenj; j++) {
      const key = keys[j];
      if (!Object.prototype.hasOwnProperty.call(result, key)) result[key] = langData[key];
    }
  }

  return result;
};

i18n.prototype.set = function(lang, data) {
  if (typeof lang !== 'string') throw new TypeError('lang must be a string!');
  if (typeof data !== 'object') throw new TypeError('data is required!');

  this.data[lang] = flattenObject(data);

  return this;
};

i18n.prototype.remove = function(lang) {
  if (typeof lang !== 'string') throw new TypeError('lang must be a string!');

  delete this.data[lang];

  return this;
};

i18n.prototype.list = function() {
  return Object.keys(this.data);
};

function flattenObject(data, obj = {}, parent = '') {
  const keys = Object.keys(data);

  for (let i = 0, len = keys.length; i < len; i++) {
    const key = keys[i];
    const item = data[key];

    if (typeof item === 'object') {
      flattenObject(item, obj, parent + key + '.');
    } else {
      obj[parent + key] = item;
    }
  }

  return obj;
}

i18n.prototype.__ = function(lang) {
  const data = this.get(lang);

  return (key, ...args) => {
    if (!key) return '';

    const str = data[key] || key;

    return vsprintf(str, args);
  };
};

i18n.prototype._p = function(lang) {
  const data = this.get(lang);

  return (key, ...args) => {
    if (!key) return '';

    const number = args.length ? +args[0] : 0;
    let str = key;

    if (!number && Object.prototype.hasOwnProperty.call(data, key + '.zero')) {
      str = data[key + '.zero'];
    } else if (number === 1 && Object.prototype.hasOwnProperty.call(data, key + '.one')) {
      str = data[key + '.one'];
    } else if (Object.prototype.hasOwnProperty.call(data, key + '.other')) {
      str = data[key + '.other'];
    } else if (Object.prototype.hasOwnProperty.call(data, key)) {
      str = data[key];
    }

    return vsprintf(str, args);
  };
};

module.exports = i18n;
