'use strict';

const yaml = require('js-yaml');
const rPrefixSep = /^(-{3,}|;{3,})/;
const rFrontMatter = /^(-{3,}|;{3,})\n([\s\S]+?)\n\1\n?([\s\S]*)/;
const rFrontMatterNew = /^([\s\S]+?)\n(-{3,}|;{3,})\n?([\s\S]*)/;

function split(str) {
  if (typeof str !== 'string') throw new TypeError('str is required!');

  const matchOld = str.match(rFrontMatter);
  if (matchOld) {
    return {
      data: matchOld[2],
      content: matchOld[3] || '',
      separator: matchOld[1],
      prefixSeparator: true
    };
  }

  if (rPrefixSep.test(str)) return { content: str };

  const matchNew = str.match(rFrontMatterNew);

  if (matchNew) {
    return {
      data: matchNew[1],
      content: matchNew[3] || '',
      separator: matchNew[2],
      prefixSeparator: false
    };
  }

  return { content: str };
}

function parse(str, options) {
  if (typeof str !== 'string') throw new TypeError('str is required!');

  const splitData = split(str);
  const raw = splitData.data;

  if (!raw) return { _content: str };

  let data;

  if (splitData.separator.startsWith(';')) {
    data = parseJSON(raw);
  } else {
    data = parseYAML(raw, options);
  }

  if (!data) return { _content: str };

  // Convert timezone
  Object.keys(data).forEach(key => {
    const item = data[key];

    if (item instanceof Date) {
      data[key] = new Date(item.getTime() + (item.getTimezoneOffset() * 60 * 1000));
    }
  });

  data._content = splitData.content;
  return data;
}

function parseYAML(str, options) {
  const result = yaml.load(escapeYAML(str), options);
  if (typeof result !== 'object') return;

  return result;
}

function parseJSON(str) {
  try {
    return JSON.parse(`{${str}}`);
  } catch (err) {
    return; // eslint-disable-line
  }
}

function escapeYAML(str) {
  if (typeof str !== 'string') throw new TypeError('str is required!');

  return str.replace(/\n(\t+)/g, (match, tabs) => {
    let result = '\n';

    for (let i = 0, len = tabs.length; i < len; i++) {
      result += '  ';
    }

    return result;
  });
}

function stringify(obj, options = {}) {
  if (!obj) throw new TypeError('obj is required!');

  const { _content: content = '' } = obj;
  delete obj._content;

  if (!Object.keys(obj).length) return content;

  const { mode, prefixSeparator } = options;
  const separator = options.separator || (mode === 'json' ? ';;;' : '---');
  let result = '';

  if (prefixSeparator) result += `${separator}\n`;

  if (mode === 'json') {
    result += stringifyJSON(obj);
  } else {
    result += stringifyYAML(obj, options);
  }

  result += `${separator}\n${content}`;

  return result;
}

function stringifyYAML(obj, options) {
  const keys = Object.keys(obj);
  const data = {};
  const nullKeys = [];
  const dateKeys = [];
  let key, value, i, len;

  for (i = 0, len = keys.length; i < len; i++) {
    key = keys[i];
    value = obj[key];

    if (value == null) {
      nullKeys.push(key);
    } else if (value instanceof Date) {
      dateKeys.push(key);
    } else {
      data[key] = value;
    }
  }

  let result = yaml.dump(data, options);

  if (dateKeys.length) {
    for (i = 0, len = dateKeys.length; i < len; i++) {
      key = dateKeys[i];
      result += `${key}: ${formatDate(obj[key])}\n`;
    }
  }

  if (nullKeys.length) {
    for (i = 0, len = nullKeys.length; i < len; i++) {
      result += `${nullKeys[i]}:\n`;
    }
  }

  return result;
}

function stringifyJSON(obj) {
  return JSON.stringify(obj, null, '  ')
    // Remove indention
    .replace(/\n {2}/g, () => '\n')
    // Remove prefixing and trailing braces
    .replace(/^{\n|}$/g, '');
}

function doubleDigit(num) {
  return num.toString().padStart(2, '0');
}

function formatDate(date) {
  return `${date.getFullYear()}-${doubleDigit(date.getMonth() + 1)}-${doubleDigit(date.getDate())} ${doubleDigit(date.getHours())}:${doubleDigit(date.getMinutes())}:${doubleDigit(date.getSeconds())}`;
}

module.exports = {
  parse,
  split,
  escape: escapeYAML,
  stringify
};
