'use strict';

const Promise = require('bluebird');
const fs = require('graceful-fs');
const { dirname, join, extname, basename } = require('path');
const fsPromises = fs.promises;
const chokidar = require('chokidar');
const { escapeRegExp } = require('hexo-util');

const rEOL = /\r\n/g;

function exists(path) {
  if (!path) throw new TypeError('path is required!');

  const promise = fsPromises.access(path).then(() => true, err => {
    if (err.code !== 'ENOENT') throw err;
    return false;
  });

  return Promise.resolve(promise);
}

function existsSync(path) {
  if (!path) throw new TypeError('path is required!');

  try {
    fs.accessSync(path);
  } catch (err) {
    if (err.code !== 'ENOENT') throw err;
    return false;
  }

  return true;
}

function mkdirs(path) {
  if (!path) throw new TypeError('path is required!');

  return Promise.resolve(fsPromises.mkdir(path, { recursive: true }));
}

function mkdirsSync(path) {
  if (!path) throw new TypeError('path is required!');

  fs.mkdirSync(path, { recursive: true });
}

function checkParent(path) {
  return Promise.resolve(fsPromises.mkdir(dirname(path), { recursive: true }));
}

function writeFile(path, data, options = {}) {
  if (!path) throw new TypeError('path is required!');

  if (!data) data = '';

  return checkParent(path).then(() => fsPromises.writeFile(path, data, options));
}

function writeFileSync(path, data, options) {
  if (!path) throw new TypeError('path is required!');

  fs.mkdirSync(dirname(path), { recursive: true });
  fs.writeFileSync(path, data, options);
}

function appendFile(path, data, options = {}) {
  if (!path) throw new TypeError('path is required!');

  return checkParent(path).then(() => fsPromises.appendFile(path, data, options));
}

function appendFileSync(path, data, options) {
  if (!path) throw new TypeError('path is required!');

  fs.mkdirSync(dirname(path), { recursive: true });
  fs.appendFileSync(path, data, options);
}

function copyFile(src, dest, flags) {
  if (!src) throw new TypeError('src is required!');
  if (!dest) throw new TypeError('dest is required!');

  return checkParent(dest).then(() => fsPromises.copyFile(src, dest, flags));
}

function trueFn() {
  return true;
}

function ignoreHiddenFiles(ignore) {
  if (!ignore) return trueFn;

  return ({ name }) => !name.startsWith('.');
}

function ignoreFilesRegex(regex) {
  if (!regex) return trueFn;

  return ({ name }) => !regex.test(name);
}

function ignoreExcludeFiles(arr, parent) {
  if (!arr || !arr.length) return trueFn;

  const set = new Set(arr);

  return ({ name }) => !set.has(join(parent, name));
}

async function _readAndFilterDir(path, options) {
  const { ignoreHidden = true, ignorePattern } = options;
  return (await fsPromises.readdir(path, { ...options, withFileTypes: true }))
    .filter(ignoreHiddenFiles(ignoreHidden))
    .filter(ignoreFilesRegex(ignorePattern));
}

function _readAndFilterDirSync(path, options) {
  const { ignoreHidden = true, ignorePattern } = options;
  return fs.readdirSync(path, { ...options, withFileTypes: true })
    .filter(ignoreHiddenFiles(ignoreHidden))
    .filter(ignoreFilesRegex(ignorePattern));
}

function _copyDirWalker(src, dest, results, parent, options) {
  return Promise.map(_readAndFilterDir(src, options), item => {
    const childSrc = join(src, item.name);
    const childDest = join(dest, item.name);
    const currentPath = join(parent, item.name);

    if (item.isDirectory()) {
      return _copyDirWalker(childSrc, childDest, results, currentPath, options);
    }
    results.push(currentPath);
    return copyFile(childSrc, childDest);
  });
}

function copyDir(src, dest, options = {}) {
  if (!src) throw new TypeError('src is required!');
  if (!dest) throw new TypeError('dest is required!');

  const results = [];

  return checkParent(dest).then(() => _copyDirWalker(src, dest, results, '', options)).return(results);
}

async function _listDirWalker(path, results, parent, options) {
  const promises = [];

  for (const item of await _readAndFilterDir(path, options)) {
    const currentPath = join(parent, item.name);

    if (item.isDirectory()) {
      promises.push(_listDirWalker(join(path, item.name), results, currentPath, options));
    } else {
      results.push(currentPath);
    }
  }

  await Promise.all(promises);
}

function listDir(path, options = {}) {
  if (!path) throw new TypeError('path is required!');

  const results = [];

  return Promise.resolve(_listDirWalker(path, results, '', options)).return(results);
}

function _listDirSyncWalker(path, results, parent, options) {
  for (const item of _readAndFilterDirSync(path, options)) {
    const currentPath = join(parent, item.name);

    if (item.isDirectory()) {
      _listDirSyncWalker(join(path, item.name), results, currentPath, options);
    } else {
      results.push(currentPath);
    }
  }
}

function listDirSync(path, options = {}) {
  if (!path) throw new TypeError('path is required!');

  const results = [];

  _listDirSyncWalker(path, results, '', options);

  return results;
}

function escapeEOL(str) {
  return str.replace(rEOL, '\n');
}

function escapeBOM(str) {
  return str.charCodeAt(0) === 0xFEFF ? str.substring(1) : str;
}

function escapeFileContent(content) {
  return escapeBOM(escapeEOL(content));
}

async function _readFile(path, options) {
  if (!Object.prototype.hasOwnProperty.call(options, 'encoding')) options.encoding = 'utf8';

  const content = await fsPromises.readFile(path, options);

  if (options.escape == null || options.escape) {
    return escapeFileContent(content);
  }

  return content;
}

function readFile(path, options = {}) {
  if (!path) throw new TypeError('path is required!');

  return Promise.resolve(_readFile(path, options));
}

function readFileSync(path, options = {}) {
  if (!path) throw new TypeError('path is required!');

  if (!Object.prototype.hasOwnProperty.call(options, 'encoding')) options.encoding = 'utf8';

  const content = fs.readFileSync(path, options);

  if (options.escape == null || options.escape) {
    return escapeFileContent(content);
  }

  return content;
}

async function _emptyDir(path, parent, options) {
  const entrys = (await _readAndFilterDir(path, options))
    .filter(ignoreExcludeFiles(options.exclude, parent));
  const results = [];

  await Promise.map(entrys, item => {
    const fullPath = join(path, item.name);
    const currentPath = join(parent, item.name);

    if (item.isDirectory()) {
      return _emptyDir(fullPath, currentPath, options).then(files => {
        if (!files.length) {
          return fsPromises.rmdir(fullPath);
        }
        results.push(...files);
      });
    }
    results.push(currentPath);
    return fsPromises.unlink(fullPath);
  });

  return results;
}

function emptyDir(path, options = {}) {
  if (!path) throw new TypeError('path is required!');

  return Promise.resolve(_emptyDir(path, '', options));
}

function _emptyDirSync(path, options, parent) {
  const entrys = _readAndFilterDirSync(path, options)
    .filter(ignoreExcludeFiles(options.exclude, parent));

  const results = [];

  for (const item of entrys) {
    const childPath = join(path, item.name);
    const currentPath = join(parent, item.name);

    if (item.isDirectory()) {
      const removed = _emptyDirSync(childPath, options, currentPath);

      if (!fs.readdirSync(childPath).length) {
        rmdirSync(childPath);
      }

      results.push(...removed);
    } else {
      fs.unlinkSync(childPath);
      results.push(currentPath);
    }
  }

  return results;
}

function emptyDirSync(path, options = {}) {
  if (!path) throw new TypeError('path is required!');

  return _emptyDirSync(path, options, '');
}

async function _rmdir(path) {
  const files = fsPromises.readdir(path, { withFileTypes: true });
  await Promise.map(files, item => {
    const childPath = join(path, item.name);

    return item.isDirectory() ? _rmdir(childPath) : fsPromises.unlink(childPath);
  });
  return fsPromises.rmdir(path);
}

function rmdir(path) {
  if (!path) throw new TypeError('path is required!');

  return Promise.resolve(_rmdir(path));
}

function _rmdirSync(path) {
  const files = fs.readdirSync(path, { withFileTypes: true });

  for (const item of files) {
    const childPath = join(path, item.name);

    if (item.isDirectory()) {
      _rmdirSync(childPath);
    } else {
      fs.unlinkSync(childPath);
    }
  }

  fs.rmdirSync(path);
}

function rmdirSync(path) {
  if (!path) throw new TypeError('path is required!');

  _rmdirSync(path);
}

function watch(path, options = {}) {
  if (!path) throw new TypeError('path is required!');

  const watcher = chokidar.watch(path, options);

  return new Promise((resolve, reject) => {
    watcher.on('ready', resolve);
    watcher.on('error', reject);
  }).thenReturn(watcher);
}

function _findUnusedPath(path, files) {
  const ext = extname(path);
  const base = basename(path, ext);
  const regex = new RegExp(`^${escapeRegExp(base)}(?:-(\\d+))?${escapeRegExp(ext)}$`);
  let num = -1;

  for (let i = 0, len = files.length; i < len; i++) {
    const item = files[i];
    const match = item.match(regex);

    if (match == null) continue;
    const matchNum = match[1] ? parseInt(match[1], 10) : 0;

    if (matchNum > num) {
      num = matchNum;
    }
  }

  return join(dirname(path), `${base}-${num + 1}${ext}`);
}

async function _ensurePath(path) {
  if (!await exists(path)) return path;

  const files = await fsPromises.readdir(dirname(path));
  return _findUnusedPath(path, files);
}

function ensurePath(path) {
  if (!path) throw new TypeError('path is required!');

  return Promise.resolve(_ensurePath(path));
}

function ensurePathSync(path) {
  if (!path) throw new TypeError('path is required!');
  if (!fs.existsSync(path)) return path;

  const files = fs.readdirSync(dirname(path));

  return _findUnusedPath(path, files);
}

function ensureWriteStream(path, options = {}) {
  if (!path) throw new TypeError('path is required!');

  return checkParent(path).then(() => fs.createWriteStream(path, options));
}

function ensureWriteStreamSync(path, options) {
  if (!path) throw new TypeError('path is required!');

  fs.mkdirSync(dirname(path), { recursive: true });
  return fs.createWriteStream(path, options);
}

// access
['F_OK', 'R_OK', 'W_OK', 'X_OK'].forEach(key => {
  Object.defineProperty(exports, key, {
    enumerable: true,
    value: fs.constants[key],
    writable: false
  });
});

exports.access = Promise.promisify(fs.access);
exports.accessSync = fs.accessSync;

// appendFile
exports.appendFile = appendFile;
exports.appendFileSync = appendFileSync;

// chmod
exports.chmod = Promise.promisify(fs.chmod);
exports.chmodSync = fs.chmodSync;
exports.fchmod = Promise.promisify(fs.fchmod);
exports.fchmodSync = fs.fchmodSync;
exports.lchmod = Promise.promisify(fs.lchmod);
exports.lchmodSync = fs.lchmodSync;

// chown
exports.chown = Promise.promisify(fs.chown);
exports.chownSync = fs.chownSync;
exports.fchown = Promise.promisify(fs.fchown);
exports.fchownSync = fs.fchownSync;
exports.lchown = Promise.promisify(fs.lchown);
exports.lchownSync = fs.lchownSync;

// close
exports.close = Promise.promisify(fs.close);
exports.closeSync = fs.closeSync;

// copy
exports.copyDir = copyDir;
exports.copyFile = copyFile;

// createStream
exports.createReadStream = fs.createReadStream;
exports.createWriteStream = fs.createWriteStream;

// emptyDir
exports.emptyDir = emptyDir;
exports.emptyDirSync = emptyDirSync;

// ensurePath
exports.ensurePath = ensurePath;
exports.ensurePathSync = ensurePathSync;

// ensureWriteStream
exports.ensureWriteStream = ensureWriteStream;
exports.ensureWriteStreamSync = ensureWriteStreamSync;

// exists
exports.exists = exists;
exports.existsSync = existsSync;

// fsync
exports.fsync = Promise.promisify(fs.fsync);
exports.fsyncSync = fs.fsyncSync;

// link
exports.link = Promise.promisify(fs.link);
exports.linkSync = fs.linkSync;

// listDir
exports.listDir = listDir;
exports.listDirSync = listDirSync;

// mkdir
exports.mkdir = Promise.promisify(fs.mkdir);
exports.mkdirSync = fs.mkdirSync;

// mkdirs
exports.mkdirs = mkdirs;
exports.mkdirsSync = mkdirsSync;

// open
exports.open = Promise.promisify(fs.open);
exports.openSync = fs.openSync;

// symlink
exports.symlink = Promise.promisify(fs.symlink);
exports.symlinkSync = fs.symlinkSync;

// read
exports.read = Promise.promisify(fs.read);
exports.readSync = fs.readSync;

// readdir
exports.readdir = Promise.promisify(fs.readdir);
exports.readdirSync = fs.readdirSync;

// readFile
exports.readFile = readFile;
exports.readFileSync = readFileSync;

// readlink
exports.readlink = Promise.promisify(fs.readlink);
exports.readlinkSync = fs.readlinkSync;

// realpath
exports.realpath = Promise.promisify(fs.realpath);
exports.realpathSync = fs.realpathSync;

// rename
exports.rename = Promise.promisify(fs.rename);
exports.renameSync = fs.renameSync;

// rmdir
exports.rmdir = rmdir;
exports.rmdirSync = rmdirSync;

// stat
exports.stat = Promise.promisify(fs.stat);
exports.statSync = fs.statSync;
exports.fstat = Promise.promisify(fs.fstat);
exports.fstatSync = fs.fstatSync;
exports.lstat = Promise.promisify(fs.lstat);
exports.lstatSync = fs.lstatSync;

// truncate
exports.truncate = Promise.promisify(fs.truncate);
exports.truncateSync = fs.truncateSync;
exports.ftruncate = Promise.promisify(fs.ftruncate);
exports.ftruncateSync = fs.ftruncateSync;

// unlink
exports.unlink = Promise.promisify(fs.unlink);
exports.unlinkSync = fs.unlinkSync;

// utimes
exports.utimes = Promise.promisify(fs.utimes);
exports.utimesSync = fs.utimesSync;
exports.futimes = Promise.promisify(fs.futimes);
exports.futimesSync = fs.futimesSync;

// watch
exports.watch = watch;
exports.watchFile = fs.watchFile;
exports.unwatchFile = fs.unwatchFile;

// write
exports.write = Promise.promisify(fs.write);
exports.writeSync = fs.writeSync;

// writeFile
exports.writeFile = writeFile;
exports.writeFileSync = writeFileSync;

// Static classes
exports.Stats = fs.Stats;
exports.ReadStream = fs.ReadStream;
exports.WriteStream = fs.WriteStream;
exports.FileReadStream = fs.FileReadStream;
exports.FileWriteStream = fs.FileWriteStream;

// util
exports.escapeBOM = escapeBOM;
exports.escapeEOL = escapeEOL;
