'use strict';

const JSONStream = require('JSONStream');
const Promise = require('bluebird');
const fs = require('graceful-fs');
const Model = require('./model');
const Schema = require('./schema');
const SchemaType = require('./schematype');
const WarehouseError = require('./error');
const pkg = require('../package.json');

const uncorkAsync = stream => new Promise(resolve => {
  process.nextTick(() => {
    stream.uncork();
    resolve();
  });
});

const exportAsync = Promise.coroutine(function* (database, path) {
  const stream = fs.createWriteStream(path);

  // Start body & Meta & Start models
  stream.write(`{"meta":${JSON.stringify({
    version: database.options.version,
    warehouse: pkg.version
  })},"models":{`);

  const models = database._models;
  const keys = Object.keys(models);
  const { length } = keys;

  // support stream backpressure
  const writeAsync = Promise.promisify(stream.write, { context: stream });

  // models body
  for (let i = 0; i < length; i++) {
    const key = keys[i];

    if (!models[key]) continue;

    stream.cork();
    if (i) stream.write(',', 'ascii');

    stream.write(`"${key}":`);

    const modelExportPromise = writeAsync(models[key]._export());

    yield uncorkAsync(stream);
    yield modelExportPromise;
  }

  // End models
  stream.end('}}', 'ascii');
  return new Promise((resolve, reject) => {
    stream.on('error', reject)
      .on('finish', resolve);
  });
});

class Database {

  /**
   * Database constructor.
   *
   * @param {object} [options]
   *   @param {number} [options.version=0] Database version
   *   @param {string} [options.path] Database path
   *   @param {function} [options.onUpgrade] Triggered when the database is upgraded
   *   @param {function} [options.onDowngrade] Triggered when the database is downgraded
   */
  constructor(options) {
    this.options = Object.assign({
      version: 0,
      onUpgrade() {},

      onDowngrade() {}
    }, options);

    this._models = {};

    class _Model extends Model {}

    this.Model = _Model;

    _Model.prototype._database = this;
  }

  /**
   * Creates a new model.
   *
   * @param {string} name
   * @param {Schema|object} [schema]
   * @return {Model}
   */
  model(name, schema) {
    if (this._models[name]) {
      return this._models[name];
    }

    this._models[name] = new this.Model(name, schema);
    const model = this._models[name];
    return model;
  }

  /**
   * Loads database.
   *
   * @param {function} [callback]
   * @return {Promise}
   */
  load(callback) {
    const { path, onUpgrade, onDowngrade, version: newVersion } = this.options;

    if (!path) throw new WarehouseError('options.path is required');

    let oldVersion = 0;

    const getMetaCallBack = data => {
      if (data.meta && data.meta.version) {
        oldVersion = data.meta.version;
      }
    };

    // data event arg0 wrap key/value pair.
    const parseStream = JSONStream.parse('models.$*');

    parseStream.once('header', getMetaCallBack);
    parseStream.once('footer', getMetaCallBack);

    parseStream.on('data', data => {
      this.model(data.key)._import(data.value);
    });

    const rs = fs.createReadStream(path, 'utf8');

    const promise = new Promise((resolve, reject) => {
      parseStream.once('error', reject);
      parseStream.once('end', resolve);

      rs.once('error', reject);
    }).then(() => {
      if (newVersion > oldVersion) {
        return onUpgrade(oldVersion, newVersion);
      } else if (newVersion < oldVersion) {
        return onDowngrade(oldVersion, newVersion);
      }
    }).asCallback(callback);

    rs.pipe(parseStream);

    return promise;
  }

  /**
   * Saves database.
   *
   * @param {function} [callback]
   * @return {Promise}
   */
  save(callback) {
    const { path } = this.options;

    if (!path) throw new WarehouseError('options.path is required');
    return exportAsync(this, path).asCallback(callback);
  }

  toJSON() {
    const models = Object.keys(this._models)
      .reduce((obj, key) => {
        const value = this._models[key];
        if (value != null) obj[key] = value;
        return obj;
      }, {});

    return {
      meta: {
        version: this.options.version,
        warehouse: pkg.version
      }, models
    };
  }
}

Database.prototype.Schema = Schema;
Database.Schema = Database.prototype.Schema;
Database.prototype.SchemaType = SchemaType;
Database.SchemaType = Database.prototype.SchemaType;
Database.version = pkg.version;

module.exports = Database;
