'use strict';

class Mutex {
  constructor() {
    this._locked = false;
    this._queue = [];
  }

  lock(fn) {
    if (this._locked) {
      this._queue.push(fn);
      return;
    }

    this._locked = true;
    fn();
  }

  unlock() {
    if (!this._locked) return;

    const next = this._queue.shift();

    if (next) {
      next();
    } else {
      this._locked = false;
    }
  }
}

module.exports = Mutex;
