'use strict';

const crypto = require('crypto');

const ALGORITHM = 'sha1';

function createSha1Hash() {
  return crypto.createHash(ALGORITHM);
}

exports.hash = content => {
  const hash = createSha1Hash();
  hash.update(content);
  return hash.digest();
};

exports.createSha1Hash = createSha1Hash;
