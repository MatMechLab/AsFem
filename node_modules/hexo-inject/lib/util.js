'use strict';

Object.defineProperty(exports, "__esModule", {
  value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }();

exports.camalize = camalize;
exports.callsite = callsite;

var _const = require('./const');

var _path = require('path');

var _path2 = _interopRequireDefault(_path);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function camalize(str) {
  return str.split('_').filter(function (s) {
    return s.length > 0;
  }).map(function (s, i) {
    return i === 0 ? s : s[0].toUpperCase() + s.substr(1);
  }).join('');
}

function callsite() {
  function parse(t) {
    var _REGEX$stack_trace$ex = _const.REGEX.stack_trace.exec(t);

    var _REGEX$stack_trace$ex2 = _slicedToArray(_REGEX$stack_trace$ex, 6);

    var functionName = _REGEX$stack_trace$ex2[1];
    var alias = _REGEX$stack_trace$ex2[2];
    var filePath = _REGEX$stack_trace$ex2[3];
    var line = _REGEX$stack_trace$ex2[4];
    var col = _REGEX$stack_trace$ex2[5];

    var file = _path2.default.parse(filePath);

    line = parseInt(line);
    col = parseInt(col);

    return {
      functionName: functionName, alias: alias,
      filePath: filePath, file: file,
      line: line, col: col
    };
  }
  var stack = new Error().stack.split('\n').slice(2) // First line is 'Error', second line is this function
  .map(parse);

  return stack;
}
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbInV0aWwuZXM2Il0sIm5hbWVzIjpbXSwibWFwcGluZ3MiOiI7Ozs7Ozs7O1FBR2dCO1FBT0E7O0FBVmhCOztBQUNBOzs7Ozs7QUFFTyxTQUFTLFFBQVQsQ0FBbUIsR0FBbkIsRUFBd0I7QUFDN0IsU0FBTyxJQUFJLEtBQUosQ0FBVSxHQUFWLEVBQ0osTUFESSxDQUNHLFVBQUMsQ0FBRDtXQUFPLEVBQUUsTUFBRixHQUFXLENBQVg7R0FBUCxDQURILENBRUosR0FGSSxDQUVBLFVBQUMsQ0FBRCxFQUFJLENBQUo7V0FBVSxNQUFNLENBQU4sR0FBVSxDQUFWLEdBQWUsRUFBRSxDQUFGLEVBQUssV0FBTCxLQUFxQixFQUFFLE1BQUYsQ0FBUyxDQUFULENBQXJCO0dBQXpCLENBRkEsQ0FHSixJQUhJLENBR0MsRUFIRCxDQUFQLENBRDZCO0NBQXhCOztBQU9BLFNBQVMsUUFBVCxHQUFxQjtBQUMxQixXQUFTLEtBQVQsQ0FBZ0IsQ0FBaEIsRUFBbUI7Z0NBQ2tDLGFBQU0sV0FBTixDQUFrQixJQUFsQixDQUF1QixDQUF2QixFQURsQzs7OztRQUNWLHlDQURVO1FBQ0ksa0NBREo7UUFDVyxxQ0FEWDtRQUNxQixpQ0FEckI7UUFDMkIsZ0NBRDNCOztBQUVqQixRQUFJLE9BQU8sZUFBSyxLQUFMLENBQVcsUUFBWCxDQUFQLENBRmE7O0FBSWpCLFdBQU8sU0FBUyxJQUFULENBQVAsQ0FKaUI7QUFLakIsVUFBTSxTQUFTLEdBQVQsQ0FBTixDQUxpQjs7QUFPakIsV0FBTztBQUNMLGdDQURLLEVBQ1MsWUFEVDtBQUVMLHdCQUZLLEVBRUssVUFGTDtBQUdMLGdCQUhLLEVBR0MsUUFIRDtLQUFQLENBUGlCO0dBQW5CO0FBYUEsTUFBSSxRQUFRLElBQUksS0FBSixHQUFZLEtBQVosQ0FDVCxLQURTLENBQ0gsSUFERyxFQUNHLEtBREgsQ0FDUyxDQURUO0dBRVQsR0FGUyxDQUVMLEtBRkssQ0FBUixDQWRzQjs7QUFrQjFCLFNBQU8sS0FBUCxDQWxCMEI7Q0FBckIiLCJmaWxlIjoidXRpbC5qcyIsInNvdXJjZXNDb250ZW50IjpbImltcG9ydCB7IElOSkVDVElPTl9QT0lOVFMsIFJFR0VYIH0gZnJvbSAnLi9jb25zdCdcclxuaW1wb3J0IHBhdGggZnJvbSAncGF0aCdcclxuXHJcbmV4cG9ydCBmdW5jdGlvbiBjYW1hbGl6ZSAoc3RyKSB7XHJcbiAgcmV0dXJuIHN0ci5zcGxpdCgnXycpXHJcbiAgICAuZmlsdGVyKChzKSA9PiBzLmxlbmd0aCA+IDApXHJcbiAgICAubWFwKChzLCBpKSA9PiBpID09PSAwID8gcyA6IChzWzBdLnRvVXBwZXJDYXNlKCkgKyBzLnN1YnN0cigxKSkpXHJcbiAgICAuam9pbignJylcclxufVxyXG5cclxuZXhwb3J0IGZ1bmN0aW9uIGNhbGxzaXRlICgpIHtcclxuICBmdW5jdGlvbiBwYXJzZSAodCkge1xyXG4gICAgbGV0IFssIGZ1bmN0aW9uTmFtZSwgYWxpYXMsIGZpbGVQYXRoLCBsaW5lLCBjb2xdID0gUkVHRVguc3RhY2tfdHJhY2UuZXhlYyh0KVxyXG4gICAgbGV0IGZpbGUgPSBwYXRoLnBhcnNlKGZpbGVQYXRoKVxyXG5cclxuICAgIGxpbmUgPSBwYXJzZUludChsaW5lKVxyXG4gICAgY29sID0gcGFyc2VJbnQoY29sKVxyXG5cclxuICAgIHJldHVybiB7XHJcbiAgICAgIGZ1bmN0aW9uTmFtZSwgYWxpYXMsXHJcbiAgICAgIGZpbGVQYXRoLCBmaWxlLFxyXG4gICAgICBsaW5lLCBjb2xcclxuICAgIH1cclxuICB9XHJcbiAgbGV0IHN0YWNrID0gbmV3IEVycm9yKCkuc3RhY2tcclxuICAgIC5zcGxpdCgnXFxuJykuc2xpY2UoMikgLy8gRmlyc3QgbGluZSBpcyAnRXJyb3InLCBzZWNvbmQgbGluZSBpcyB0aGlzIGZ1bmN0aW9uXHJcbiAgICAubWFwKHBhcnNlKVxyXG5cclxuICByZXR1cm4gc3RhY2tcclxufVxyXG4iXSwic291cmNlUm9vdCI6Ii9zb3VyY2UvIn0=
