'use strict';

Object.defineProperty(exports, "__esModule", {
  value: true
});

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

var Router = function () {
  function Router(hexo) {
    _classCallCheck(this, Router);

    this.hexo = hexo;
    this._routes = [];
  }

  _createClass(Router, [{
    key: 'register',
    value: function register() {
      var _this = this;

      var generator = this.hexo.extend.generator;

      generator.register('inject', function (locals) {
        return _this._routes;
      });
    }
  }, {
    key: 'serve',
    value: function serve(module, opts) {
      var src = opts.src || '/injected/' + module.name + module.ext;
      var content = module.content;
      this._routes.push({ path: src, data: function data() {
          return content;
        } });
      return src;
    }
  }]);

  return Router;
}();

exports.default = Router;
module.exports = exports['default'];
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbInJvdXRlci5lczYiXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6Ijs7Ozs7Ozs7OztJQUFxQjtBQUNuQixXQURtQixNQUNuQixDQUFhLElBQWIsRUFBbUI7MEJBREEsUUFDQTs7QUFDakIsU0FBSyxJQUFMLEdBQVksSUFBWixDQURpQjtBQUVqQixTQUFLLE9BQUwsR0FBZSxFQUFmLENBRmlCO0dBQW5COztlQURtQjs7K0JBS1A7OztVQUNKLFlBQWMsS0FBSyxJQUFMLENBQVUsTUFBVixDQUFkLFVBREk7O0FBRVYsZ0JBQVUsUUFBVixDQUFtQixRQUFuQixFQUE2QixVQUFDLE1BQUQ7ZUFBWSxNQUFLLE9BQUw7T0FBWixDQUE3QixDQUZVOzs7OzBCQUlMLFFBQVEsTUFBTTtBQUNuQixVQUFJLE1BQU0sS0FBSyxHQUFMLG1CQUF5QixPQUFPLElBQVAsR0FBYyxPQUFPLEdBQVAsQ0FEOUI7QUFFbkIsVUFBSSxVQUFVLE9BQU8sT0FBUCxDQUZLO0FBR25CLFdBQUssT0FBTCxDQUFhLElBQWIsQ0FBa0IsRUFBRSxNQUFNLEdBQU4sRUFBVyxNQUFNO2lCQUFNO1NBQU4sRUFBckMsRUFIbUI7QUFJbkIsYUFBTyxHQUFQLENBSm1COzs7O1NBVEYiLCJmaWxlIjoicm91dGVyLmpzIiwic291cmNlc0NvbnRlbnQiOlsiZXhwb3J0IGRlZmF1bHQgY2xhc3MgUm91dGVyIHtcclxuICBjb25zdHJ1Y3RvciAoaGV4bykge1xyXG4gICAgdGhpcy5oZXhvID0gaGV4b1xyXG4gICAgdGhpcy5fcm91dGVzID0gW11cclxuICB9XHJcbiAgcmVnaXN0ZXIgKCkge1xyXG4gICAgbGV0IHsgZ2VuZXJhdG9yIH0gPSB0aGlzLmhleG8uZXh0ZW5kXHJcbiAgICBnZW5lcmF0b3IucmVnaXN0ZXIoJ2luamVjdCcsIChsb2NhbHMpID0+IHRoaXMuX3JvdXRlcylcclxuICB9XHJcbiAgc2VydmUgKG1vZHVsZSwgb3B0cykge1xyXG4gICAgbGV0IHNyYyA9IG9wdHMuc3JjIHx8IGAvaW5qZWN0ZWQvJHttb2R1bGUubmFtZX0ke21vZHVsZS5leHR9YFxyXG4gICAgbGV0IGNvbnRlbnQgPSBtb2R1bGUuY29udGVudFxyXG4gICAgdGhpcy5fcm91dGVzLnB1c2goeyBwYXRoOiBzcmMsIGRhdGE6ICgpID0+IGNvbnRlbnQgfSlcclxuICAgIHJldHVybiBzcmNcclxuICB9XHJcbn1cclxuIl0sInNvdXJjZVJvb3QiOiIvc291cmNlLyJ9
