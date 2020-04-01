'use strict';

Object.defineProperty(exports, "__esModule", {
  value: true
});
var INJECTION_POINTS = exports.INJECTION_POINTS = ['head_begin', 'head_end', 'body_begin', 'body_end'];

var REGEX = exports.REGEX = {
  head_begin: /(<head.*>[\n\r\s\t]*)/i,
  head_end: /([\n\r\s\t]*<\/head>)/i,
  body_begin: /(<body.*>[\n\r\s\t]*)/i,
  body_end: /([\n\r\s\t]*<\/body>)/i,
  injection_begin: /(<!-- hexo-inject:begin -->)/i,
  injection_end: /(<!-- hexo-inject:end -->)/i,
  stack_trace: /\s+at(?:\s(\S*))?(?:\s\[as\s(\S*)\])?\s\(?(\S*?):(\d+):(\d+)\)?/
};

var API = exports.API = ['raw', 'tag', 'script', 'style', 'link', 'require'];
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbImNvbnN0LmVzNiJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiOzs7OztBQUFPLElBQU0sOENBQW1CLENBQzlCLFlBRDhCLEVBRTlCLFVBRjhCLEVBRzlCLFlBSDhCLEVBSTlCLFVBSjhCLENBQW5COztBQU9OLElBQU0sd0JBQVE7QUFDbkIsY0FBb0Isd0JBQXBCO0FBQ0EsWUFBb0Isd0JBQXBCO0FBQ0EsY0FBb0Isd0JBQXBCO0FBQ0EsWUFBb0Isd0JBQXBCO0FBQ0EsbUJBQW9CLCtCQUFwQjtBQUNBLGlCQUFvQiw2QkFBcEI7QUFDQSxlQUFvQixpRUFBcEI7Q0FQVzs7QUFVTixJQUFNLG9CQUFNLENBQ2pCLEtBRGlCLEVBRWpCLEtBRmlCLEVBR2pCLFFBSGlCLEVBSWpCLE9BSmlCLEVBS2pCLE1BTGlCLEVBTWpCLFNBTmlCLENBQU4iLCJmaWxlIjoiY29uc3QuanMiLCJzb3VyY2VzQ29udGVudCI6WyJleHBvcnQgY29uc3QgSU5KRUNUSU9OX1BPSU5UUyA9IFtcclxuICAnaGVhZF9iZWdpbicsXHJcbiAgJ2hlYWRfZW5kJyxcclxuICAnYm9keV9iZWdpbicsXHJcbiAgJ2JvZHlfZW5kJ1xyXG5dXHJcblxyXG5leHBvcnQgY29uc3QgUkVHRVggPSB7XHJcbiAgaGVhZF9iZWdpbiAgICAgICAgOiAvKDxoZWFkLio+W1xcblxcclxcc1xcdF0qKS9pLFxyXG4gIGhlYWRfZW5kICAgICAgICAgIDogLyhbXFxuXFxyXFxzXFx0XSo8XFwvaGVhZD4pL2ksXHJcbiAgYm9keV9iZWdpbiAgICAgICAgOiAvKDxib2R5Lio+W1xcblxcclxcc1xcdF0qKS9pLFxyXG4gIGJvZHlfZW5kICAgICAgICAgIDogLyhbXFxuXFxyXFxzXFx0XSo8XFwvYm9keT4pL2ksXHJcbiAgaW5qZWN0aW9uX2JlZ2luICAgOiAvKDwhLS0gaGV4by1pbmplY3Q6YmVnaW4gLS0+KS9pLFxyXG4gIGluamVjdGlvbl9lbmQgICAgIDogLyg8IS0tIGhleG8taW5qZWN0OmVuZCAtLT4pL2ksXHJcbiAgc3RhY2tfdHJhY2UgICAgICAgOiAvXFxzK2F0KD86XFxzKFxcUyopKT8oPzpcXHNcXFthc1xccyhcXFMqKVxcXSk/XFxzXFwoPyhcXFMqPyk6KFxcZCspOihcXGQrKVxcKT8vXHJcbn1cclxuXHJcbmV4cG9ydCBjb25zdCBBUEkgPSBbXHJcbiAgJ3JhdycsXHJcbiAgJ3RhZycsXHJcbiAgJ3NjcmlwdCcsXHJcbiAgJ3N0eWxlJyxcclxuICAnbGluaycsXHJcbiAgJ3JlcXVpcmUnXHJcbl1cclxuIl0sInNvdXJjZVJvb3QiOiIvc291cmNlLyJ9
