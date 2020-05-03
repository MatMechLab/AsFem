# markdown-it-toc

> markdown-it plugin for adding a table of contents to markdown documents

## Usage

#### Enable plugin

```js
var md = require('markdown-it')({
  html: true,
  linkify: true,
  typography: true
}).use(require('markdown-it-toc')); // <-- this use(package_name) is required
```

#### Example

```md
@[toc](Title)
```

Adding this tag with add anchors to each ```<h[n]>``` tag on your document, and will add a ```<ul>``` of hyperlinks pointing to these places on the page.

The end results looks like:

```html
<p>
     <h3>Title</h3>
     <ul>
	<li><a href="...">Heading 1</a></li>
	...
	... 
     </ul> 
</p>
...
...
<h1><a href="..."></a>Heading 1</h1>
```

### Testing

To run the tests use:
```bash
make test
```