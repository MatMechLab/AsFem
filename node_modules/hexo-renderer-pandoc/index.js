var spawnSync = require('child_process').spawnSync;

var pandocRenderer = function(data, options){
  var config = hexo.config.pandoc;
  var extensions = '', filters = [], extra = [];
  // To satisfy pandoc's requirement that html5 must have a title.
  // Since the markdown file is only rendered as body part,
  // the title is never used and thus does not matter
  var meta = ['-M', 'pagetitle=dummy'];
  var math = '--mathjax';

  if(config) {
    if(config.extensions) {
      config.extensions.forEach(function(extension) {
        extensions += extension;
      });
    }

    if(config.filters) {
      config.filters.forEach(function(filter) {
        filters.push('--filter');
        filters.push(filter);
      });
    }

    if(config.extra) {
      for(var e in config.extra) {
        var eoption = config.extra[e];
        for (var key in eoption){
          extra.push('--' + key);
          if(eoption[key]!=null) {
            extra.push(eoption[key]); 
          }
        }
      }
    }

    if(config.meta) {
      config.meta.forEach(function(m) {
        meta.push('-M');
        if(m.length) {
          meta.push(m);
        } else {
          for(var m2 in m) {
            meta.push(m2 + '=' + m[m2]);
          }
        }
      });
    }

    if(config.mathEngine) {
      if(typeof config.mathEngine === 'string') {
        math = '--' + config.mathEngine;
      }
    }
  }

  var args = [ '-f', 'markdown-smart'+extensions, '-t', 'html-smart', math]
  .concat(filters)
  .concat(extra)
  .concat(meta);


  // if we are rendering a post,
  // `data` has the key `path`
  // https://github.com/hexojs/hexo/blob/2ed17cd105768df379dad8bbbe4df30964fe8f2d/lib/hexo/post.js#L269
  // otherwise (e.g., rendering a tag),
  // `path` is not present in `data`.
  // https://github.com/hexojs/hexo/blob/2ed17cd105768df379dad8bbbe4df30964fe8f2d/lib/extend/tag.js#L173
  // https://github.com/hexojs/hexo/blob/a6dc0ea28dddad1b5f1bad7c6f86f1e0627b564a/lib/plugins/tag/blockquote.js#L64

  // are we rendering a standalone post?
  if("path" in data) {
    // only apply template when rendering post, not tags
    if (config && config.template) {
      args.push("--template=" + config.template);
    }

    // do not apply `--standalone`,
    // header/footer are to be added by Hexo

    // also set a metavariable to let concerned
    // pandoc filters know
    args.push(...["-M", "standalone=True"]);
  }
  // or some thing to be embedded in a post,
  // like tags?
  else {
    args.push(...["-M", "standalone=False"]);
  }

  var src = data.text.toString();

  var res = spawnSync('pandoc', args, {
    cwd: process.cwd(),
    env: process.env,
    encoding: "utf8",
    input: src
  });

  if (res.status == 0) {
    if (res.stderr) {
      var warn_msg = ''
        + '[WARNING][hexo-renderer-pandoc] On ' + data.path + '\n'
        + '[WARNING][hexo-renderer-pandoc] ' + res.stderr;
      console.log(warn_msg);
    }

    return res.stdout;
  } else {
    var error_msg = '\n'
      + '[ERROR][hexo-renderer-pandoc] On ' + data.path + '\n'
      + '[ERROR][hexo-renderer-pandoc] pandoc exited with code '+res.status+(res.stderr ? ': ' + res.stderr : '.');
    throw Error(error_msg);
    console.log(error_msg);
    return null;
  }
}

hexo.extend.renderer.register('md', 'html', pandocRenderer, true);
hexo.extend.renderer.register('markdown', 'html', pandocRenderer, true);
hexo.extend.renderer.register('mkd', 'html', pandocRenderer, true);
hexo.extend.renderer.register('mkdn', 'html', pandocRenderer, true);
hexo.extend.renderer.register('mdwn', 'html', pandocRenderer, true);
hexo.extend.renderer.register('mdtxt', 'html', pandocRenderer, true);
hexo.extend.renderer.register('mdtext', 'html', pandocRenderer, true);
