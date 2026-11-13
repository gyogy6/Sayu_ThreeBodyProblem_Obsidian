// Preview server for developing the scene without Obsidian: node dev/serve.js
const http = require('http');
const fs = require('fs');
const path = require('path');

const root = path.join(__dirname, '..');
const types = { '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css' };

http
  .createServer((req, res) => {
    const url = req.url.split('?')[0];
    const rel = url === '/' ? 'dev/harness.html' : url.replace(/^\/+/, '');
    const file = path.join(root, rel);
    if (!file.startsWith(root) || !fs.existsSync(file) || fs.statSync(file).isDirectory()) {
      res.writeHead(404);
      res.end('not found');
      return;
    }
    res.writeHead(200, { 'Content-Type': (types[path.extname(file)] || 'text/plain') + '; charset=utf-8' });
    fs.createReadStream(file).pipe(res);
  })
  .listen(5177, '127.0.0.1', () => console.log('http://localhost:5177'));
