// Convert reviewed Markdown notes to jemdoc wrappers with build-time KaTeX.
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { createHash } from 'node:crypto';
import { Marked } from './tools/markdown/marked.esm.js';
import katex from './tools/jemdoc-cvx-1.0.0/node_modules/katex/dist/katex.mjs';

const root = path.dirname(fileURLToPath(import.meta.url));
const site = path.join(root, 'site');
const directory = path.join(site, 'notes');
const escape = value => String(value).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const cssVersion = file => createHash('sha256').update(fs.readFileSync(path.join(site, file))).digest('hex').slice(0, 12);

if (fs.existsSync(directory)) for (const file of fs.readdirSync(directory).filter(f => f.endsWith('.md'))) {
  let source = fs.readFileSync(path.join(directory, file), 'utf8');
  const title = source.match(/^# (.+)$/m)?.[1];
  if (!title) throw new Error(`Missing title in ${file}`);
  source = source.replace(/^# .+\r?\n/, '');
  // Stable per-note IDs keep existing deep links valid when sections move.
  // Each heading maps to an ID or [primary ID, ...legacy aliases].
  const anchorFile = path.join(directory, file.replace(/\.md$/, '.anchors.json'));
  const hasAnchorMap = fs.existsSync(anchorFile);
  const anchorMap = hasAnchorMap ? JSON.parse(fs.readFileSync(anchorFile, 'utf8')) : null;
  if (hasAnchorMap) {
    if (anchorMap === null || typeof anchorMap !== 'object' || Array.isArray(anchorMap)) {
      throw new Error(`Invalid heading anchor map in ${anchorFile}`);
    }
    const configuredIds = new Set();
    for (const [heading, value] of Object.entries(anchorMap)) {
      const ids = Array.isArray(value) ? value : [value];
      if (!ids.length || ids.some(id => typeof id !== 'string' || !/^[A-Za-z][A-Za-z0-9_-]*$/.test(id))) {
        throw new Error(`Invalid anchor ID for heading "${heading}" in ${anchorFile}`);
      }
      for (const id of ids) {
        if (configuredIds.has(id)) throw new Error(`Duplicate anchor ID "${id}" in ${anchorFile}`);
        configuredIds.add(id);
      }
    }
  }
  const equations = [];
  // Protect math before Markdown table parsing so TeX pipes and escapes survive.
  source = source.replace(/(`[^`\n]*`)|\$\$([\s\S]*?)\$\$|\$([^$\n]+?)\$/g, (raw, code, display, inline) => {
    if (code) return code;
    const id = `NOTEMATHPLACEHOLDER${equations.length}END`;
    const isDisplay = display !== undefined;
    const tex = (isDisplay ? display : inline).trim();
    try {
      equations.push({id, isDisplay, html: katex.renderToString(tex, {displayMode:isDisplay, throwOnError:true, trust:false, strict:'error'})});
    } catch (error) {
      throw new Error(`${file}, equation ${equations.length + 1}: ${tex}\n${error.message}`);
    }
    return isDisplay ? `\n\n${id}\n\n` : id;
  });
  const headings = [];
  const usedIds = new Set();
  const renderer = {
    heading(token) {
      if (anchorMap !== null && !Object.hasOwn(anchorMap, token.text)) {
        throw new Error(`Missing stable anchor for heading "${token.text}" in ${file}`);
      }
      const mapped = anchorMap === null ? `section-${headings.length + 1}` : anchorMap[token.text];
      const ids = Array.isArray(mapped) ? mapped : [mapped];
      for (const id of ids) {
        if (usedIds.has(id)) throw new Error(`Duplicate rendered anchor ID "${id}" in ${file}`);
        usedIds.add(id);
      }
      const [id, ...aliases] = ids;
      const label = this.parser.parseInline(token.tokens);
      headings.push({id, depth:token.depth, label});
      const legacy = aliases.map(alias => `<span id="${escape(alias)}"></span>\n`).join('');
      return legacy + `<h${token.depth} id="${escape(id)}">${label}</h${token.depth}>\n`;
    },
    link(token) {
      if (!/^(https?:\/\/|#)/.test(token.href)) throw new Error(`Unsupported note link: ${token.href}`);
      const target = token.href.startsWith('#') ? '' : ' target="_blank" rel="noopener noreferrer"';
      return `<a href="${escape(token.href)}"${target}>${this.parser.parseInline(token.tokens)}</a>`;
    },
    image(token) {
      if (!token.href.startsWith('../images/notes/')) throw new Error(`Unexpected note image: ${token.href}`);
      const src = escape(token.href);
      const data = fs.readFileSync(path.resolve(directory, token.href));
      // Reserve the PNG's aspect ratio so lazy loading does not move section links.
      const dimensions = data.subarray(0, 8).equals(Buffer.from([137,80,78,71,13,10,26,10]))
        ? ` width="${data.readUInt32BE(16)}" height="${data.readUInt32BE(20)}"` : '';
      return `<a class="note-image-link" href="${src}" target="_blank" rel="noopener noreferrer"><img src="${src}" alt="${escape(token.text)}"${dimensions} loading="lazy"></a>`;
    },
    html(token) { return escape(token.text); },
  };
  const parser = new Marked({gfm:true, breaks:false});
  parser.use({renderer});
  let html = parser.parse(source);
  html = html.replace(/<table>/g, '<div class="note-table-scroll"><table>').replace(/<\/table>/g, '</table></div>');
  const toc = '<details class="note-toc"><summary>目录</summary><ul>' + headings.filter(h => h.depth === 2).map(h => `<li><a href="#${h.id}">${h.label}</a></li>`).join('') + '</ul></details>';
  html = toc + '\n<article class="note-article">\n' + html + '</article>\n';
  for (const eq of equations) {
    if (eq.isDisplay) html = html.replace(`<p>${eq.id}</p>`, eq.html);
    html = html.replaceAll(eq.id, eq.html);
  }
  if (html.includes('NOTEMATHPLACEHOLDER')) throw new Error(`Unresolved math in ${file}`);
  const wrapper = [
    `# jemdoc: menu{../MENU}{notes.html}{../}, addcss{../custom.css?v=${cssVersion('custom.css')}}, addcss{../notes.css?v=${cssVersion('notes.css')}}, addcss{../katex/katex.min.css}`,
    '# Generated by generate_notes.mjs; edit the corresponding Markdown file.',
    '= '+title,
    '', '~~~', '{}{raw}', html.trimEnd(), '~~~', '',
  ].join('\n');
  fs.writeFileSync(path.join(directory, file.replace(/\.md$/, '.jemdoc')), wrapper);
  console.log(`Rendered ${file}: ${equations.length} equations, ${headings.length} headings.`);
}
