// test_pages.mjs -- offline integration test of the front-end against real API
// payloads captured in webui/testdata/*.json.
//
// It installs a minimal DOM/Canvas shim, imports every page module, mounts it and
// drives the interactive bits (sorting, filtering, tabs, the 3D renderer) so that
// undefined fields, wrong imports and render-time crashes are caught without a
// browser. It does NOT verify pixels -- open the console for that.
//
// Usage:
//   node webui/test_pages.mjs            (after: python webui/capture_fixtures.py)

import { readFileSync, existsSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const HERE = dirname(fileURLToPath(import.meta.url));
const TESTDATA = join(HERE, 'testdata');

/* ------------------------------------------------------------------ DOM shim */
class ClassList {
  constructor(node) { this.node = node; this.set = new Set(); }
  add(...names) { for (const n of names) if (n) this.set.add(n); }
  remove(...names) { for (const n of names) this.set.delete(n); }
  toggle(name, force) {
    if (force === undefined) {
      if (this.set.has(name)) { this.set.delete(name); return false; }
      this.set.add(name); return true;
    }
    if (force) this.set.add(name); else this.set.delete(name);
    return !!force;
  }
  contains(name) { return this.set.has(name); }
}

class FakeNode {
  constructor(tag, ns = null) {
    this.tagName = String(tag).toUpperCase();
    this.ns = ns;
    this.children = [];
    this.parentNode = null;
    this.attrs = {};
    this.style = {};
    this.dataset = {};
    this.classList = new ClassList(this);
    this._text = '';
    this._html = '';
    this._listeners = {};
    this.value = '';
    this.checked = false;
    this.disabled = false;
    this.scrollTop = 0;
    this.scrollHeight = 0;
  }
  appendChild(child) { this.children.push(child); child.parentNode = this; return child; }
  append(...kids) { for (const kid of kids) if (kid) this.appendChild(kid); }
  prepend(...kids) { for (const kid of kids.reverse()) if (kid) { this.children.unshift(kid); kid.parentNode = this; } }
  replaceChildren(...kids) { this.children = []; this._text = ''; this._html = ''; for (const kid of kids) if (kid) this.appendChild(kid); }
  removeChild(child) { this.children = this.children.filter((c) => c !== child); }
  remove() { this.parentNode?.removeChild(this); }
  setAttribute(name, value) { this.attrs[name] = value; if (name === 'class') this.className = value; }
  getAttribute(name) { return this.attrs[name]; }
  removeAttribute(name) { delete this.attrs[name]; }
  addEventListener(type, fn) { (this._listeners[type] = this._listeners[type] || []).push(fn); }
  removeEventListener(type, fn) {
    if (this._listeners[type]) this._listeners[type] = this._listeners[type].filter((f) => f !== fn);
  }
  dispatch(type, event = {}) {
    for (const fn of this._listeners[type] || []) fn(Object.assign({ target: this, preventDefault() {} }, event));
  }
  click() { this.dispatch('click'); }
  get className() { return [...this.classList.set].join(' '); }
  set className(value) {
    this.classList.set = new Set(String(value || '').split(/\s+/).filter(Boolean));
  }
  get id() { return this.attrs.id || ''; }
  set id(value) { this.attrs.id = value; }
  get href() { return this.attrs.href || ''; }
  set href(value) { this.attrs.href = value; }
  get textContent() {
    if (this.children.length) return this.children.map((c) => c.textContent).join('');
    if (!this._text && this._html) return this._html.replace(/<[^>]*>/g, '');
    return this._text;
  }
  set textContent(value) { this._text = String(value ?? ''); this.children = []; }
  get innerHTML() { return this._html; }
  set innerHTML(value) { this._html = String(value ?? ''); this.children = []; }
  get firstChild() { return this.children[0] || null; }
  get clientWidth() { return 1000; }
  get clientHeight() { return 620; }
  get offsetWidth() { return 1000; }
  getBoundingClientRect() { return { left: 0, top: 0, right: 1000, bottom: 620, width: 1000, height: 620 }; }
  focus() {}
  blur() {}
  scrollIntoView() {}
  setPointerCapture() {}
  releasePointerCapture() {}
  contains(node) { return node === this || this.children.some((c) => c.contains(node)); }
  querySelector(selector) { return this.querySelectorAll(selector)[0] || null; }
  querySelectorAll(selector) {
    const groups = selector.split(',').map((s) => s.trim()).filter(Boolean);
    const matches = (node, part) => {
      for (const piece of part.split(/\s+/).filter(Boolean)) {
        // only the last piece is matched against the node; ancestors are checked below
        if (piece.startsWith('.') && node.classList.contains(piece.slice(1))) return true;
        if (piece.startsWith('#') && node.id === piece.slice(1)) return true;
        if (!piece.startsWith('.') && !piece.startsWith('#') && node.tagName === piece.toUpperCase()) return true;
      }
      return false;
    };
    const matchesChain = (node, parts) => {
      if (!matches(node, parts[parts.length - 1])) return false;
      let idx = parts.length - 2;
      let cursor = node.parentNode;
      while (idx >= 0 && cursor) {
        if (matches(cursor, parts[idx])) idx -= 1;
        cursor = cursor.parentNode;
      }
      return idx < 0;
    };
    const out = [];
    const walk = (node) => {
      for (const child of node.children) {
        for (const group of groups) {
          const parts = group.split(/\s+/).filter(Boolean);
          if (matchesChain(child, parts)) { out.push(child); break; }
        }
        walk(child);
      }
    };
    walk(this);
    return out;
  }
  getContext(kind) {
    if (kind !== '2d') return null;
    return makeContext(this);
  }
}

function makeContext(canvas) {
  const noop = () => {};
  return {
    canvas,
    setTransform: noop, resetTransform: noop, save: noop, restore: noop,
    translate: noop, rotate: noop, scale: noop, clearRect: noop, fillRect: noop, strokeRect: noop,
    beginPath: noop, closePath: noop, moveTo: noop, lineTo: noop, arc: noop, rect: noop,
    fill: noop, stroke: noop, clip: noop, setLineDash: noop,
    fillText: noop, strokeText: noop,
    measureText: (text) => ({ width: String(text).length * 6 }),
    createRadialGradient: () => ({ addColorStop: noop }),
    createLinearGradient: () => ({ addColorStop: noop }),
    createPattern: () => null,
    drawImage: noop,
  };
}

class FakeCanvas extends FakeNode {
  constructor() { super('canvas'); this.width = 800; this.height = 600; }
  getContext(kind) { return kind === '2d' ? makeContext(this) : null; }
}

const ids = new Map();
const documentShim = {
  documentElement: new FakeNode('html'),
  body: new FakeNode('body'),
  createElement(tag) {
    if (String(tag).toLowerCase() === 'canvas') return new FakeCanvas();
    return new FakeNode(tag);
  },
  createElementNS(ns, tag) { return new FakeNode(tag, ns); },
  createTextNode(text) { const n = new FakeNode('#text'); n._text = String(text); return n; },
  getElementById(id) {
    if (!ids.has(id)) {
      const node = new FakeNode('div');
      node.id = id;
      ids.set(id, node);
      if (id === 'nav' || id === 'view' || id === 'crumb') documentShim.body.appendChild(node);
    }
    return ids.get(id);
  },
  querySelector(sel) { return documentShim.body.querySelector(sel); },
  querySelectorAll(sel) { return documentShim.body.querySelectorAll(sel); },
  addEventListener() {}, removeEventListener() {},
};

class FakePath2D {
  constructor() { this.ops = 0; }
  rect() { this.ops += 1; }
  moveTo() {} lineTo() {} arc() {} closePath() {}
}

const rafQueue = [];
globalThis.Node = FakeNode;
const winListeners = {};
globalThis.document = documentShim;
globalThis.window = {
  devicePixelRatio: 2,
  addEventListener(type, fn) { (winListeners[type] = winListeners[type] || []).push(fn); },
  removeEventListener(type, fn) {
    if (winListeners[type]) winListeners[type] = winListeners[type].filter((f) => f !== fn);
  },
  location: { hash: '#/overview', href: 'http://local/#/overview' },
  history: { replaceState() {} },
  requestAnimationFrame: (fn) => { rafQueue.push(fn); return rafQueue.length; },
  cancelAnimationFrame() {},
  getComputedStyle: () => ({ getPropertyValue: () => '' }),
  matchMedia: () => ({ matches: false, addEventListener() {} }),
};
globalThis.document.defaultView = globalThis.window;
globalThis.location = globalThis.window.location;
globalThis.history = globalThis.window.history;
globalThis.localStorage = {
  store: new Map(),
  getItem(k) { return this.store.get(k) ?? null; },
  setItem(k, v) { this.store.set(k, String(v)); },
  removeItem(k) { this.store.delete(k); },
};
Object.defineProperty(globalThis, 'navigator', {
  value: { clipboard: { writeText: async () => {} }, userAgent: 'node-test' },
  configurable: true, writable: true,
});
globalThis.requestAnimationFrame = globalThis.window.requestAnimationFrame;
globalThis.cancelAnimationFrame = () => {};
globalThis.devicePixelRatio = 2;
globalThis.ResizeObserver = class { observe() {} disconnect() {} unobserve() {} };
globalThis.Path2D = FakePath2D;
globalThis.URL.createObjectURL = () => 'blob:test';
globalThis.URL.revokeObjectURL = () => {};
globalThis.setInterval = () => 0;
globalThis.clearInterval = () => {};
globalThis.alert = () => {};

function flushRaf() {
  let guard = 0;
  while (rafQueue.length && guard < 20) {
    const fn = rafQueue.shift();
    guard += 1;
    try { fn(performance.now()); } catch (err) { throw new Error(`rAF callback failed: ${err.message}`); }
  }
}

/* --------------------------------------------------------------- fetch stub */
const FIXTURES = {
  '/api/meta': 'meta', '/api/overview': 'overview', '/api/hotspots': 'hotspots',
  '/api/pharmacophore': 'pharmacophore', '/api/pocket': 'pocket',
  '/api/scaffolds': 'scaffolds', '/api/filters': 'filters', '/api/poses': 'poses',
  '/api/evidence/report': 'evidence_report', '/api/evidence/priors': 'evidence_priors',
  '/api/evidence/census': 'evidence_census', '/api/evidence/interface': 'evidence_interface',
  '/api/docs': 'docs', '/api/jobs': 'jobs', '/api/pocket/points': 'pocket_points',
};

function fixtureFor(url) {
  const [path, query = ''] = url.split('?');
  if (FIXTURES[path]) return FIXTURES[path];
  if (path === '/api/scene') return query.includes('ref=SglM') ? 'scene_m' : 'scene_cj3';
  if (path.startsWith('/api/docs/')) return 'docs_01';
  if (path.startsWith('/api/jobs/')) return 'jobs';
  return null;
}

globalThis.fetch = async (url) => {
  const name = fixtureFor(url);
  if (!name) return { ok: false, status: 404, statusText: 'not found', text: async () => JSON.stringify({ error: `no fixture for ${url}` }) };
  const file = join(TESTDATA, `${name}.json`);
  if (!existsSync(file)) throw new Error(`missing fixture ${file} -- run: python webui/capture_fixtures.py`);
  const body = readFileSync(file, 'utf8');
  const payload = JSON.parse(body);
  if (url.includes('/api/pocket/points') && payload.status !== 'ready') {
    payload.status = 'ready';
    payload.n_voxels = 0;
    payload.xyz = [];
  }
  return { ok: true, status: 200, statusText: 'OK', text: async () => JSON.stringify(payload) };
};

/* ------------------------------------------------------------------- runner */
const ctx = {
  params: {},
  navigate: (route, params) => console.log(`      navigate -> #/${route}${params ? ` ${JSON.stringify(params)}` : ''}`),
  toast: (msg, type) => console.log(`      toast(${type || 'info'}): ${msg}`),
  meta: null,
  refresh: async () => {},
};

const results = [];
async function runPage(name, file, params = {}, extra = null) {
  const root = new FakeNode('main');
  const localCtx = { ...ctx, params };
  const started = Date.now();
  try {
    const module = await import(`./static/js/pages/${file}`);
    const cleanup = await module.mount(root, localCtx);
    flushRaf();
    await new Promise((resolve) => setTimeout(resolve, 0));
    if (extra) await extra(root, localCtx);
    const nodes = countNodes(root);
    if (!nodes) throw new Error('page rendered no DOM nodes');
    if (typeof cleanup === 'function') cleanup();
    results.push([name, 'ok', `${nodes} nodes, ${Date.now() - started} ms`]);
  } catch (err) {
    results.push([name, 'FAIL', `${err.message}\n${String(err.stack).split('\n').slice(1, 4).join('\n')}`]);
  }
}

function countNodes(node) {
  let n = 1;
  for (const child of node.children || []) n += countNodes(child);
  return n;
}

const tables = (root) => root.querySelectorAll('table');
const textOf = (root) => root.textContent;
const expectText = (root, pattern, what) => {
  if (!pattern.test(textOf(root))) throw new Error(`missing ${what} in rendered text`);
};

await runPage('overview', 'overview.js', {}, async (root) => {
  expectText(root, /20,307/, '枚举总数 20,307');
  expectText(root, /538/, '三重达标 538');
  expectText(root, /13\.0/, '最佳骨架分 13.0');
  expectText(root, /5,440|5440/, '腔体体积');
});
await runPage('hotspots', 'hotspots.js', {}, async (root) => {
  const trs = root.querySelectorAll('tr');
  if (trs.length < 20) throw new Error(`hotspot table looks empty (${trs.length} rows)`);
  const th = root.querySelectorAll('th').find((t) => t.classList.contains('sorted'));
  if (!th) throw new Error('no sorted header found');
  th.click();                          // re-sort by another column
  const first = root.querySelectorAll('tbody tr')[0];
  first.click();                       // open the residue detail panel
  expectText(root, /接触肽残基数|收敛度/, '逐残基详情');
  expectText(root, /dSASA/, 'dSASA 列');
});
await runPage('pharmacophore', 'pharmacophore.js', {}, async (root) => {
  if (!root.querySelectorAll('svg').length) throw new Error('scatter svg missing');
  expectText(root, /A1-H3|A1-R2/, '设计轴编号');
  expectText(root, /99\.0|9.*18/, '设计窗口');
});
await runPage('pocket', 'pocket.js', {}, async (root) => {
  if (root.querySelectorAll('table').length < 1) throw new Error('pocket table missing');
  expectText(root, /S1_basic_bowl/, '亚口袋名');
  expectText(root, /0\.94|0\.73/, '埋藏度');
});
await runPage('scaffolds', 'scaffolds.js', {}, async (root) => {
  const rows = root.querySelectorAll('tbody tr');
  if (rows.length < 25) throw new Error(`scaffold table too small (${rows.length})`);
  expectText(root, /OC\(=O\)c1ccc/, 'SMILES');
  expectText(root, /318/, 'MW 318');
  rows[0].click();                     // opens the drawer
  expectText(root, /位阻|盐桥|氢键/, '打分分解');
  const search = root.querySelectorAll('input').find((i) => i.attrs && i.attrs.type === 'search');
  if (!search) throw new Error('search input missing');
  search.value = 'K368';
  search.dispatch('input');
  const filtered = root.querySelectorAll('tbody tr').length;
  if (filtered >= 200) throw new Error(`search did not filter (${filtered} rows)`);
  search.value = 'ZZZZ';
  search.dispatch('input');
  expectText(root, /没有满足|没有匹配|无匹配|没有数据/, '空结果提示');
  search.value = '';
  search.dispatch('input');
  console.log(`      filter "K368" -> ${filtered} rows, empty state ok`);
});
await runPage('evidence', 'evidence.js', {}, async (root) => {
  const tabs = root.querySelectorAll('.tab');
  if (tabs.length < 5) throw new Error(`expected 5 tabs, got ${tabs.length}`);
  expectText(root, /界面规模|先验/, '报告正文');
  tabs[1].click(); flushRaf();
  expectText(root, /PASS/, '先验核对状态');
  expectText(root, /SglM:Asp18/, 'FLAG 条目');
  tabs[2].click(); flushRaf();
  expectText(root, /496/, 'MurJ 496 残基');
  tabs[3].click(); flushRaf();
  expectText(root, /dSASA/, '界面测量 dSASA');
  tabs[4].click(); flushRaf();
  expectText(root, /K46|R53|K368/, '锚点残基');
});
await runPage('docs', 'docs.js', {}, async (root) => {
  const md = root.querySelectorAll('.md');
  if (!md.length) throw new Error('markdown body missing');
  expectText(root, /WEDGE-Design/, '文档正文');
});
await runPage('run', 'run.js', {}, async (root) => {
  const cards = root.querySelectorAll('.card');
  if (cards.length < 5) throw new Error(`launcher cards missing (${cards.length})`);
  expectText(root, /证据链|设计管线/, '任务卡片');
  expectText(root, /运行历史|控制台/, '日志面板');
});
await runPage('structure', 'structure.js', { ligand: '0' }, async (root) => {
  const canvas = root.querySelectorAll('canvas');
  if (!canvas.length) throw new Error('viewer canvas missing');
  const checkboxes = root.querySelectorAll('input');
  if (checkboxes.length < 6) throw new Error('viewer controls missing');
  expectText(root, /9NU8/, '参考结构标签');
  expectText(root, /OC\(=O\)/, '配体 SMILES');
  for (const box of checkboxes.slice(0, 3)) {
    box.checked = !box.checked;
    box.dispatch('change');
  }
  flushRaf();
});

/* --------------------------------------------------------------- unit checks */
async function unit(name, fn) {
  try { await fn(); results.push([name, 'ok', '']); } catch (err) { results.push([name, 'FAIL', err.message]); }
}

await unit('ui.markdown (tables, lists, code, inline)', async () => {
  const { markdown } = await import('./static/js/ui.js');
  const html = markdown('# 标题\n\n- a\n- b\n\n| x | y |\n|---|---|\n| 1 | 2 |\n\n`code` **bold** [l](http://e)\n\n> quote\n');
  for (const needle of ['<h1>标题</h1>', '<li>a</li>', '<table>', '<td style="text-align:left">1</td>', '<code>code</code>', '<strong>bold</strong>', '<blockquote>quote</blockquote>']) {
    if (!html.includes(needle)) throw new Error(`missing ${needle} in markdown output`);
  }
});

await unit('ui.toCSV escaping', async () => {
  const { toCSV } = await import('./static/js/ui.js');
  const csv = toCSV([{ key: 'a', label: 'A' }, { key: 'b', label: 'B' }], [{ a: 'x,y', b: 'q"z' }]);
  if (!csv.includes('"x,y"') || !csv.includes('"q""z"')) throw new Error(`bad csv: ${csv}`);
});

await unit('markdown: 真实报告（8 张表、表格列数一致）', async () => {
  const { markdown } = await import('./static/js/ui.js');
  const report = JSON.parse(readFileSync(join(TESTDATA, 'evidence_report.json'), 'utf8')).markdown;
  const html = markdown(report);
  const tables = html.match(/<table>/g) || [];
  if (tables.length !== 8) throw new Error(`expected 8 tables, got ${tables.length}`);
  if (/<p>\s*\|/.test(html)) throw new Error('a table row leaked into a paragraph');
  if (!/<blockquote>/.test(html)) throw new Error('blockquotes not rendered');
  if (!/<h2>/.test(html)) throw new Error('headings not rendered');
  for (const chunk of html.split('<table>').slice(1)) {
    const body = chunk.slice(0, chunk.indexOf('</table>'));
    const head = (body.match(/<th[\s>]/g) || []).length;
    if (!head) continue;                       // title-less fallback table
    for (const row of body.split('<tr>').slice(2)) {
      const cells = (row.match(/<td[\s>]/g) || []).length;
      if (cells && cells !== head) throw new Error(`ragged table: ${head} headers vs ${cells} cells`);
    }
  }
});

await unit('markdown: 方法学文档（代码块、列表、行内代码）', async () => {
  const { markdown } = await import('./static/js/ui.js');
  const doc = JSON.parse(readFileSync(join(TESTDATA, 'docs_01.json'), 'utf8')).markdown;
  const html = markdown(doc);
  if (!/<pre><code/.test(html)) throw new Error('fenced code not rendered');
  if (!/<ul>|<ol>/.test(html)) throw new Error('lists not rendered');
  if (!/<code>/.test(html)) throw new Error('inline code not rendered');
  if (/^\s*```/m.test(html.replace(/<pre>[\s\S]*?<\/pre>/g, ''))) throw new Error('unclosed fence');
});

await unit('viewer projection / rotation / fit', async () => {
  const { StructureViewer } = await import('./static/js/viewer.js');
  const scene = JSON.parse(readFileSync(join(TESTDATA, 'scene_cj3.json'), 'utf8'));
  const host = new FakeNode('div');
  const viewer = new StructureViewer(host, {});
  viewer.setScene(scene);
  const before = viewer.view.rot.slice();
  viewer.view.rot = viewer.view.rot.map((v, i) => (i === 0 ? v + 0.1 : v));
  const p = viewer._project([0, 0, 0]);
  if (!p.every((v) => Number.isFinite(v))) throw new Error('projection produced non-finite values');
  if (!before.some((v) => Number.isFinite(v))) throw new Error('initial rotation invalid');
  viewer.setPocket({ xyz: [0, 0, 0, 1, 1, 1, 2, 2, 2], chem: [0, 1, 2], sub: [[0, 1, 2], ['A', 'B']] });
  viewer.setOptions({ showPocket: true, pocketMode: 'sub' });
  viewer.render();
  viewer.setOptions({ showPocket: false, showPeptides: false, spin: false });
  viewer.render();
  viewer.destroy();
});

/* ------------------------------------------------------------- shell checks */
const settle = async (ms = 30) => {
  for (let i = 0; i < 4; i += 1) {
    await new Promise((resolve) => setTimeout(resolve, ms));
    flushRaf();
  }
};

const fireHashChange = () => {
  for (const fn of winListeners.hashchange || []) fn({});
};

await unit('shell: boot, nav, routing, theme', async () => {
  const nav = documentShim.getElementById('nav');
  const view = documentShim.getElementById('view');
  const crumb = documentShim.getElementById('crumb');
  await import('./static/js/app.js');
  await settle();
  const items = nav.querySelectorAll('.nav-item');
  if (items.length !== 9) throw new Error(`expected 9 nav items, got ${items.length}`);
  if (!items.some((i) => i.classList.contains('active'))) throw new Error('no active nav item');
  if (view.children.length === 0) throw new Error('view is empty after boot');
  if (!crumb.textContent.includes('概览')) throw new Error(`crumb not set: ${crumb.textContent}`);

  globalThis.location.hash = '#/pocket';
  fireHashChange();
  await settle();
  if (!crumb.textContent.includes('口袋')) throw new Error(`routing failed: ${crumb.textContent}`);
  if (!nav.querySelectorAll('.nav-item').some((i) => i.dataset.id === 'pocket' && i.classList.contains('active'))) {
    throw new Error('nav highlight did not move to 口袋');
  }

  globalThis.location.hash = '#/structure?ligand=3&pocket=1';
  fireHashChange();
  await settle();
  if (!view.querySelectorAll('canvas').length) throw new Error('structure route did not mount the viewer');
  if (!crumb.textContent.includes('3D 结构')) throw new Error(`params route failed: ${crumb.textContent}`);

  const themeBtn = documentShim.getElementById('btn-theme');
  const before = documentShim.documentElement.dataset.theme;
  themeBtn.click();
  const after = documentShim.documentElement.dataset.theme;
  if (before === after) throw new Error('theme toggle had no effect');
  themeBtn.click();
});

/* ------------------------------------------------------------------ report */
let failures = 0;
console.log('');
for (const [name, status, detail] of results) {
  if (status !== 'ok') failures += 1;
  console.log(`${status === 'ok' ? 'PASS' : 'FAIL'}  ${name.padEnd(46)} ${detail.split('\n')[0]}`);
  if (status !== 'ok' && detail) console.log(detail.split('\n').map((l) => `        ${l}`).join('\n'));
}
console.log(`\n${results.length - failures}/${results.length} checks passed`);
process.exit(failures ? 1 : 0);
