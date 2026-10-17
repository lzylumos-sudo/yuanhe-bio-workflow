/* ui.js -- tiny DOM toolkit: hyperscript, formatting, tables, charts, markdown.
   No framework, no build step: everything here is plain ES modules. */

/* ---------------- DOM ---------------- */
export function el(tag, props, ...children) {
  const node = document.createElement(tag);
  applyProps(node, props);
  append(node, children);
  return node;
}

export function svg(tag, attrs, ...children) {
  const node = document.createElementNS('http://www.w3.org/2000/svg', tag);
  for (const [key, value] of Object.entries(attrs || {})) {
    if (value === undefined || value === null) continue;
    if (key === 'class') node.setAttribute('class', value);
    else if (key === 'text') node.textContent = value;
    else node.setAttribute(key, value);
  }
  for (const child of children.flat(3)) {
    if (child === null || child === undefined || child === false) continue;
    node.appendChild(typeof child === 'string' ? document.createTextNode(child) : child);
  }
  return node;
}

function applyProps(node, props) {
  if (!props) return;
  for (const [key, value] of Object.entries(props)) {
    if (value === undefined || value === null || value === false) continue;
    if (key === 'class' || key === 'className') node.className = Array.isArray(value) ? value.filter(Boolean).join(' ') : value;
    else if (key === 'text') node.textContent = value;
    else if (key === 'html') node.innerHTML = value;
    else if (key === 'style' && typeof value === 'object') Object.assign(node.style, value);
    else if (key === 'dataset') Object.assign(node.dataset, value);
    else if (key.startsWith('on') && typeof value === 'function') node.addEventListener(key.slice(2).toLowerCase(), value);
    else if (key === 'value') node.value = value;
    else if (key === 'checked' || key === 'disabled' || key === 'selected') node[key] = !!value;
    else if (key === 'ref' && typeof value === 'function') value(node);
    else node.setAttribute(key, value === true ? '' : value);
  }
}

function append(node, children) {
  for (const child of children.flat(6)) {
    if (child === null || child === undefined || child === false) continue;
    node.appendChild(child instanceof Node ? child : document.createTextNode(String(child)));
  }
}

export const h = el;
export const div = (props, ...kids) => el('div', props, ...kids);
export const span = (props, ...kids) => el('span', props, ...kids);

/* ---------------- escaping & formatting ---------------- */
export function escapeHtml(text) {
  return String(text ?? '')
    .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;').replace(/'/g, '&#39;');
}

export function fmt(value, digits = 2) {
  if (value === null || value === undefined || value === '') return '—';
  if (typeof value !== 'number' || !Number.isFinite(value)) return String(value);
  if (Number.isInteger(value) && digits === 0) return value.toLocaleString('en-US');
  if (Math.abs(value) >= 10000) return value.toLocaleString('en-US', { maximumFractionDigits: 0 });
  return value.toFixed(digits);
}

export function fmtInt(value) {
  if (value === null || value === undefined || value === '') return '—';
  const n = Number(value);
  return Number.isFinite(n) ? n.toLocaleString('en-US') : String(value);
}

export function fmtPct(value, digits = 0) {
  const n = Number(value);
  if (!Number.isFinite(n)) return '—';
  return `${(n * 100).toFixed(digits)}%`;
}

export function fmtBytes(bytes) {
  const n = Number(bytes);
  if (!Number.isFinite(n)) return '—';
  const units = ['B', 'KB', 'MB', 'GB'];
  let i = 0;
  let v = n;
  while (v >= 1024 && i < units.length - 1) { v /= 1024; i += 1; }
  return `${v.toFixed(i ? 1 : 0)} ${units[i]}`;
}

export function timeAgo(seconds) {
  if (!seconds) return '—';
  const diff = Date.now() / 1000 - Number(seconds);
  if (!Number.isFinite(diff)) return '—';
  if (diff < 60) return 'Just now';
  if (diff < 3600) return `${Math.floor(diff / 60)} minutes ago`;
  if (diff < 86400) return `${Math.floor(diff / 3600)} hours ago`;
  return `${Math.floor(diff / 86400)} days ago`;
}

export function fmtDate(seconds) {
  if (!seconds) return '—';
  const d = new Date(Number(seconds) * 1000);
  if (Number.isNaN(d.getTime())) return '—';
  const pad = (n) => String(n).padStart(2, '0');
  return `${d.getFullYear()}-${pad(d.getMonth() + 1)}-${pad(d.getDate())} ${pad(d.getHours())}:${pad(d.getMinutes())}`;
}

/* ---------------- components ---------------- */
export function card(opts, ...children) {
  const { title, right, note, class: cls = '' } = opts || {};
  const head = (title || right)
    ? div({ class: 'card-head' },
        title ? el('h3', { text: title }) : null,
        right ? div({ class: 'right' }, right) : null)
    : null;
  return div({ class: `card ${cls}` },
    head,
    ...children,
    note ? div({ class: 'card-note', text: note }) : null);
}

export function stat({ label, value, unit, foot, tone = '', small = false }) {
  return div({ class: `stat ${tone}` },
    div({ class: 'label', text: label }),
    div({ class: `value ${small ? 'sm' : ''}` },
      String(value),
      unit ? el('span', { class: 'unit', text: unit }) : null),
    foot ? div({ class: 'foot' }, foot) : null);
}

export function chip(text, cls = '') {
  return span({ class: `chip ${cls}`, text });
}

export function swatch(color) {
  return span({ class: 'swatch', style: { background: color } });
}

const STATUS_CLASS = {
  PASS: 'pass', FLAG: 'flag', FAIL: 'fail', WEAK: 'weak',
  MISMATCH: 'mismatch', NO_RESIDUE: 'fail', OK: 'pass',
};
export function statusChip(status) {
  const key = String(status || '').toUpperCase();
  return chip(key || '—', STATUS_CLASS[key] || 'info');
}

export function residChip(label, color) {
  return span({ class: 'resid' },
    color ? span({ class: 'swatch', style: { background: color } }) : null,
    color ? ' ' : null,
    label);
}

export function kvList(entries) {
  const rows = [];
  for (const [key, value] of entries) {
    if (value === null || value === undefined || value === '') continue;
    rows.push(el('div', { class: 'k', text: key }), el('div', { class: 'v' }, value));
  }
  return div({ class: 'kv' }, ...rows);
}

export function barList(items, { max, format = (v) => fmt(v, 2) } = {}) {
  const peak = max || Math.max(...items.map((it) => Number(it.value) || 0), 1e-9);
  return div({ class: 'bar-list' },
    ...items.map((it) => div({ class: 'row' },
      div({ class: 'k', text: it.label, title: it.title || it.label }),
      div({ class: 'track' },
        div({
          class: `fill ${it.tone || ''}`,
          style: { width: `${Math.max(1.5, (Number(it.value) / peak) * 100)}%`, background: it.color || undefined },
        })),
      div({ class: 'num', text: format(it.value) }))));
}

export function funnel(rows) {
  const peak = Math.max(...rows.map((r) => Number(r.value) || 0), 1);
  return div({ class: 'funnel' },
    ...rows.map((r) => div({ class: `funnel-row ${r.final ? 'final' : ''}` },
      div({ class: 'k', text: r.label, title: r.hint || '' }),
      div({ class: 'track' },
        div({ class: 'fill', style: { width: `${Math.max(1.2, (Number(r.value) / peak) * 100)}%` } })),
      div({ class: 'n', text: fmtInt(r.value) }))));
}

export function meterRow(label, value, max, color) {
  const pct = Math.max(0, Math.min(1, Number(value) / (max || 1)));
  return div({ class: 'meter-row' },
    div({ class: 'k', text: label }),
    div({ class: 'meter' }, el('i', { style: { width: `${pct * 100}%`, background: color } })),
    div({ class: 'v', text: fmt(value, 2) }));
}

export function legend(items) {
  return div({ class: 'legend' },
    ...items.map((it) => div({ class: 'item' },
      span({ class: 'swatch', style: { background: it.color } }), it.label)));
}

export function emptyState(text) {
  return div({ class: 'empty-state', text });
}

export function loadingPane(text = 'Loading...') {
  return div({ class: 'loading-pane' }, span({ class: 'spinner' }), text);
}

/* ---------------- toast ---------------- */
export function toast(message, type = 'info', ms = 3600) {
  const stack = document.getElementById('toasts');
  if (!stack) return;
  const node = el('div', { class: `toast ${type}`, text: message });
  stack.appendChild(node);
  setTimeout(() => {
    node.style.transition = 'opacity .2s, transform .2s';
    node.style.opacity = '0';
    node.style.transform = 'translateY(6px)';
    setTimeout(() => node.remove(), 220);
  }, ms);
}

/* ---------------- data table ---------------- */
export class DataTable {
  constructor({
    columns, rows = [], pageSize = 50, empty = 'No matching records',
    onRowClick = null, initialSort = null, keyField = null, pageSizes = [25, 50, 100, 250],
  }) {
    this.columns = columns;
    this.allRows = rows;
    this.viewRows = rows;
    this.pageSize = pageSize;
    this.page = 0;
    this.emptyText = empty;
    this.onRowClick = onRowClick;
    this.keyField = keyField;
    this.sort = initialSort ? { ...initialSort } : null;
    this.selectedKey = null;
    this.pageSizes = pageSizes;
    this.el = this._build();
  }

  _build() {
    this.tbody = el('tbody');
    this.headRow = el('tr');
    for (const col of this.columns) {
      const th = el('th', {
        class: [col.num ? 'num' : '', col.sortable === false ? 'no-sort' : ''].filter(Boolean).join(' '),
        style: col.width ? { width: col.width } : undefined,
        title: col.title || col.label,
      }, col.label,
        el('span', { class: 'arrow', text: '▲' }));
      if (col.sortable !== false) {
        th.addEventListener('click', () => this.sortBy(col.key));
      }
      col._th = th;
      this.headRow.appendChild(th);
    }
    this.table = el('table', null, el('thead', null, this.headRow), this.tbody);
    this.wrap = div({ class: 'table-wrap' }, this.table);
    this.foot = div({ class: 'table-foot' });
    this.body = div({ class: 'table-body' }, this.wrap, this.foot);
    this._renderHead();
    this._render();
    return this.body;
  }

  _renderHead() {
    for (const col of this.columns) {
      const active = this.sort && this.sort.key === col.key;
      col._th?.classList.toggle('sorted', !!active);
      const arrow = col._th?.querySelector('.arrow');
      if (arrow) arrow.textContent = active ? (this.sort.dir === 'asc' ? '▲' : '▼') : '▲';
    }
  }

  _sortRows(rows) {
    if (!this.sort) return rows;
    const { key, dir } = this.sort;
    const col = this.columns.find((c) => c.key === key);
    const sign = dir === 'asc' ? 1 : -1;
    return [...rows].sort((a, b) => {
      const va = col?.sortValue ? col.sortValue(a) : a[key];
      const vb = col?.sortValue ? col.sortValue(b) : b[key];
      const na = Number(va);
      const nb = Number(vb);
      if (Number.isFinite(na) && Number.isFinite(nb)) return (na - nb) * sign;
      const sa = String(va ?? '');
      const sb = String(vb ?? '');
      return sa.localeCompare(sb, 'zh-Hans-CN') * sign;
    });
  }

  _render() {
    const sorted = this._sortRows(this.viewRows);
    const total = sorted.length;
    const pages = Math.max(1, Math.ceil(total / this.pageSize));
    this.page = Math.min(this.page, pages - 1);
    const start = this.page * this.pageSize;
    const slice = sorted.slice(start, start + this.pageSize);
    this.tbody.replaceChildren();
    if (!slice.length) {
      const td = el('td', { colspan: this.columns.length }, emptyState(this.emptyText));
      this.tbody.appendChild(el('tr', null, td));
    } else {
      for (const row of slice) {
        const tr = el('tr', { class: this.onRowClick ? 'clickable' : '' });
        if (this.keyField) tr.dataset.key = String(row[this.keyField]);
        if (this.selectedKey !== null && this.keyField && String(row[this.keyField]) === String(this.selectedKey)) {
          tr.classList.add('selected');
        }
        for (const col of this.columns) {
          const value = col.render ? col.render(row) : row[col.key];
          const td = el('td', {
            class: [col.num ? 'num' : '', col.mono ? 'mono' : '', col.wrap ? 'wrap' : ''].filter(Boolean).join(' '),
          });
          if (value instanceof Node) td.appendChild(value);
          else if (col.html && typeof value === 'string') td.innerHTML = value;
          else td.textContent = value === null || value === undefined || value === '' ? '—' : String(value);
          tr.appendChild(td);
        }
        if (this.onRowClick) tr.addEventListener('click', () => {
          this.select(this.keyField ? row[this.keyField] : null);
          this.onRowClick(row);
        });
        this.tbody.appendChild(tr);
      }
    }
    this._renderFoot(total, pages, start, slice.length);
  }

  _renderFoot(total, pages, start, shown) {
    const range = total ? `${start + 1}–${start + shown} / ${total}` : '0 / 0';
    const pager = el('div', { class: 'split' },
      el('button', { class: 'btn sm ghost', text: '‹', disabled: this.page === 0, onclick: () => { this.page = Math.max(0, this.page - 1); this._render(); } }),
      span({ class: 'mono-sm', text: `${this.page + 1} / ${pages}` }),
      el('button', { class: 'btn sm ghost', text: '›', disabled: this.page >= pages - 1, onclick: () => { this.page = Math.min(pages - 1, this.page + 1); this._render(); } }));
    const sizeSel = el('select', { class: 'sm', onchange: (e) => { this.pageSize = Number(e.target.value); this.page = 0; this._render(); } },
      ...this.pageSizes.map((n) => el('option', { value: n, selected: n === this.pageSize }, `${n} entries/page`)));
    this.foot.replaceChildren(
      span({ class: 'grow', text: range }),
      pager,
      sizeSel,
    );
  }

  sortBy(key, dir) {
    const nextDir = dir || (this.sort && this.sort.key === key ? (this.sort.dir === 'asc' ? 'desc' : 'asc') : 'desc');
    this.sort = { key, dir: nextDir };
    this.page = 0;
    this._renderHead();
    this._render();
  }

  setRows(rows) {
    this.viewRows = rows || [];
    this.page = 0;
    this._render();
  }

  select(key) {
    this.selectedKey = key;
    for (const tr of this.tbody.querySelectorAll('tr')) {
      tr.classList.toggle('selected', this.keyField && tr.dataset.key === String(key));
    }
  }

  get rows() { return this.viewRows; }
}

/* ---------------- markdown ---------------- */
function inline(text) {
  let out = escapeHtml(text);
  const codes = [];
  out = out.replace(/`([^`]+)`/g, (_, code) => {
    codes.push(code);
    return `\u0000${codes.length - 1}\u0000`;
  });
  out = out.replace(/\*\*([^*]+)\*\*/g, '<strong>$1</strong>');
  out = out.replace(/(^|[^*\w])\*([^*\s][^*]*)\*/g, '$1<em>$2</em>');
  out = out.replace(/~~([^~]+)~~/g, '<del>$1</del>');
  out = out.replace(/\[([^\]]+)\]\(([^)\s]+)\)/g, '<a href="$2" target="_blank" rel="noreferrer">$1</a>');
  out = out.replace(/\u0000(\d+)\u0000/g, (_, i) => `<code>${escapeHtml(codes[Number(i)])}</code>`);
  return out;
}

export function markdown(text) {
  const lines = String(text || '').replace(/\r\n?/g, '\n').split('\n');
  const out = [];
  let i = 0;
  let listStack = null;

  const closeList = () => {
    if (listStack) { out.push(`</${listStack}>`); listStack = null; }
  };

  while (i < lines.length) {
    const line = lines[i];
    if (/^\s*$/.test(line)) { closeList(); i += 1; continue; }

    if (/^\s*```/.test(line)) {
      closeList();
      const lang = line.replace(/^\s*```/, '').trim();
      const buf = [];
      i += 1;
      while (i < lines.length && !/^\s*```/.test(lines[i])) { buf.push(lines[i]); i += 1; }
      i += 1;
      out.push(`<pre><code${lang ? ` class="lang-${escapeHtml(lang)}"` : ''}>${escapeHtml(buf.join('\n'))}</code></pre>`);
      continue;
    }

    const heading = /^(#{1,4})\s+(.*)$/.exec(line);
    if (heading) {
      closeList();
      const level = heading[1].length;
      out.push(`<h${level}>${inline(heading[2].trim())}</h${level}>`);
      i += 1;
      continue;
    }

    if (/^\s*([-*_])\1{2,}\s*$/.test(line)) { closeList(); out.push('<hr>'); i += 1; continue; }

    if (/^\s*\|/.test(line)) {
      closeList();
      const rows = [];
      while (i < lines.length && /^\s*\|/.test(lines[i])) {
        const raw = lines[i].trim().replace(/^\|/, '').replace(/\|$/, '');
        rows.push(raw.split('|').map((cell) => cell.trim()));
        i += 1;
      }
      if (rows.length >= 2 && rows[1].every((c) => /^:?-{2,}:?$/.test(c) || c === '')) {
        const aligns = rows[1].map((c) => (c.startsWith(':') && c.endsWith(':') ? 'center' : c.endsWith(':') ? 'right' : 'left'));
        const head = rows[0];
        out.push('<table><thead><tr>' + head.map((c, idx) => `<th style="text-align:${aligns[idx] || 'left'}">${inline(c)}</th>`).join('') + '</tr></thead><tbody>');
        for (const row of rows.slice(2)) {
          out.push('<tr>' + row.map((c, idx) => `<td style="text-align:${aligns[idx] || 'left'}">${inline(c)}</td>`).join('') + '</tr>');
        }
        out.push('</tbody></table>');
      } else {
        out.push('<table><tbody>' + rows.map((row) => '<tr>' + row.map((c) => `<td>${inline(c)}</td>`).join('') + '</tr>').join('') + '</tbody></table>');
      }
      continue;
    }

    if (/^\s*>\s?/.test(line)) {
      closeList();
      const buf = [];
      while (i < lines.length && /^\s*>\s?/.test(lines[i])) {
        buf.push(lines[i].replace(/^\s*>\s?/, ''));
        i += 1;
      }
      out.push(`<blockquote>${buf.map((l) => inline(l)).join('<br>')}</blockquote>`);
      continue;
    }

    const ol = /^\s*(\d+)[.)]\s+(.*)$/.exec(line);
    const ul = /^\s*[-*+]\s+(.*)$/.exec(line);
    if (ol || ul) {
      const kind = ol ? 'ol' : 'ul';
      if (listStack !== kind) {
        closeList();
        out.push(`<${kind}>`);
        listStack = kind;
      }
      out.push(`<li>${inline(ol ? ol[2] : ul[1])}</li>`);
      i += 1;
      continue;
    }

    closeList();
    const para = [line];
    i += 1;
    while (i < lines.length) {
      const next = lines[i];
      if (/^\s*$/.test(next) || /^\s*(#{1,4}\s|\||>|```)/.test(next) || /^\s*[-*+]\s/.test(next) || /^\s*\d+[.)]\s/.test(next)) break;
      para.push(next);
      i += 1;
    }
    out.push(`<p>${para.map((l) => inline(l)).join('<br>')}</p>`);
  }
  closeList();
  return out.join('\n');
}

export function markdownNode(text) {
  return div({ class: 'md', html: markdown(text) });
}

/* ---------------- drawer ---------------- */
export function openDrawer({ title, sub, body, footer, width } = {}) {
  const backdrop = div({ class: 'drawer-backdrop' });
  const drawer = div({ class: 'drawer', style: width ? { width } : undefined },
    div({ class: 'drawer-head' },
      div({ style: { flex: '1', minWidth: '0' } },
        el('h2', { text: title || '' }),
        sub ? div({ class: 'sub' }, sub) : null),
      el('button', { class: 'icon-btn', title: 'Close', text: '✕', onclick: () => close() })),
    div({ class: 'drawer-body' }, body),
    footer ? div({ class: 'drawer-foot' }, footer) : null);
  document.body.append(backdrop, drawer);
  requestAnimationFrame(() => { backdrop.classList.add('show'); drawer.classList.add('show'); });
  const onKey = (event) => { if (event.key === 'Escape') close(); };
  document.addEventListener('keydown', onKey);
  backdrop.addEventListener('click', () => close());
  function close() {
    document.removeEventListener('keydown', onKey);
    backdrop.classList.remove('show');
    drawer.classList.remove('show');
    setTimeout(() => { backdrop.remove(); drawer.remove(); }, 220);
  }
  return { close, drawer };
}

/* ---------------- download helpers ---------------- */
export function downloadText(filename, text, mime = 'text/plain;charset=utf-8') {
  const blob = new Blob([text], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = el('a', { href: url, download: filename });
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1500);
}

export function toCSV(columns, rows) {
  const esc = (v) => {
    const s = v === null || v === undefined ? '' : String(v);
    return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
  };
  const head = columns.map((c) => esc(c.label)).join(',');
  const body = rows.map((row) => columns.map((c) => esc(c.csv ? c.csv(row) : row[c.key])).join(',')).join('\n');
  return `${head}\n${body}\n`;
}

export function formatNow() {
  const d = new Date();
  const pad = (n) => String(n).padStart(2, '0');
  return `${d.getFullYear()}${pad(d.getMonth() + 1)}${pad(d.getDate())}-${pad(d.getHours())}${pad(d.getMinutes())}`;
}
