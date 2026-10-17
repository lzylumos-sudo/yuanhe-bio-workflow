/* flow/widgets.js -- small shared pieces: funnel, molecule table, 3D panel, 2D card. */
import { get } from '../api.js';
import { el, div, span } from '../ui.js';
import { StructureViewer } from '../viewer.js';

const nf = (v) => (v === null || v === undefined ? '—' : Number(v).toLocaleString('en-US'));

/* ---------------- funnel: bars with the drop between consecutive rows --------------- */
export function funnel(rows, { foot } = {}) {
  const max = Math.max(1, ...rows.map((r) => r.count || 0));
  const box = div({ class: 'funnel' });
  rows.forEach((r, i) => {
    const has = r.count !== null && r.count !== undefined;
    const prev = i > 0 ? rows[i - 1].count : null;
    const drop = has && prev && r.count < prev ? Math.round((1 - r.count / prev) * 100) : null;
    box.appendChild(div({ class: ['f-row', !has && 'empty'] },
      div({ class: 'f-top' },
        span({ class: 'lb', text: r.label }),
        span({ class: 'drop', text: drop ? `−${drop}%` : '' }),
        span({ class: 'num', text: nf(r.count) })),
      div({ class: 'f-bar' }, el('i', { style: { width: has ? `${Math.max(1.5, (r.count / max) * 100)}%` : '100%' } })),
      r.note ? div({ class: 'f-note', text: r.note }) : null));
  });
  if (foot) box.appendChild(div({ class: 'f-foot', text: foot }));
  return box;
}

/* ---------------- 3D panel --------------- */
const LEGEND = [
  ['#ecc94b', 'SglM'], ['#dd6b20', 'SglPP7'], ['#38a169', 'SglCJ3'],
];

export async function viewerPanel({ ligandHint = true, peptides = true, tip } = {}) {
  const scene0 = await get('/api/scene?ref=SglCJ3&radius=13');
  const scene = { ...scene0, ligand: null };
  const host = div({ class: 'vbox' });
  const over = div({ class: 'v-over' });
  host.appendChild(over);
  host.appendChild(div({ class: 'v-hint', text: '拖动旋转 · 滚轮缩放 · Shift 拖动平移 · 双击复位' }));
  const viewer = new StructureViewer(host, {
    onHover: (item, ev) => {
      if (!tip) return;
      if (!item) { tip.hidden = true; return; }
      tip.hidden = false;
      tip.textContent = item.label;
      tip.style.left = `${(ev?.clientX ?? 0) + 14}px`;
      tip.style.top = `${(ev?.clientY ?? 0) + 12}px`;
    },
  });
  viewer.setScene(scene);
  viewer.setOptions({ showPeptides: peptides, showFeatures: false, showHotspots: false, showLabels: true });

  const opts = {};
  const mk = (key, label, on) => {
    opts[key] = on;
    return el('label', { class: 'chk' },
      el('input', { type: 'checkbox', checked: on, onchange: (e) => viewer.setOptions({ [key]: e.target.checked }) }),
      label);
  };
  const bar = div({ class: 'vbar' },
    mk('showPocket', '口袋', true),
    mk('showPeptides', '三条肽', peptides),
    mk('showFeatures', '药效团', false),
    mk('showHotspots', '热点残基', false),
    mk('spin', '自转', false),
    el('button', { class: 'btn', text: '复位', onclick: () => viewer.fit() }),
    div({ class: 'legend' },
      ...LEGEND.map(([c, t]) => span(null, el('i', { style: { background: c } }), t)),
      ligandHint ? span(null, el('i', { style: { background: '#e0b64a' } }), '候选分子') : null));
  viewer.setOptions({ showPocket: true });

  // cavity voxels: computed on demand by the server, poll briefly
  let alive = true;
  (async () => {
    for (let i = 0; i < 40 && alive; i += 1) {
      try {
        const p = await get('/api/pocket/points?step=0.8', { force: i > 0 });
        if (p && p.status === 'ready') { viewer.setPocket(p); return; }
      } catch (e) { return; }
      await new Promise((r) => setTimeout(r, 1500));
    }
  })();

  return {
    root: div(null, host, bar),
    viewer,
    scene,
    setLigand(lig, caption) {
      scene.ligand = lig || null;
      viewer.scene = scene;
      over.replaceChildren(
        div({ text: `参考框架  ${scene.ref_pdbid} (SglCJ3)` }),
        caption ? div(null, span({ text: '配体  ' }), el('b', { text: caption })) : null);
      viewer.requestRender();
    },
    destroy() { alive = false; viewer.destroy(); if (tip) tip.hidden = true; },
  };
}

/* ---------------- molecule table --------------- */
export function moleculeTable(data, { onSelect, limit = 200 } = {}) {
  const cols = data.columns;
  let sortKey = null; let asc = false; let q = '';
  let rows = data.rows;
  const body = el('tbody');
  const head = el('tr');
  const more = div({ class: 'more' });
  let selected = null;

  const cell = (c, r) => {
    const v = r[c.key];
    if (typeof v === 'boolean') return el('td', { class: v ? 'yes' : 'no', text: v ? '✓' : '·' });
    const td = el('td');
    if (c.key === 'smiles') { td.className = 'smi'; td.title = v ?? ''; td.textContent = v ?? ''; return td; }
    if (c.key === 'id') td.className = 'id';
    td.textContent = v === null || v === undefined ? '—' : (typeof v === 'number' && !Number.isInteger(v) ? v.toFixed(2) : v);
    return td;
  };

  function draw() {
    let list = rows;
    if (q) list = list.filter((r) => `${r.id} ${r.smiles}`.toLowerCase().includes(q));
    if (sortKey) {
      list = [...list].sort((a, b) => {
        const x = a[sortKey]; const y = b[sortKey];
        if (x === y) return 0;
        if (x === null || x === undefined) return 1;
        if (y === null || y === undefined) return -1;
        return (x > y ? 1 : -1) * (asc ? 1 : -1);
      });
    }
    body.replaceChildren();
    for (const r of list.slice(0, limit)) {
      const tr = el('tr', { class: selected === r.id ? 'sel' : '', onclick: () => { select(r); } }, ...cols.map((c) => cell(c, r)));
      tr.dataset.id = r.id;
      body.appendChild(tr);
    }
    more.textContent = list.length > limit ? `仅显示前 ${limit} 行，共 ${nf(list.length)} 行，可搜索或排序` : `共 ${nf(list.length)} 行`;
    for (const th of head.children) {
      th.className = th.dataset.key === sortKey ? `sorted${asc ? ' asc' : ''}` : '';
    }
  }
  function select(r) {
    selected = r.id;
    for (const tr of body.children) tr.classList.toggle('sel', tr.dataset.id === String(r.id));
    onSelect?.(r);
  }
  for (const c of cols) {
    head.appendChild(el('th', {
      text: c.label,
      dataset: { key: c.key },
      onclick: () => { if (sortKey === c.key) asc = !asc; else { sortKey = c.key; asc = false; } draw(); },
    }));
  }
  const search = el('input', {
    type: 'search', placeholder: '搜索 ID / SMILES', 'aria-label': '搜索分子',
    oninput: (e) => { q = e.target.value.trim().toLowerCase(); draw(); },
  });
  const root = div(null,
    div({ class: 'tools' }, search),
    div({ class: 'tbl-wrap' }, el('table', { class: 'tbl' }, el('thead', null, head), body)),
    more);
  draw();
  return { root, select, first: () => rows[0], count: () => rows.length, rows: () => rows };
}

/* ---------------- 2D structure + properties --------------- */
export function molCard() {
  const draw = div({ class: 'draw', text: '—' });
  const info = div();
  const root = div({ class: 'mol' }, draw, info);
  return {
    root,
    async show(row) {
      if (!row) { draw.textContent = '—'; info.replaceChildren(div({ class: 'note', text: '从列表中选择一个分子' })); return; }
      const kv = el('dl', { class: 'kv' });
      for (const [k, v] of Object.entries(row)) {
        if (['id', 'smiles'].includes(k) || k.startsWith('_') || v === null || v === undefined || v === '') continue;
        kv.append(el('dt', { text: k }), el('dd', { text: typeof v === 'boolean' ? (v ? '是' : '否') : (typeof v === 'number' && !Number.isInteger(v) ? v.toFixed(2) : v) }));
      }
      info.replaceChildren(div({ class: 'id', text: row.id }), div({ class: 'smi', text: row.smiles || '' }), kv);
      draw.textContent = '…';
      try {
        const r = await get(`/api/mol2d?smiles=${encodeURIComponent(row.smiles || '')}`);
        if (r.svg) draw.innerHTML = r.svg; else draw.textContent = r.reason || '无 2D 图';
      } catch (e) { draw.textContent = '无 2D 图'; }
    },
  };
}
