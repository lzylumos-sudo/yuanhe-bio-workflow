/* pages/pharmacophore.js -- consensus pharmacophore features, feature pairs and
   the design axes the scaffold assembler enumerates along. */

import { api } from '../api.js';
import {
  el, div, span, card, stat, chip, DataTable, svg, fmt, loadingPane,
} from '../ui.js';

const TYPE_COLOR = {
  anion: '#e55353', cation: '#4f8ef7', aromatic: '#c99a3a', donor: '#38a169', hydrophobic: '#8b97a8',
};
const TYPE_LABEL = {
  anion: '阴离子（羧酸/磷酸）', cation: '阳离子（碱性胺）', aromatic: '芳香/π',
  donor: '氢键给体', hydrophobic: '疏水',
};

export async function mount(root, ctx) {
  root.appendChild(loadingPane('正在读取药效团与设计轴…'));
  const data = await api.pharmacophore();
  const features = data.features || [];
  const pairs = data.pairs || [];
  const axes = data.axes || [];

  root.replaceChildren();
  root.appendChild(div({ class: 'page-head' },
    div(null,
      el('h1', { text: '共识药效团与设计轴' }),
      div({ class: 'sub', text: '每条肽的功能原子（丙氨酸扫描 + 自发耐药突变 + 三肽趋同接触所定义）被投到同一 MurJ 坐标系里聚类；特征之间的两两距离就是骨架必须满足的几何约束。' })),
    div({ class: 'actions' },
      el('a', { class: 'btn', href: '/files/results/pharmacophore.json?download=1', text: '下载 JSON' }),
      el('button', { class: 'btn', text: '在 3D 中查看特征', onclick: () => ctx.navigate('structure', { ligand: '0' }) }))));

  /* ---- explanation ---- */
  root.appendChild(card({ class: 'section', title: '为什么是“轴”而不是“形状匹配”' },
    div({ class: 'grid cols-3', style: { gap: '12px' } },
      div(null,
        div({ class: 'tiny muted', text: '1 · 药效团对' }),
        div({ class: 'tiny', text: '任取两个权重足够高的特征，得到一条有方向、有长度的设计轴（9–18 Å）。' })),
      div(null,
        div({ class: 'tiny muted', text: '2 · 一维长度约束' }),
        div({ class: 'tiny', text: '骨架是“头帽–桥–环–桥–…–尾帽”的刚性单元链，每个单元跨距固定；能否同时够到两端，变成一个一维区间问题。' })),
      div(null,
        div({ class: 'tiny muted', text: '3 · 确定性求解' }),
        div({ class: 'tiny', text: '区间查找精确求解，因此 20,307 个骨架可穷举、可复现、可归因——不是随机采样。' })))));

  /* ---- feature cards ---- */
  const cards = div({ class: 'grid cols-auto', style: { gap: '10px' } },
    ...features.map((f) => card({ class: 'tight' },
      div({ class: 'split' },
        span({ class: 'chip mono', style: { borderColor: TYPE_COLOR[f.type], color: TYPE_COLOR[f.type] }, text: f.fid }),
        span({ class: 'tiny muted', text: `${f.type}` }),
        span({ class: 'right chip mono', style: { marginLeft: 'auto' }, text: `w ${fmt(f.weight, 2)}` })),
      div({ class: 'split tiny muted', style: { marginTop: '6px' } },
        span({ text: `趋同 ${fmt(f.convergence, 2)}` }),
        span({ text: `容差 ±${fmt(f.tol, 1)} Å` }),
        span({ text: `spread ${fmt(f.spread, 2)}` })),
      div({ class: 'chips', style: { marginTop: '6px' } },
        ...(f.murj_contacts || []).slice(0, 6).map((c) => span({ class: 'resid', text: c })),
        (f.murj_contacts || []).length > 6 ? span({ class: 'tiny muted', text: `+${f.murj_contacts.length - 6}` }) : null),
      div({ class: 'tiny muted', style: { marginTop: '6px' }, text: (f.sources || []).join(' · ') }))));

  root.appendChild(card({ title: `共识特征（${features.length}）`, class: 'section' }, cards,
    div({ class: 'hr' }),
    div({ class: 'legend' },
      ...Object.entries(TYPE_LABEL).map(([key, label]) => div({ class: 'item' },
        span({ class: 'swatch', style: { background: TYPE_COLOR[key] } }), `${label}`)))));

  /* ---- distance vs weight scatter ---- */
  const scatter = buildScatter(pairs, axes);
  root.appendChild(card({ title: '特征对的距离—权重分布', class: 'section' },
    div({ class: 'tiny muted', style: { marginBottom: '6px' }, text: '灰点：全部特征对；彩色点：被 14 条设计轴采用的（阴离子–X，距离 9–18 Å，权重 ≥ 0.45）。' }),
    scatter));

  /* ---- pairs table ---- */
  const pairTable = new DataTable({
    columns: [
      { key: 'a', label: '特征 A', width: '70px' },
      { key: 'b', label: '特征 B', width: '70px' },
      { key: 'types', label: '类型', mono: true },
      { key: 'distance', label: '距离 Å', num: true, render: (r) => fmt(r.distance, 2) },
      { key: 'weight', label: '权重', num: true, render: (r) => fmt(r.weight, 2) },
    ],
    rows: pairs,
    pageSize: 12,
    initialSort: { key: 'weight', dir: 'desc' },
  });

  const axisTable = new DataTable({
    columns: [
      { key: 'axis', label: '设计轴', width: '90px' },
      { key: 'types', label: '类型', mono: true },
      { key: 'distance', label: '轴长 Å', num: true, render: (r) => fmt(r.distance, 2) },
      { key: 'weight', label: '权重', num: true, render: (r) => fmt(r.weight, 2) },
      {
        key: 'candidates', label: '枚举骨架', num: true,
        render: (r) => (r.candidates ? String(r.candidates) : span({ class: 'chip weak sm', text: '无候选' })),
      },
      { key: 'scored', label: '性质内打分', num: true },
      { key: 'seconds', label: '耗时 s', num: true, render: (r) => fmt(r.seconds, 1) },
    ],
    rows: axes,
    pageSize: 20,
    initialSort: { key: 'weight', dir: 'desc' },
  });
  const withCandidates = axes.filter((a) => a.candidates > 0).length;
  const deadAxes = axes.filter((a) => !a.candidates);

  root.appendChild(div({ class: 'section grid cols-2' },
    card({ title: `特征对（${pairs.length}）` }, pairTable.el),
    card({ title: `设计轴（${axes.length} 条，${withCandidates} 条有候选）` },
      axisTable.el,
      deadAxes.length
        ? div({ class: 'card-note', style: { marginTop: '8px' }, text: `${deadAxes.map((a) => a.axis).join('、')} 在现有刚性单元库下没有候选：它要求的跨距落在可枚举区间之外。这是负结果，保留在表里而不是删除——说明“哪些几何是被骨架库本身排除的”。` })
        : null)));

  return null;
}

function buildScatter(pairs, axes) {
  const W = 720;
  const H = 260;
  const pad = { l: 44, r: 16, t: 14, b: 30 };
  if (!pairs.length) return div({ class: 'empty-state', text: '无特征对数据' });
  const xs = pairs.map((p) => Number(p.distance));
  const ys = pairs.map((p) => Number(p.weight));
  const xmin = Math.min(...xs);
  const xmax = Math.max(...xs);
  const ymax = Math.max(...ys, 0.1);
  const x = (v) => pad.l + ((v - xmin) / (xmax - xmin || 1)) * (W - pad.l - pad.r);
  const y = (v) => H - pad.b - (v / (ymax * 1.12)) * (H - pad.t - pad.b);
  const axisSet = new Set(axes.map((a) => a.axis));

  const ticks = [];
  for (let t = Math.ceil(xmin / 2) * 2; t <= xmax; t += 2) {
    ticks.push(svg('line', { x1: x(t), y1: H - pad.b, x2: x(t), y2: H - pad.t, stroke: '#1d2937', 'stroke-width': 1 }));
    ticks.push(svg('text', { x: x(t), y: H - pad.b + 14, fill: '#67788c', 'font-size': 9, 'text-anchor': 'middle', text: String(t) }));
  }
  const gridY = [];
  for (let i = 1; i <= 4; i += 1) {
    const v = (ymax * 1.06) * (i / 4);
    gridY.push(svg('line', { x1: pad.l, y1: y(v), x2: W - pad.r, y2: y(v), stroke: '#16212e', 'stroke-width': 1 }));
    gridY.push(svg('text', { x: pad.l - 6, y: y(v) + 3, fill: '#67788c', 'font-size': 9, 'text-anchor': 'end', text: v.toFixed(2) }));
  }

  const dots = pairs.map((p) => {
    const isAxis = axisSet.has(`${p.a}-${p.b}`) || axisSet.has(`${p.b}-${p.a}`);
    return svg('circle', {
      cx: x(Number(p.distance)), cy: y(Number(p.weight)),
      r: isAxis ? 3.6 : 2,
      fill: isAxis ? '#38bdf8' : '#3c4c5f',
      opacity: isAxis ? 0.95 : 0.55,
    });
  });

  return svg('svg', { viewBox: `0 0 ${W} ${H}`, style: 'width:100%;height:auto;display:block' },
    svg('rect', { x: x(9), y: pad.t, width: Math.max(4, x(18) - x(9)), height: H - pad.t - pad.b, fill: 'rgba(56,189,248,.06)' }),
    svg('text', { x: (x(9) + x(18)) / 2, y: pad.t + 11, fill: '#38bdf8', 'font-size': 9, 'text-anchor': 'middle', text: '设计窗口 9–18 Å' }),
    ...ticks, ...gridY,
    svg('line', { x1: pad.l, y1: H - pad.b, x2: W - pad.r, y2: H - pad.b, stroke: '#2a3a4c' }),
    svg('line', { x1: pad.l, y1: pad.t, x2: pad.l, y2: H - pad.b, stroke: '#2a3a4c' }),
    ...dots,
    svg('text', { x: W / 2, y: H - 3, fill: '#67788c', 'font-size': 10, 'text-anchor': 'middle', text: '特征对距离 (Å)' }),
    svg('text', { x: 12, y: H / 2, fill: '#67788c', 'font-size': 10, transform: `rotate(-90 12 ${H / 2})`, 'text-anchor': 'middle', text: '权重' }));
}
