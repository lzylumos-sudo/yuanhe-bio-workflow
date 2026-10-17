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
  anion: 'Anion (carboxylate//phosphate)', cation: 'Cation (basic amine)', aromatic: 'Aromatic/π',
  donor: 'Hydrogen-bond donor', hydrophobic: 'Hydrophobic',
};

export async function mount(root, ctx) {
  root.appendChild(loadingPane('Loading pharmacophores and design axes...'));
  const data = await api.pharmacophore();
  const features = data.features || [];
  const pairs = data.pairs || [];
  const axes = data.axes || [];

  root.replaceChildren();
  root.appendChild(div({ class: 'page-head' },
    div(null,
      el('h1', { text: 'Consensus pharmacophores and design axes' }),
      div({ class: 'sub', text: 'Peptide functional atoms defined by alanine scanning, + spontaneous resistance mutations, + and convergent contacts are clustered in a shared MurJ coordinate frame; pairwise feature distances define scaffold constraints.' })),
    div({ class: 'actions' },
      el('a', { class: 'btn', href: '/files/results/pharmacophore.json?download=1', text: 'Download JSON' }),
      el('button', { class: 'btn', text: 'in 3D view features', onclick: () => ctx.navigate('structure', { ligand: '0' }) }))));

  /* ---- explanation ---- */
  root.appendChild(card({ class: 'section', title: 'Why axes instead of shape matching?' },
    div({ class: 'grid cols-3', style: { gap: '12px' } },
      div(null,
        div({ class: 'tiny muted', text: '1 · Pharmacophore pair' }),
        div({ class: 'tiny', text: 'Two sufficiently weighted features define a directed design axis with a length (9–18 Å). ' })),
      div(null,
        div({ class: 'tiny muted', text: '2 · One-dimensional length constraint' }),
        div({ class: 'tiny', text: 'Scaffolds form rigid chains of head cap,–bridge,–Ring–bridge,–…–and tail cap units; fixed spans reduce endpoint reachability to an interval problem.' })),
      div(null,
        div({ class: 'tiny muted', text: '3 · Deterministic solution' }),
        div({ class: 'tiny', text: 'Interval lookup solves the constraint exactly, making 20,307 scaffolds enumerable, reproducible and attributable rather than randomly sampled.' })))));

  /* ---- feature cards ---- */
  const cards = div({ class: 'grid cols-auto', style: { gap: '10px' } },
    ...features.map((f) => card({ class: 'tight' },
      div({ class: 'split' },
        span({ class: 'chip mono', style: { borderColor: TYPE_COLOR[f.type], color: TYPE_COLOR[f.type] }, text: f.fid }),
        span({ class: 'tiny muted', text: `${f.type}` }),
        span({ class: 'right chip mono', style: { marginLeft: 'auto' }, text: `w ${fmt(f.weight, 2)}` })),
      div({ class: 'split tiny muted', style: { marginTop: '6px' } },
        span({ text: `Convergence ${fmt(f.convergence, 2)}` }),
        span({ text: `Tolerance ±${fmt(f.tol, 1)} Å` }),
        span({ text: `spread ${fmt(f.spread, 2)}` })),
      div({ class: 'chips', style: { marginTop: '6px' } },
        ...(f.murj_contacts || []).slice(0, 6).map((c) => span({ class: 'resid', text: c })),
        (f.murj_contacts || []).length > 6 ? span({ class: 'tiny muted', text: `+${f.murj_contacts.length - 6}` }) : null),
      div({ class: 'tiny muted', style: { marginTop: '6px' }, text: (f.sources || []).join(' · ') }))));

  root.appendChild(card({ title: `Consensus features (${features.length})`, class: 'section' }, cards,
    div({ class: 'hr' }),
    div({ class: 'legend' },
      ...Object.entries(TYPE_LABEL).map(([key, label]) => div({ class: 'item' },
        span({ class: 'swatch', style: { background: TYPE_COLOR[key] } }), `${label}`)))));

  /* ---- distance vs weight scatter ---- */
  const scatter = buildScatter(pairs, axes);
  root.appendChild(card({ title: 'Feature-pair distance versus weight', class: 'section' },
    div({ class: 'tiny muted', style: { marginBottom: '6px' }, text: 'Gray: all feature pairs; colored: pairs selected for 14 design axes (anionic–X, Distance 9–18 Å, Weight ≥ 0.45). ' }),
    scatter));

  /* ---- pairs table ---- */
  const pairTable = new DataTable({
    columns: [
      { key: 'a', label: 'Feature A', width: '70px' },
      { key: 'b', label: 'Feature B', width: '70px' },
      { key: 'types', label: 'Type', mono: true },
      { key: 'distance', label: 'Distance Å', num: true, render: (r) => fmt(r.distance, 2) },
      { key: 'weight', label: 'Weight', num: true, render: (r) => fmt(r.weight, 2) },
    ],
    rows: pairs,
    pageSize: 12,
    initialSort: { key: 'weight', dir: 'desc' },
  });

  const axisTable = new DataTable({
    columns: [
      { key: 'axis', label: 'Design axis', width: '90px' },
      { key: 'types', label: 'Type', mono: true },
      { key: 'distance', label: 'Axis length Å', num: true, render: (r) => fmt(r.distance, 2) },
      { key: 'weight', label: 'Weight', num: true, render: (r) => fmt(r.weight, 2) },
      {
        key: 'candidates', label: 'Enumerated scaffolds', num: true,
        render: (r) => (r.candidates ? String(r.candidates) : span({ class: 'chip weak sm', text: 'No candidates' })),
      },
      { key: 'scored', label: 'Scored within property envelope', num: true },
      { key: 'seconds', label: 'Elapsed time s', num: true, render: (r) => fmt(r.seconds, 1) },
    ],
    rows: axes,
    pageSize: 20,
    initialSort: { key: 'weight', dir: 'desc' },
  });
  const withCandidates = axes.filter((a) => a.candidates > 0).length;
  const deadAxes = axes.filter((a) => !a.candidates);

  root.appendChild(div({ class: 'section grid cols-2' },
    card({ title: `Feature pairs (${pairs.length})` }, pairTable.el),
    card({ title: `Design axes (${axes.length} entries,${withCandidates} with candidates)` },
      axisTable.el,
      deadAxes.length
        ? div({ class: 'card-note', style: { marginTop: '8px' }, text: `${deadAxes.map((a) => a.axis).join(', ')} No candidate in the current rigid-unit library: the requested span lies outside the enumerable intervals. This negative result is retained to show geometries excluded by the library.` })
        : null)));

  return null;
}

function buildScatter(pairs, axes) {
  const W = 720;
  const H = 260;
  const pad = { l: 44, r: 16, t: 14, b: 30 };
  if (!pairs.length) return div({ class: 'empty-state', text: 'No feature-pair data' });
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
    svg('text', { x: (x(9) + x(18)) / 2, y: pad.t + 11, fill: '#38bdf8', 'font-size': 9, 'text-anchor': 'middle', text: 'Design window 9–18 Å' }),
    ...ticks, ...gridY,
    svg('line', { x1: pad.l, y1: H - pad.b, x2: W - pad.r, y2: H - pad.b, stroke: '#2a3a4c' }),
    svg('line', { x1: pad.l, y1: pad.t, x2: pad.l, y2: H - pad.b, stroke: '#2a3a4c' }),
    ...dots,
    svg('text', { x: W / 2, y: H - 3, fill: '#67788c', 'font-size': 10, 'text-anchor': 'middle', text: 'Feature-pair distance (Å)' }),
    svg('text', { x: 12, y: H / 2, fill: '#67788c', 'font-size': 10, transform: `rotate(-90 12 ${H / 2})`, 'text-anchor': 'middle', text: 'Weight' }));
}
