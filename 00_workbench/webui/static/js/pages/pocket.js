/* pages/pocket.js -- the buried cavity: volume, chemistry, sub-pocket anatomy. */

import { api } from '../api.js';
import {
  el, div, span, card, stat, chip, DataTable, kvList, meterRow,
  fmt, fmtInt, loadingPane, legend,
} from '../ui.js';

const CHEM = [
  { key: 'pos', label: 'Basic (K/R/H)', color: '#4f8ef7' },
  { key: 'neg', label: 'Acidic (D/E)', color: '#e55353' },
  { key: 'polar', label: 'Polar (S/T/N/Q)', color: '#38a169' },
  { key: 'aromatic', label: 'Aromatic (Y/W/F)', color: '#c99a3a' },
  { key: 'hydrophobic', label: 'Hydrophobic', color: '#8b97a8' },
];

export async function mount(root, ctx) {
  root.appendChild(loadingPane('Loading pocket voxels...'));
  const data = await api.pocket();
  const pocket = data.pocket || {};
  const subs = pocket.subpockets || [];
  const total = Number(pocket.total_volume_A3) || 0;
  const deepest = subs.reduce((acc, s) => Math.max(acc, Number(s.mean_buriedness) || 0), 0);

  root.replaceChildren();
  root.appendChild(div({ class: 'page-head' },
    div(null,
      el('h1', { text: 'Buried cavity and subpockets' }),
      div({ class: 'sub', text: 'The peptide defines a negative image of the cavity.0.8 Å After voxelization, 50 rays estimate burial per voxel (SiteMap approximation); nearest subpocket centroids define chemical regions.' })),
    div({ class: 'actions' },
      el('a', { class: 'btn', href: '/files/results/pocket.json?download=1', text: 'Download JSON' }),
      el('button', {
        class: 'btn primary', text: 'in 3D view cavity',
        onclick: () => ctx.navigate('structure', { ligand: '0', pocket: '1' }),
      }))));

  root.appendChild(div({ class: 'grid cols-4', style: { marginBottom: '12px' } },
    card({ class: 'tight' }, stat({ label: 'Cavity volume', value: fmtInt(total), unit: 'Å³', small: true, tone: 'accent' })),
    card({ class: 'tight' }, stat({ label: 'Subpocket', value: subs.length, unit: 'items', small: true })),
    card({ class: 'tight' }, stat({ label: 'Maximum burial', value: fmt(deepest, 2), small: true, tone: 'good', foot: 'S2 Carboxylate gate' })),
    card({ class: 'tight' }, stat({
      label: 'Mean subpocket volume', value: subs.length ? fmt(subs.reduce((a, s) => a + Number(s.volume_A3), 0) / subs.length, 0) : '—', unit: 'Å³', small: true,
    }))));

  /* ---- sub-pocket anatomy ---- */
  const cards = div({ class: 'grid cols-2', style: { gap: '12px' } },
    ...subs.map((sp) => {
      const chem = sp.chemistry || {};
      const chemTotal = CHEM.reduce((acc, c) => acc + (chem[c.key] || 0), 0) || 1;
      return card({ class: 'tight' },
        div({ class: 'split' },
          span({ class: 'chip mono', text: sp.subpocket }),
          span({ class: 'right tiny muted', style: { marginLeft: 'auto' }, text: `${fmtInt(sp.volume_A3)} Å³ · ${fmtInt(sp.n_voxels)} Voxels` })),
        div({ style: { marginTop: '8px' } },
          div({ class: 'tiny muted', style: { marginBottom: '4px' }, text: 'Chemical composition by nearest residue' }),
          div({ style: { display: 'flex', height: '12px', borderRadius: '4px', overflow: 'hidden', border: '1px solid var(--line)' } },
            ...CHEM.map((c) => div({
              style: {
                width: `${((chem[c.key] || 0) / chemTotal) * 100}%`,
                background: c.color,
              },
              title: `${c.label}: ${chem[c.key] || 0}`,
            }))),
          div({ class: 'legend', style: { marginTop: '6px' } },
            ...CHEM.map((c) => div({ class: 'item' },
              span({ class: 'swatch', style: { background: c.color } }),
              span({ class: 'tiny', text: `${c.label} ${Math.round(((chem[c.key] || 0) / chemTotal) * 100)}%` }))))),
        div({ style: { marginTop: '10px', display: 'grid', gap: '6px' } },
          meterRow('Mean burial', sp.mean_buriedness, 1, 'var(--accent-2)'),
          div({ class: 'kv tiny' },
            el('div', { class: 'k', text: 'Centroid' }),
            el('div', { class: 'v mono-sm', text: `(${(sp.centroid || []).map((v) => fmt(v, 1)).join(', ')})` }),
            el('div', { class: 'k', text: 'Range Å' }),
            el('div', { class: 'v mono-sm', text: (sp.extent_A || []).map((v) => fmt(v, 1)).join(' × ') }))),
        div({ class: 'chips', style: { marginTop: '10px' } },
          ...(sp.residues || []).map((r) => span({ class: 'resid', text: r }))));
    }));

  root.appendChild(card({ title: 'Four subpockets', class: 'section' }, cards));

  /* ---- table + method ---- */
  const table = new DataTable({
    columns: [
      { key: 'subpocket', label: 'Subpocket', mono: true },
      { key: 'volume_A3', label: 'Volume Å³', num: true, render: (r) => fmt(r.volume_A3, 1) },
      { key: 'n_voxels', label: 'Voxels', num: true },
      { key: 'mean_buriedness', label: 'Mean burial', num: true, render: (r) => fmt(r.mean_buriedness, 3) },
      { key: 'chemistry', label: 'Dominant chemistry', render: (r) => dominant(r.chemistry) },
      { key: 'residues', label: 'Lining residues', wrap: true, render: (r) => (r.residues || []).join(' ') },
    ],
    rows: subs,
    pageSize: 10,
    initialSort: { key: 'volume_A3', dir: 'desc' },
  });

  const method = card({ title: 'How this value is calculated' },
    div({ class: 'tiny', style: { display: 'grid', gap: '8px' } },
      div(null, span({ class: 'muted', text: '1 · Voxelization:' }), 'Within the peptide envelope (±7 Å)using a 0.8 Å grid, remove voxels whose distance to MurJ heavy atoms is < 3.0 Å .'),
      div(null, span({ class: 'muted', text: '2 · Burial:' }), 'Cast 50 rays per voxel (spherical Fabonacci distribution), then count the fraction 13 Å blocked by receptor atoms (probe 1.8 Å)within the specified range.'),
      div(null, span({ class: 'muted', text: '3 · Retention:' }), 'Keep only voxels with burial ≥ 0.5 ; their summed volume gives 5,440 Å³. '),
      div(null, span({ class: 'muted', text: '4 · Chemical labeling:' }), 'Assign each voxel the nearest MurJ residue type; assign subpockets by nearest residue centroid.'),
      div({ class: 'card-note', text: 'These geometric quantities describe cavity size, depth and chemistry; they do not independently establish druggability.' })));

  root.appendChild(div({ class: 'section grid wide-left' },
    card({ title: 'Subpocket data' }, table.el),
    method));

  return null;
}

function dominant(chemistry) {
  if (!chemistry) return '—';
  const best = CHEM.map((c) => ({ ...c, n: chemistry[c.key] || 0 })).sort((a, b) => b.n - a.n)[0];
  return span({ class: 'split' },
    span({ class: 'swatch', style: { background: best.color, width: '8px', height: '8px', borderRadius: '2px' } }),
    span({ class: 'tiny', text: best.label }));
}
