/* pages/scaffolds.js -- the ranked scaffold library: search, filter, inspect. */

import { api } from '../api.js';
import {
  el, div, span, card, stat, chip, DataTable, kvList, meterRow, toast,
  fmt, fmtInt, loadingPane, openDrawer, downloadText, toCSV, formatNow, escapeHtml,
} from '../ui.js';

const TRIPLE = (row) => (row.property_penalty ?? 9) <= 0.35
  && (row.peptide_overlap ?? 0) >= 0.5
  && (row.clash ?? 9) < 0.25;

export async function mount(root, ctx) {
  root.appendChild(loadingPane('Loading 200 ranked scaffolds...'));
  const data = await api.scaffolds();
  const modes = await api.poses().catch(() => ({ top: [] }));
  const rows = data.rows || [];
  const poseBySmiles = new Map((modes.top || []).map((m, i) => [m.smiles, i]));
  for (const row of rows) {
    if (row.pose_index === null || row.pose_index === undefined) row.pose_index = poseBySmiles.get(row.smiles) ?? null;
    row.triple = TRIPLE(row);
  }

  const axes = [...new Set(rows.map((r) => r.axis).filter(Boolean))].sort();
  const state = {
    q: '',
    axis: '',
    onlyTriple: false,
    onlyKeep: false,
    maxMw: 500,
  };

  const table = new DataTable({
    columns: [
      { key: 'rank', label: '#', num: true, width: '48px' },
      { key: 'smiles', label: 'SMILES', mono: true, html: true, sortable: false, width: '300px', render: (r) => `<span style="color:#ffd7a8">${escapeHtml(r.smiles)}</span>` },
      { key: 'axis', label: 'Design axis', width: '76px' },
      { key: 'MW', label: 'MW', num: true, render: (r) => fmt(r.MW, 1) },
      { key: 'cLogP', label: 'cLogP', num: true, render: (r) => fmt(r.cLogP, 2) },
      { key: 'TPSA', label: 'TPSA', num: true, render: (r) => fmt(r.TPSA, 0) },
      { key: 'heavy_atoms', label: 'Heavy atoms', num: true },
      { key: 'peptide_overlap', label: 'Peptide mimicry', num: true, render: (r) => fmt(r.peptide_overlap, 2) },
      { key: 'clash', label: 'Steric clash', num: true, render: (r) => fmt(r.clash, 2) },
      {
        key: 'total', label: 'Total score', num: true, render: (r) => r.total,
      },
      {
        key: 'f_keep', label: 'Structure filters', sortable: false,
        render: (r) => (r.f_keep === true ? chip('Retained', 'pass') : r.f_keep === false ? chip('Not passed', 'weak') : chip('—')),
      },
      {
        key: 'n_contacts', label: 'Contact residues', sortable: false,
        render: (r) => (r.murj_contacts ? String(r.murj_contacts).split(/\s+/).length : 0),
      },
    ],
    rows,
    keyField: 'smiles',
    pageSize: 50,
    initialSort: { key: 'total', dir: 'desc' },
    empty: 'No scaffolds match the current filters',
    onRowClick: (row) => openDetail(row),
  });

  /* ---- toolbar ---- */
  const search = el('input', {
    type: 'search', class: 'search', placeholder: 'Search SMILES / Design axis / Contact residue (e.g. K368)',
    oninput: (event) => { state.q = event.target.value.trim(); apply(); },
  });
  const axisSel = el('select', { onchange: (event) => { state.axis = event.target.value; apply(); } },
    el('option', { value: '' }, 'All design axes'),
    ...axes.map((axis) => el('option', { value: axis }, axis)));
  const tripleBox = el('input', {
    type: 'checkbox', onchange: (event) => { state.onlyTriple = event.target.checked; apply(); },
  });
  const keepBox = el('input', {
    type: 'checkbox', onchange: (event) => { state.onlyKeep = event.target.checked; apply(); },
  });
  const mwInput = el('input', {
    type: 'number', min: '150', max: '600', step: '10', value: '500', style: { width: '74px' },
    onchange: (event) => { state.maxMw = Number(event.target.value) || 500; apply(); },
  });
  const countLabel = span({ class: 'tiny muted' });

  const toolbar = div({ class: 'toolbar' },
    search,
    axisSel,
    label(tripleBox, 'Only triple-pass scaffolds'),
    label(keepBox, 'Only structure-filter survivors'),
    label(mwInput, 'MW ≤'),
    div({ class: 'grow' }),
    countLabel,
    el('button', { class: 'btn sm', text: 'Export current view CSV', onclick: exportCsv }));

  function label(input, text) {
    return el('label', { class: 'check' }, input, text);
  }

  function apply() {
    const q = state.q.toLowerCase();
    const filtered = rows.filter((row) => {
      if (state.axis && row.axis !== state.axis) return false;
      if (state.onlyTriple && !row.triple) return false;
      if (state.onlyKeep && row.f_keep !== true) return false;
      if (state.maxMw && Number(row.MW) > state.maxMw) return false;
      if (!q) return true;
      const haystack = `${row.smiles} ${row.axis} ${row.family} ${row.murj_contacts}`.toLowerCase();
      return haystack.includes(q);
    });
    table.setRows(filtered);
    countLabel.textContent = `${filtered.length} / ${rows.length} entries`;
  }

  function exportCsv() {
    const columns = [
      { key: 'rank', label: 'rank' }, { key: 'smiles', label: 'smiles' },
      { key: 'axis', label: 'axis' }, { key: 'family', label: 'family' },
      { key: 'MW', label: 'MW' }, { key: 'cLogP', label: 'cLogP' },
      { key: 'TPSA', label: 'TPSA' }, { key: 'HBD', label: 'HBD' },
      { key: 'HBA', label: 'HBA' }, { key: 'RotB', label: 'RotB' },
      { key: 'heavy_atoms', label: 'heavy_atoms' },
      { key: 'total', label: 'total' }, { key: 'pose_score', label: 'pose_score' },
      { key: 'peptide_overlap', label: 'peptide_overlap' }, { key: 'clash', label: 'clash' },
      { key: 'property_penalty', label: 'property_penalty' },
      { key: 'murj_contacts', label: 'murj_contacts', csv: (r) => String(r.murj_contacts || '').replace(/\s+/g, ';') },
      { key: 'f_keep', label: 'structure_filter_pass' },
    ];
    downloadText(`wedge-scaffolds-${formatNow()}.csv`, toCSV(columns, table.rows), 'text/csv;charset=utf-8');
    toast(`Exported ${table.rows.length} rows`, 'ok');
  }

  /* ---- detail drawer ---- */
  function openDetail(row) {
    const poseIndex = row.pose_index === null || row.pose_index === undefined ? null : Number(row.pose_index);
    const contacts = String(row.murj_contacts || '').split(/\s+/).filter(Boolean);
    const f = row.f_scores || {};
    const filters = row.f_keep === undefined
      ? div({ class: 'card-note', text: 'This scaffold was not evaluated by the structure-derived filters; only the triple-pass subset was assessed.' })
      : div({ style: { display: 'grid', gap: '8px' } },
        div({ class: 'split' },
          chip(row.f_keep ? 'F1–F4 All passed' : 'Structure filters failed', row.f_keep ? 'pass' : 'weak'),
          row.f_anchor_ids ? chip(`Anchor ${row.f_anchor_ids}`, 'info') : null,
          chip(`Bridging ${row.f_bridge ? 'represent' : 'No'}`, row.f_bridge ? 'pass' : ''),
          chip(`Cross-structure score difference ${fmt(row.f_ens_spread, 2)}`, Math.abs(row.f_ens_spread ?? 99) <= 5 ? 'pass' : 'weak')),
        div({ class: 'kv tiny' },
          el('div', { class: 'k', text: 'Scores across three structures' }),
          el('div', { class: 'v mono-sm', text: Object.entries(f).map(([k, v]) => `${k.replace('Sgl', 'Sgl-')} ${fmt(v, 2)}`).join('   ') }),
          ...(row.f_ens_rmsd !== undefined ? [
            el('div', { class: 'k', text: 'Replay RMSD' }),
            el('div', { class: 'v mono-sm', text: `${fmt(row.f_ens_rmsd, 2)} Å` }),
          ] : [])));

    const body = div({ style: { display: 'grid', gap: '14px' } },
      div({ class: 'smiles-box', text: row.smiles }),
      div({ class: 'split' },
        chip(`Rank #${row.rank}`, 'info'),
        chip(`Axis ${row.axis}`, 'mono'),
        chip(`Family ${row.family}`, 'mono'),
        row.triple ? chip('Three checks passed', 'pass') : null),
      div({ class: 'grid cols-2', style: { gap: '10px' } },
        div({ class: 'card flat' },
          div({ class: 'card-head' }, el('h3', { text: 'Property envelope' })),
          kvList([
            ['Molecular weight MW', `${fmt(row.MW, 1)} Da`],
            ['cLogP', fmt(row.cLogP, 2)],
            ['TPSA', `${fmt(row.TPSA, 1)} Å²`],
            ['H-bond donors / acceptors', `${row.HBD} / ${row.HBA}`],
            ['Rotatable bonds', row.RotB],
            ['Ring count', row.Rings],
            ['Heavy atoms', row.heavy_atoms],
          ])),
        div({ class: 'card flat' },
          div({ class: 'card-head' }, el('h3', { text: 'Score breakdown' })),
          div({ style: { display: 'grid', gap: '6px' } },
            meterRow('Total score', row.total, 16, 'var(--accent)'),
            meterRow('Pose score', row.pose_score, 12, 'var(--accent-2)'),
            meterRow('Peptide mimicry', row.peptide_overlap, 1, '#a78bfa'),
            meterRow('Burial', row.buried, 1, '#38bdf8'),
            meterRow('Rigidity', row.rigidity, 1, '#38a169'),
            meterRow('Salt bridge', row.salt, 3, '#4f8ef7'),
            meterRow('Hydrogen bond', row.hb, 12, '#c99a3a'),
            meterRow('Hydrophobic', row.hyd, 25, '#8b97a8'),
            meterRow('Steric penalty', row.clash, 1, 'var(--danger)'),
            meterRow('Property penalty', row.property_penalty, 1, 'var(--warn)')))),

      div({ class: 'card flat' },
        div({ class: 'card-head' }, el('h3', { text: 'Structure-derived filters (F1–F4)' })),
        filters),

      div({ class: 'card flat' },
        div({ class: 'card-head' },
          el('h3', { text: 'Contacting MurJ Residue' }),
          span({ class: 'right tiny muted', text: `${contacts.length} items,4.2 Å within` })),
        div({ class: 'chips' }, ...contacts.map((c) => span({ class: 'resid', text: c })))),

      div({ class: 'card flat' },
        div({ class: 'card-head' }, el('h3', { text: 'Construction' })),
        kvList([
          ['Head cap', row.head], ['Ring system', row.rings], ['Linker', row.links], ['Tail cap', row.tail],
          ['Design axis length', `${fmt(row.axis_length_A, 2)} Å`],
        ])));

    const footer = [
      el('button', {
        class: 'btn primary',
        text: poseIndex === null ? 'This scaffold has no 3D Pose' : 'in 3D view pose',
        disabled: poseIndex === null,
        onclick: () => {
          if (poseIndex === null) return;
          close();
          ctx.navigate('structure', { ligand: String(poseIndex) });
        },
      }),
      el('button', {
        class: 'btn', text: 'Copy SMILES',
        onclick: async () => {
          try {
            await navigator.clipboard.writeText(row.smiles);
            toast('SMILES Copied', 'ok', 1500);
          } catch (err) {
            toast('Clipboard access denied; select and copy manually', 'err');
          }
        },
      }),
      el('button', {
        class: 'btn ghost', text: 'Download row JSON',
        onclick: () => downloadText(`scaffold-${row.rank}.json`, JSON.stringify(row, null, 2), 'application/json'),
      }),
    ];

    const { close } = openDrawer({
      title: `Scaffold #${row.rank}`,
      sub: `${row.axis} · ${row.family} · Total score ${fmt(row.total, 3)}`,
      body, footer,
    });
  }

  /* ---- page ---- */
  const nTriple = rows.filter((r) => r.triple).length;
  const nKeepHere = rows.filter((r) => r.f_keep === true).length;
  const fsum = data.filters || {};
  const nKeepAll = fsum.survivors ?? null;
  const nPoses = data.poses ?? 0;

  root.replaceChildren();
  root.appendChild(div({ class: 'page-head' },
    div(null,
      el('h1', { text: 'Scaffold library' }),
      div({ class: 'sub', text: `Top scaffolds ranked by total score: ${rows.length} scaffolds, including ${nTriple} passing the triple checks (sterics, + Drug-like + peptide mimicry); of all ${fsum.n_input ?? 538} triple-pass scaffolds, ${nKeepAll ?? 32} passed structure-derived filters; this page covers ${nKeepHere} . Poses are rigid placements, not docked conformations; downstream MD/MM-PBSA refinement is required.` })),
    div({ class: 'actions' },
      el('a', { class: 'btn', href: '/files/results/scaffolds_ranked.csv?download=1', text: 'Download CSVs' }),
      el('a', { class: 'btn', href: '/files/results/top_scaffold_poses.sdf?download=1', text: 'Download SDF' }))));

  root.appendChild(div({ class: 'grid cols-4', style: { marginBottom: '12px' } },
    card({ class: 'tight' }, stat({ label: 'Scaffolds on this page', value: rows.length, unit: 'items', small: true, foot: 'Ranked by total score' })),
    card({ class: 'tight' }, stat({ label: 'Triple-pass subset', value: nTriple, unit: 'items', small: true, tone: 'accent', foot: 'Steric clash + Drug-like + Peptide mimicry' })),
    card({ class: 'tight' }, stat({
      label: 'Structure filters passed', value: nKeepAll ?? '—', unit: 'items', small: true, tone: 'good',
      foot: `Among all triple-pass scaffolds; on this page: ${nKeepHere} items`,
    })),
    card({ class: 'tight' }, stat({ label: 'Available 3D Pose', value: nPoses, unit: 'items', small: true, foot: 'top_scaffold_poses.sdf' }))));

  root.appendChild(card({ title: 'Ranking and filters' }, toolbar, table.el));

  apply();

  if (ctx.params.focus) {
    const hit = rows.find((row) => row.smiles === ctx.params.focus);
    if (hit) {
      table.select(hit.smiles);
      openDetail(hit);
    }
  }

  return null;
}
