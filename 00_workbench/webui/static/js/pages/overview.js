/* pages/overview.js -- one-screen summary of the whole design pipeline. */

import { api } from '../api.js';
import {
  el, div, span, card, stat, chip, barList, funnel, legend, kvList,
  fmt, fmtInt, fmtDate, timeAgo, loadingPane, statusChip, markdownNode,
} from '../ui.js';

function pctOf(value, total) {
  if (!total) return '—';
  return `${((value / total) * 100).toFixed(value / total < 0.1 ? 1 : 0)}%`;
}

export async function mount(root, ctx) {
  root.appendChild(loadingPane('Loading pipeline results...'));
  const [overview, meta] = await Promise.all([api.overview(), api.meta()]);
  const summary = overview.summary || {};
  const fun = overview.funnel || {};
  const filters = overview.filters || {};
  const ev = overview.evidence || {};
  const best = overview.best || [];
  const hot = overview.hotspots || [];
  const feats = overview.features || [];
  const axes = overview.axes || [];
  const pocket = overview.pocket || {};

  const nDrug = fun.drug_like ?? 0;
  const nKeep = filters.survivors ?? 0;
  const nEnumerated = fun.enumerated ?? 0;
  const bestTotal = best.length ? best[0].total : null;

  const typeCount = {};
  for (const f of feats) typeCount[f.type] = (typeCount[f.type] || 0) + 1;

  root.replaceChildren();

  root.appendChild(div({ class: 'page-head' },
    div(null,
      el('h1', { text: 'Design pipeline overview' }),
      div({ class: 'sub', text: 'Starting from three cryo-EM structures (9NU4 / 9NU5 / 9NU8)), translate convergent phage-peptide inhibition hypotheses into enumerable, interpretable scaffold candidates. Geometric values are calculated locally from coordinates; synthesis feasibility requires separate assessment.' })),
    div({ class: 'actions' },
      el('button', { class: 'btn', text: '3D Observation', onclick: () => ctx.navigate('structure', { ligand: '0' }) }),
      el('button', { class: 'btn primary', text: 'Run pipeline', onclick: () => ctx.navigate('run') }))));

  /* ---- status strip ---- */
  const byStatus = ev.by_status || {};
  const statusChips = Object.entries(byStatus)
    .sort((a, b) => b[1] - a[1])
    .map(([key, n]) => span({ class: 'split' }, statusChip(key), span({ class: 'tiny muted', text: `×${n}` })));

  root.appendChild(div({ class: 'grid cols-auto', style: { marginBottom: '12px' } },
    card({ class: 'tight' },
      div({ class: 'split' },
        span({ class: 'tiny muted', text: 'Evidence checks' }),
        span({ class: 'right tiny muted', style: { marginLeft: 'auto' }, text: timeAgo(ev.generated || 0) })),
      div({ class: 'split', style: { marginTop: '6px' } }, ...statusChips)),
    card({ class: 'tight' },
      div({ class: 'tiny muted', text: 'Structure filters (F1–F4)' }),
      div({ class: 'split', style: { marginTop: '6px' } },
        chip(`${filters.n_input ?? 0} Input`, 'info'),
        span({ class: 'tiny muted', text: '→' }),
        chip(`${nKeep} Retained`, 'pass'),
        span({ class: 'tiny muted', text: `Maximum score difference across structures ${filters.max_spread ?? '—'} Å·score` }))),
    card({ class: 'tight' },
      div({ class: 'tiny muted', text: 'Structural input' }),
      div({ class: 'split', style: { marginTop: '6px' } },
        ...['9NU4', '9NU5', '9NU8'].map((id) => chip(id, 'mono info')),
        span({ class: 'tiny muted', text: '3.6–3.7 Å' }))),
    card({ class: 'tight' },
      div({ class: 'tiny muted', text: 'Result generation time' }),
      div({ class: 'mono-sm', style: { marginTop: '6px' }, text: fmtDate(Math.max(...Object.values(meta.files || { x: 0 }))) }))));

  /* ---- KPI row ---- */
  root.appendChild(div({ class: 'grid cols-4' },
    card(null, stat({
      label: 'Enumerated scaffolds', value: fmtInt(nEnumerated), unit: 'items',
      foot: `${axes.length} design axes · Within property envelope ${fmtInt(fun.property_pass)} entered pose scoring`,
    })),
    card(null, stat({
      label: 'Three checks passed', value: fmtInt(nDrug), unit: 'items', tone: 'accent',
      foot: `No steric clash ${fmtInt(fun.no_clash)} · Peptide mimicry ${fmtInt(fun.peptide_mimic)}`,
    })),
    card(null, stat({
      label: 'Structure filters passed', value: fmtInt(nKeep), unit: 'items', tone: 'good',
      foot: 'Anchor + Cross-subpocket + Three-structure consistency (F4 is the main bottleneck)',
    })),
    card(null, stat({
      label: 'Best scaffold total score', value: bestTotal === null ? '—' : fmt(bestTotal, 2), tone: 'warn',
      foot: best.length ? `${best[0].props?.MW} Da · cLogP ${best[0].props?.cLogP} · TPSA ${best[0].props?.TPSA}` : '',
    }))));

  /* ---- funnel + best ---- */
  const funnelCard = card({
    title: 'Screening funnel',
    right: span({ class: 'tiny muted', text: 'Deterministic enumeration along one axis' }),
  },
    funnel([
      { label: 'Scaffold enumeration', value: nEnumerated, hint: 'Exact length-constraint solution' },
      { label: 'Property envelope passed', value: fun.property_pass ?? 0, hint: 'MW ≤ 500, cLogP/TPSA Envelope' },
      { label: 'Completed 3D Pose', value: fun.posed ?? 0, hint: 'Rigid placement + 12 rotations' },
      { label: 'No steric clashes', value: fun.no_clash ?? 0, hint: 'clash < 0.25' },
      { label: 'Peptide mimicry threshold passed', value: fun.peptide_mimic ?? 0, hint: 'Spatial overlap with three peptides ≥ 0.5' },
      { label: 'Three checks passed', value: nDrug, hint: 'Steric clash + Drug-like + Peptide mimicry', final: true },
    ]),
    div({ class: 'hr' }),
    div({ class: 'split tiny muted' },
      span({ text: `Structure-derived filters:${filters.n_input ?? 0} → ` }),
      span({ style: { color: 'var(--accent-2)' }, text: `${nKeep}` }),
      span({ text: `(Anchor ${filters.by_stage?.anchor ?? '—'} · Bridging ${filters.by_stage?.bridge ?? '—'} · Three-structure agreement ${filters.by_stage?.ensemble ?? '—'})` })));

  const bestCard = card({
    title: 'Top ten scaffolds',
    right: el('button', { class: 'btn sm ghost', text: 'View all →', onclick: () => ctx.navigate('scaffolds') }),
  },
    div({ class: 'list-reset', style: { display: 'grid', gap: '6px' } },
      ...best.slice(0, 10).map((row, idx) => div({
        class: 'split', style: {
          padding: '6px 8px', borderRadius: '8px', background: idx === 0 ? 'rgba(56,189,248,.07)' : 'transparent',
          border: '1px solid ' + (idx === 0 ? 'rgba(56,189,248,.25)' : 'transparent'), cursor: 'pointer',
        },
        onclick: () => ctx.navigate('scaffolds', { focus: row.smiles }),
      },
        span({ class: 'mono-sm', style: { color: 'var(--text-mute)', width: '22px' }, text: `${idx + 1}.` }),
        el('code', { style: { flex: '1', wordBreak: 'break-all', fontSize: '11.5px' }, text: row.smiles }),
        span({ class: 'tiny muted nowrap', text: `MW ${row.props?.MW} · cLogP ${row.props?.cLogP}` }),
        span({ class: 'chip mono', style: { color: 'var(--accent)' }, text: fmt(row.total, 2) })))));

  root.appendChild(div({ class: 'section grid wide-left' }, funnelCard, bestCard));

  /* ---- hotspots + pharmacophore ---- */
  const hotCard = card({
    title: 'Interface hotspots (top 10)',
    right: el('button', { class: 'btn sm ghost', text: 'All 52 items →', onclick: () => ctx.navigate('hotspots') }),
  },
    barList(hot.slice(0, 10).map((h) => ({
      label: `${h.resname}${h.murj_resid}`,
      value: h.score,
      color: h.evidence ? '#ff8a5b' : '#3f7fb5',
      title: h.evidence || `${h.resname}${h.murj_resid}`,
    })), { format: (v) => fmt(v, 3) }),
    div({ class: 'hr' }),
    legend([
      { color: '#ff8a5b', label: 'Genetic evidence and/resistance evidence' },
      { color: '#3f7fb5', label: 'Geometry only' },
    ]));

  const pharmCard = card({
    title: 'Pharmacophore consensus',
    right: el('button', { class: 'btn sm ghost', text: 'Details →', onclick: () => ctx.navigate('pharmacophore') }),
  },
    div({ class: 'grid cols-2', style: { gap: '8px' } },
      stat({ label: 'Consensus features', value: feats.length, unit: 'items', small: true }),
      stat({ label: 'Design axis', value: axes.length, unit: 'entries', small: true })),
    div({ class: 'hr' }),
    div({ class: 'chips' },
      ...Object.entries(typeCount).map(([type, n]) => span({
        class: 'chip', style: { borderColor: 'rgba(255,255,255,.12)' },
      }, span({ class: 'swatch', style: { background: FEATURE_COLOR[type] || '#8b97a8' } }), `${type} ×${n}`))),
    div({ class: 'card-note', style: { marginTop: '8px' }, text: 'Features are clustered from functional atoms in the shared frame of three peptide structures;weight weighted by convergence.' }));

  root.appendChild(div({ class: 'section grid cols-2' }, hotCard, pharmCard));

  /* ---- pocket + figures ---- */
  const subs = pocket.subpockets || [];
  const pocketCard = card({
    title: 'Buried cavity and subpockets',
    right: el('button', { class: 'btn sm ghost', text: 'Pocket details →', onclick: () => ctx.navigate('pocket') }),
  },
    div({ class: 'split', style: { marginBottom: '8px' } },
      stat({ label: 'Cavity volume', value: fmtInt(pocket.total_volume_A3 || 0), unit: 'Å³', small: true }),
      stat({ label: 'Subpocket', value: subs.length, unit: 'items', small: true })),
    div({ class: 'bar-list' },
      ...subs.map((sp) => div({ class: 'row' },
        div({ class: 'k', text: sp.subpocket.replace(/^S\d_/, ''), title: sp.subpocket }),
        div({ class: 'track' },
          div({ class: 'fill', style: { width: `${(sp.volume_A3 / Math.max(...subs.map((x) => x.volume_A3), 1)) * 100}%`, background: SUBPOCKET_COLOR[sp.subpocket] || '#8b97a8' } })),
        div({ class: 'num', text: fmtInt(sp.volume_A3) })))),
    div({ class: 'card-note', style: { marginTop: '8px' }, text: 'Volume = Voxels × 0.8³ Å³; S2 Carboxylate gate burial 0.94, is the deepest site.' }));

  const fileCard = card({ title: 'Artifact inventory' },
    div({ class: 'kv tiny' },
      ...Object.entries(meta.files || {}).map(([name, mtime]) => [
        el('div', { class: 'k mono-sm', style: { wordBreak: 'break-all' }, text: name }),
        div({ class: 'v tiny muted', text: mtime ? timeAgo(mtime) : 'Missing' }),
      ]).flat()));

  root.appendChild(div({ class: 'section grid cols-2' }, pocketCard, fileCard));

  const figs = (meta.figures || []).filter((f) => f.name.endsWith('.png'));
  if (figs.length) {
    root.appendChild(card({ title: 'Pipeline diagram', class: 'section' },
      div({ class: 'figure-grid' },
        ...figs.slice(0, 4).map((f) => el('figure', null,
          el('img', { src: `/files/figures/${f.name}`, alt: f.name, loading: 'lazy' }),
          el('figcaption', null,
            el('a', { href: `/files/figures/${f.name}`, target: '_blank', text: f.name }),
            span({ class: 'muted', text: ` · ${(f.size / 1024).toFixed(0)} KB` })))))));
  }

  return null;
}

const FEATURE_COLOR = {
  anion: '#e55353', cation: '#4f8ef7', aromatic: '#c99a3a', donor: '#38a169', hydrophobic: '#8b97a8',
};
const SUBPOCKET_COLOR = {
  S1_basic_bowl: '#4f8ef7', S2_carboxylate_gate: '#e55353',
  S3_TM2_TM7_groove: '#38a169', S4_TM7_floor: '#c99a3a',
};
