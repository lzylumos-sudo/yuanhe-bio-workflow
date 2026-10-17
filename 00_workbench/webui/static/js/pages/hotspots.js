/* pages/hotspots.js -- MurJ interface hotspots and the per-structure contacts. */

import { api } from '../api.js';
import {
  el, div, span, card, stat, chip, barList, DataTable, kvList, legend,
  fmt, fmtInt, loadingPane, statusChip,
} from '../ui.js';

export async function mount(root, ctx) {
  root.appendChild(loadingPane('Summarizing interface residues...'));
  const [hot, iface] = await Promise.all([api.hotspots(), api.evidenceInterface()]);
  const rows = (hot.rows || []).map((r) => ({
    ...r,
    label: `${r.resname}${r.murj_resid}`,
    has_evidence: !!r.evidence,
  })).sort((a, b) => Number(b.score) - Number(a.score));
  const nCommon = rows.filter((r) => Number(r.convergence) >= 1).length;
  const nEvidence = rows.filter((r) => r.has_evidence).length;
  const interfaceRows = (iface.murj || []);

  root.replaceChildren();
  root.appendChild(div({ class: 'page-head' },
    div(null,
      el('h1', { text: 'Interface hotspots' }),
      div({ class: 'sub', text: 'Contact surfaces of three independently evolved phage peptides on MurJ .score = 0.30·Cross-peptide convergence + 0.25·Burial + 0.15·Contact density + 0.20·Functionally required + 0.10·Resistance evidence; all geometric values are recomputable, not manually entered.' })),
    div({ class: 'actions' },
      el('a', { class: 'btn', href: '/files/results/interface_hotspots.csv?download=1', text: 'Download CSV' }),
      el('button', { class: 'btn', text: 'in 3D view', onclick: () => ctx.navigate('structure', { ligand: '0' }) }))));

  root.appendChild(div({ class: 'grid cols-4', style: { marginBottom: '12px' } },
    card({ class: 'tight' }, stat({ label: 'Interface residues', value: rows.length, unit: 'items', small: true })),
    card({ class: 'tight' }, stat({ label: 'Contacted by all three peptides', value: nCommon, unit: 'items', small: true, tone: 'accent' })),
    card({ class: 'tight' }, stat({ label: 'Genetic evidence available', value: nEvidence, unit: 'items', small: true, tone: 'warn' })),
    card({ class: 'tight' }, stat({
      label: 'Maximum dSASA', unit: 'Å²', small: true,
      value: fmt(Math.max(...rows.map((r) => Number(r.dSASA_A2) || 0)), 0),
    }))));

  /* ---- chart + common residues ---- */
  const chartCard = card({ title: 'Hotspot ranking (top 14)' },
    barList(rows.slice(0, 14).map((r) => ({
      label: r.label, value: r.score, title: r.evidence || 'Geometric evidence',
      color: r.has_evidence ? '#ff8a5b' : '#3f7fb5',
    })), { format: (v) => fmt(v, 3) }),
    div({ class: 'hr' }),
    card({ class: 'flat tight' },
      div({ class: 'tiny muted', style: { marginBottom: '6px' } },
        'All three peptides in contact (geometric convergence = 1.0, n = ', String(nCommon), ')'),
      div({ class: 'chips' },
        ...rows.filter((r) => Number(r.convergence) >= 1).map((r) => span({ class: 'resid', text: `${r.resname}${r.murj_resid}` })))),
    div({ class: 'card-note', style: { marginTop: '8px' }, text: 'Convergence here means contact in all three structures, not sequence conservation. These quantities differ.' }));

  /* ---- detail panel ---- */
  const detail = div({ class: 'card' },
    div({ class: 'card-head' }, el('h3', { text: 'Residue details' })),
    div({ class: 'card-note', text: 'Click a residue for its contact counts, minimum distances and buried areas across three structures.' }));

  function showDetail(row) {
    const perComplex = interfaceRows.filter((r) => Number(r.resid) === Number(row.murj_resid));
    detail.replaceChildren(
      div({ class: 'card-head' },
        el('h3', { text: `${row.resname}${row.murj_resid}` }),
        span({ class: 'right chip mono', text: `score ${fmt(row.score, 3)}` })),
      div({ class: 'grid cols-2', style: { gap: '8px', marginBottom: '10px' } },
        stat({ label: 'Convergence', value: fmt(row.convergence, 2), small: true, foot: 'Fraction of three structures' }),
        stat({ label: 'Buried area', value: fmt(row.dSASA_A2, 0), unit: 'Å²', small: true })),
      kvList([
        ['Complex', (row.complexes || '').split(';').filter(Boolean).join(' · ') || '—'],
        ['Contact atom pairs', fmtInt(row.n_contacts)],
        ['Contacting peptide residues', fmtInt(row.n_peptide_residues)],
        ['Nearest heavy-atom distance', row.dmin_A ? `${fmt(row.dmin_A, 2)} Å` : '—'],
      ]),
      row.evidence ? div({ class: 'card flat tight', style: { marginTop: '10px' } },
        div({ class: 'tiny muted', text: 'Experimental evidence' }),
        div({ class: 'tiny', style: { marginTop: '2px' }, text: row.evidence })) : null,
      div({ class: 'hr' }),
      div({ class: 'tiny muted', style: { marginBottom: '4px' }, text: 'Measurements by structure (evidence/01_interface_murj.csv)' }),
      div({ class: 'table-wrap', style: { border: 'none' } },
        el('table', null,
          el('thead', null, el('tr', null,
            ...['Complex', 'Contacting atoms', 'Peptide residues', 'dmin Å', 'dSASA Å²'].map((t, i) => el('th', { class: i ? 'num no-sort' : 'no-sort', text: t })))),
          el('tbody', null,
            ...(perComplex.length ? perComplex : [{ complex: '—' }]).map((r) => el('tr', null,
              el('td', { text: r.complex || '—' }),
              el('td', { class: 'num', text: fmtInt(r.n_atom_contacts) }),
              el('td', { class: 'num', text: fmtInt(r.n_peptide_residues) }),
              el('td', { class: 'num', text: fmt(r.dmin_A, 2) }),
              el('td', { class: 'num', text: fmt(r.dSASA_A2, 1) })))))),
      div({ class: 'split', style: { marginTop: '10px' } },
        el('button', {
          class: 'btn sm', text: 'Locate in pocket →',
          onclick: () => ctx.navigate('structure', { ligand: '0', focus: String(row.murj_resid) }),
        })));
    table.select(row.murj_resid);
  }

  /* ---- table ---- */
  const table = new DataTable({
    columns: [
      { key: 'murj_resid', label: 'Residue number', num: true, width: '70px' },
      { key: 'label', label: 'Residue', width: '60px' },
      { key: 'score', label: 'score', num: true, render: (r) => fmt(r.score, 3) },
      { key: 'convergence', label: 'Convergence', num: true, render: (r) => fmt(r.convergence, 2) },
      { key: 'dSASA_A2', label: 'dSASA Å²', num: true, render: (r) => fmt(r.dSASA_A2, 1) },
      { key: 'n_contacts', label: 'Contacting atoms', num: true },
      { key: 'n_peptide_residues', label: 'Peptide residues', num: true },
      { key: 'dmin_A', label: 'dmin Å', num: true, render: (r) => fmt(r.dmin_A, 2) },
      { key: 'complexes', label: 'Complex', mono: true },
      {
        key: 'evidence', label: 'Experimental evidence', sortable: false, wrap: true,
        render: (r) => (r.evidence ? chip(r.evidence.length > 56 ? `${r.evidence.slice(0, 56)}…` : r.evidence, 'warn') : span({ class: 'muted tiny', text: '—' })),
      },
    ],
    rows,
    keyField: 'murj_resid',
    pageSize: 25,
    initialSort: { key: 'score', dir: 'desc' },
    onRowClick: showDetail,
  });

  root.appendChild(div({ class: 'section grid wide-left' }, chartCard, detail));
  root.appendChild(card({ title: 'All interface residues', class: 'section' }, table.el));
  showDetail(rows[0]);

  return null;
}
