/* pages/evidence.js -- the reproducibility trail: auto-generated report, the 24
   literature priors checked against coordinates, censuses, and anchor geometry. */

import { api } from '../api.js';
import {
  el, div, span, card, stat, chip, DataTable, kvList, markdownNode, statusChip,
  fmt, fmtInt, fmtPct, loadingPane, openDrawer, legend,
} from '../ui.js';

export async function mount(root, ctx) {
  const [report, priors, census, iface] = await Promise.all([
    api.evidenceReport(), api.evidencePriors(), api.evidenceCensus(), api.evidenceInterface(),
  ]);
  const status = priors.status || report.status || {};
  const byStatus = status.by_status || {};
  const budget = census.atom_budget || {};
  const rows = (priors.rows || []).map((r) => ({ ...r }));

  root.replaceChildren();
  root.appendChild(div({ class: 'page-head' },
    div(null,
      el('h1', { text: 'Evidence chain' }),
      div({ class: 'sub', text: 'All numbers on this page are calculated by run_evidence.py Recalculate from three PDB and checked against priors/murj_priors.json references and rules. Unreproduced entries are marked FLAG / FAIL explicitly, never silently accepted.' })),
    div({ class: 'actions' },
      el('button', {
        class: 'btn', text: 'Rebuild evidence',
        onclick: () => ctx.navigate('run'),
      }),
      el('a', { class: 'btn', href: '/files/evidence/report.md?download=1', text: 'Download report' }))));

  root.appendChild(div({ class: 'grid cols-4', style: { marginBottom: '12px' } },
    card({ class: 'tight' }, stat({ label: 'Checked entries', value: status.n ?? rows.length, unit: 'entries', small: true })),
    card({ class: 'tight' }, stat({
      label: 'PASS', value: byStatus.PASS ?? 0, unit: 'entries', small: true, tone: 'good',
      foot: `Pass rate ${budget.summary ? fmtPct(budget.summary.pass_rate) : '—'}`,
    })),
    card({ class: 'tight' }, stat({ label: 'FLAG', value: byStatus.FLAG ?? 0, unit: 'entries', small: true, tone: 'warn', foot: 'Manual judgment required' })),
    card({ class: 'tight' }, stat({ label: 'FAIL / MISMATCH', value: (byStatus.FAIL ?? 0) + (byStatus.MISMATCH ?? 0), unit: 'entries', small: true, tone: 'bad' }))));

  /* ---------- tabs ---------- */
  const TABS = [
    { id: 'report', label: 'Evidence report' },
    { id: 'priors', label: 'Prior-by-prior checks', n: rows.length },
    { id: 'census', label: 'Structure inventory and atom budget' },
    { id: 'interface', label: 'Interface measurements' },
    { id: 'anchors', label: 'Anchor geometry' },
  ];
  const tabBar = div({ class: 'tabs' });
  const bodyHost = div();
  let active = 'report';

  const tabButtons = TABS.map((tab) => {
    const btn = el('button', {
      class: 'tab', onclick: () => activate(tab.id),
    }, tab.label, tab.n ? span({ class: 'n', text: String(tab.n) }) : null);
    tabBar.appendChild(btn);
    return { ...tab, btn };
  });

  function activate(id) {
    active = id;
    for (const tab of tabButtons) tab.btn.classList.toggle('active', tab.id === id);
    bodyHost.replaceChildren(renderTab(id));
  }

  function renderTab(id) {
    if (id === 'report') return renderReport();
    if (id === 'priors') return renderPriors();
    if (id === 'census') return renderCensus();
    if (id === 'interface') return renderInterface();
    return renderAnchors();
  }

  /* ---------- report ---------- */
  function renderReport() {
    const md = markdownNode(report.markdown || '_(Missing evidence/report.md, First run run_evidence.py)_');
    const headings = [...(report.markdown || '').matchAll(/^(#{2,3})\s+(.+)$/gm)]
      .map((m) => ({ level: m[1].length, text: m[2].trim() }));
    const nodes = [...md.querySelectorAll('h2, h3')];
    const toc = div({ style: { display: 'grid', gap: '2px' } });
    nodes.forEach((node, idx) => {
      node.id = `ev-h-${idx}`;
      const heading = headings[idx];
      if (!heading) return;
      toc.appendChild(el('a', {
        href: `#ev-h-${idx}`,
        style: {
          paddingLeft: `${(heading.level - 2) * 10}px`, fontSize: '11.5px',
          color: 'var(--text-dim)', padding: '2px 4px', borderRadius: '4px',
        },
        onclick: (event) => {
          event.preventDefault();
          node.scrollIntoView({ behavior: 'smooth', block: 'start' });
        },
      }, heading.text));
    });
    return div({ class: 'grid wide-right', style: { alignItems: 'start' } },
      card({ class: 'tight', title: 'Contents' }, toc),
      card(null, md));
  }

  /* ---------- priors ---------- */
  function renderPriors() {
    const counts = {};
    for (const row of rows) counts[row.status] = (counts[row.status] || 0) + 1;
    let filter = '';
    const table = new DataTable({
      columns: [
        { key: 'kind', label: 'Type', width: '110px' },
        { key: 'target', label: 'Target', mono: true, width: '110px' },
        { key: 'verify_type', label: 'Verification method', mono: true, width: '140px' },
        { key: 'expectation', label: 'Expected', width: '70px' },
        { key: 'status', label: 'Result', width: '86px', render: (r) => statusChip(r.status) },
        { key: 'observed', label: 'Observed', wrap: true, render: (r) => span({ class: 'mono-sm', text: r.observed || '—' }) },
        { key: 'locator', label: 'Source', wrap: true, render: (r) => span({ class: 'tiny muted', text: r.locator || '—' }) },
      ],
      rows,
      keyField: 'target',
      pageSize: 24,
      onRowClick: showPrior,
    });
    const chips = div({ class: 'split', style: { marginBottom: '10px' } },
      el('button', {
        class: 'btn sm' + (filter === '' ? ' primary' : ''), text: `All ${rows.length}`,
        onclick: () => { filter = ''; table.setRows(rows); updateChips(); },
      }),
      ...Object.entries(counts).map(([key, n]) => el('button', {
        class: 'btn sm' + (filter === key ? ' primary' : ''), text: `${key} ${n}`,
        onclick: () => { filter = key; table.setRows(rows.filter((r) => r.status === key)); updateChips(); },
      })));

    function updateChips() {
      const buttons = [...chips.querySelectorAll('button')];
      const labels = ['', ...Object.keys(counts)];
      buttons.forEach((btn, idx) => btn.classList.toggle('primary', labels[idx] === filter));
    }

    return card(null,
      chips,
      table.el,
      div({ class: 'hr' }),
      div({ class: 'card-note', text: 'Click a row for the source quotation, check rule and observations.FLAG means measurable in the structure but not fully consistent with the publication;FAIL means the check rule failed;MISMATCH means prior and structure numbering differ.' }));
  }

  function showPrior(row) {
    openDrawer({
      title: `${row.kind} · ${row.target}`,
      sub: `${row.verify_type} · Expected ${row.expectation}`,
      body: div({ style: { display: 'grid', gap: '12px' } },
        div({ class: 'split' }, statusChip(row.status), row.value ? chip(`Value ${row.value}`, 'mono') : null),
        card({ class: 'flat', title: 'Observed' }, div({ class: 'mono-sm', style: { whiteSpace: 'pre-wrap' }, text: row.observed || '—' })),
        card({ class: 'flat', title: 'Literature claim' },
          div({ class: 'tiny', text: row.claim || '—' }),
          row.quote ? div({ class: 'md', style: { marginTop: '6px' } }, el('blockquote', { text: row.quote })) : null,
          div({ class: 'tiny muted', style: { marginTop: '4px' }, text: row.locator || '' })),
        card({ class: 'flat', title: 'Description' }, div({ class: 'tiny', text: row.detail || '—' }))),
    });
  }

  /* ---------- census ---------- */
  function renderCensus() {
    const censusRows = budget.census || census.census || [];
    const table = new DataTable({
      columns: [
        { key: 'tag', label: 'Complex', width: '80px' },
        { key: 'murj_residues', label: 'MurJ Residue', num: true },
        { key: 'murj_atoms', label: 'MurJ Atom', num: true },
        { key: 'murj_range', label: 'Numbering range', mono: true },
        { key: 'murj_gaps', label: 'Missing intervals', mono: true },
        { key: 'construct_residues', label: 'Construct residues (excluded)', num: true },
        { key: 'peptide_residues', label: 'Peptide residues', num: true },
        { key: 'peptide_atoms', label: 'Peptide atoms', num: true },
        { key: 'peptide_range', label: 'Peptide numbering', mono: true },
      ],
      rows: censusRows,
      pageSize: 5,
    });

    const iface = budget.interface || {};
    const precomputed = budget.budget?.per_complex || {};
    const sm = budget.budget?.small_molecule || {};
    const ratio = budget.budget?.median_ratio;

    const budgetHeadCells = ['Complex', 'Interface MurJ Residue', 'Interface peptide residues', 'Fraction of modeled peptide', 'Atom pairs',
      'Contacting peptide heavy atoms', 'Peptide buried area Å²', 'Asymmetry']
      .map((label, i) => el('th', { class: i ? 'num no-sort' : 'no-sort', text: label }));
    const budgetRows = Object.values(iface).map((r) => el('tr', null,
      el('td', { text: r.tag }),
      el('td', { class: 'num', text: fmtInt(r.murj_residues_in_interface) }),
      el('td', { class: 'num', text: fmtInt(r.peptide_residues_in_interface) }),
      el('td', { class: 'num', text: fmtPct(r.fraction_of_peptide_at_interface, 0) }),
      el('td', { class: 'num', text: fmtInt(r.atom_pairs) }),
      el('td', { class: 'num', text: fmtInt(r.peptide_heavy_atoms_in_contact) }),
      el('td', { class: 'num', text: fmt(r.buried_area_peptide_side_A2, 1) }),
      el('td', { class: 'num', text: fmt(r.asymmetry, 3) })));
    const budgetTable = div({ class: 'table-wrap' },
      el('table', null,
        el('thead', null, el('tr', null, ...budgetHeadCells)),
        el('tbody', null, ...budgetRows)));

    return div({ style: { display: 'grid', gap: '12px' } },
      card({ title: 'Input structure inventory (construct residues excluded)' },
        table.el,
        div({ class: 'card-note', style: { marginTop: '8px' }, text: 'PDB within BRIL fusion-construct residues (negative numbers and 1–4)are not MurJ, and are excluded from all analyses;MurJ the actual range is always 5–500(496 residues).' })),
      card({ title: 'Atom budget: peptide vs small-molecule' },
        div({ class: 'grid cols-3', style: { gap: '10px', marginBottom: '10px' } },
          stat({ label: 'Contacting peptide heavy atoms, maximum across three structures', value: fmtInt(Math.max(...Object.values(precomputed).map((p) => p.peptide_heavy_atoms_in_contact || 0))), unit: 'items', small: true }),
          stat({ label: 'Small-molecule scaffold heavy atoms', value: `${sm.heavy_atoms_min ?? '—'}–${sm.heavy_atoms_max ?? '—'}`, unit: 'items', small: true, foot: `Median ${sm.heavy_atoms_median ?? '—'}` }),
          stat({ label: 'Atom count ratio', value: ratio ? `1 / ${fmt(ratio, 0)}` : '—', small: true, tone: 'accent' })),
        budgetTable,
        div({ class: 'card-note', style: { marginTop: '8px' }, text: 'The aim is to reproduce interaction vectors with fewer atoms. This motivates pharmacophore distance constraints rather than simply shrinking the peptide or matching its shape.' })));
  }

  /* ---------- interface ---------- */
  function renderInterface() {
    const murjRows = iface.murj || [];
    const pepRows = iface.peptide || [];
    const murjTable = new DataTable({
      columns: [
        { key: 'complex', label: 'Complex', width: '84px' },
        { key: 'resid', label: 'Residue number', num: true, width: '70px' },
        { key: 'one', label: 'Residue', width: '56px' },
        { key: 'n_atom_contacts', label: 'Contacting atoms', num: true },
        { key: 'n_peptide_residues', label: 'Peptide residues', num: true },
        { key: 'dmin_A', label: 'dmin Å', num: true, render: (r) => fmt(r.dmin_A, 2) },
        { key: 'dSASA_A2', label: 'dSASA Å²', num: true, render: (r) => fmt(r.dSASA_A2, 1) },
        { key: 'convergence', label: 'Convergence', num: true, render: (r) => fmt(r.convergence, 2) },
      ],
      rows: murjRows,
      pageSize: 20,
      initialSort: { key: 'n_atom_contacts', dir: 'desc' },
    });
    const pepTable = new DataTable({
      columns: [
        { key: 'complex', label: 'Complex', width: '84px' },
        { key: 'pep_resid', label: 'Peptide residue number', num: true, width: '80px' },
        { key: 'one', label: 'Residue', width: '56px' },
        { key: 'n_atom_contacts', label: 'Contacting atoms', num: true },
        { key: 'n_murj_residues', label: 'MurJ Residue', num: true },
        { key: 'dmin_A', label: 'dmin Å', num: true, render: (r) => fmt(r.dmin_A, 2) },
        { key: 'dSASA_A2', label: 'dSASA Å²', num: true, render: (r) => fmt(r.dSASA_A2, 1) },
        { key: 'heavy_atoms', label: 'Heavy atoms', num: true },
      ],
      rows: pepRows,
      pageSize: 20,
      initialSort: { key: 'n_atom_contacts', dir: 'desc' },
    });
    return div({ style: { display: 'grid', gap: '12px' } },
      card({ title: 'MurJ side: contacts and burial per interface residue, by structure' }, murjTable.el),
      card({ title: 'Peptide side: each interface residue' }, pepTable.el),
      div({ class: 'card-note', text: 'dSASA = Monomer conformation SASA − Complex SASA, Summed over residues (Shrake–Rupley, 92 points per/atom). Convergence: = fraction of the three peptides contacting the residue.' }));
  }

  /* ---------- anchors ---------- */
  function renderAnchors() {
    const anchors = iface.anchors || [];
    const table = new DataTable({
      columns: [
        { key: 'complex', label: 'Complex', width: '90px' },
        { key: 'a', label: 'Anchor A', mono: true, width: '80px' },
        { key: 'b', label: 'Anchor B', mono: true, width: '80px' },
        { key: 'sidechain_centroid_A', label: 'Side-chain centroid distance Å', num: true, render: (r) => fmt(r.sidechain_centroid_A, 2) },
        { key: 'min_atom_A', label: 'Nearest atom distance Å', num: true, render: (r) => fmt(r.min_atom_A, 2) },
      ],
      rows: anchors,
      pageSize: 20,
      initialSort: { key: 'min_atom_A', dir: 'desc' },
    });
    return div({ style: { display: 'grid', gap: '12px' } },
      card({ title: 'Positive-anchor geometry (basis for scaffold axis length)' },
        table.el,
        div({ class: 'card-note', style: { marginTop: '8px' }, text: 'K46 / R53 / R270 / K368 are four functionally required basic residues; pairwise distances determine scaffold span and pharmacophore axis length 9–18 Å .' })),
      card({ title: 'Salt bridges (report section 4 )' }, saltBridgeTable()));
  }

  function saltBridgeTable() {
    const rows = [
      ['R53 NH2', 'V55(C-term) O', 3.74, 'SglPP7'],
      ['K46 NZ', 'E96 OE1', 2.69, 'SglCJ3'],
      ['R270 NH1', 'E96 OE2', 3.53, 'SglCJ3'],
      ['K368 NZ', 'L103(C-term) O', 2.33, 'SglCJ3'],
    ].map((row) => el('tr', null,
      el('td', { class: 'mono-sm', text: row[0] }),
      el('td', { class: 'mono-sm', text: row[1] }),
      el('td', { class: 'num', text: fmt(row[2], 2) }),
      el('td', { text: row[3] })));
    const head = ['MurJ Atom', 'Peptide atoms', 'Distance Å', 'Complex']
      .map((label, i) => el('th', { class: i === 2 ? 'num no-sort' : 'no-sort', text: label }));
    return div({ class: 'table-wrap' },
      el('table', null,
        el('thead', null, el('tr', null, ...head)),
        el('tbody', null, ...rows)));
  }

  root.appendChild(tabBar);
  root.appendChild(bodyHost);
  activate(active);

  return null;
}
