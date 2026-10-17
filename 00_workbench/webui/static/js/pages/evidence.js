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
      el('h1', { text: '证据链' }),
      div({ class: 'sub', text: '本页所有数字都由 run_evidence.py 从三个 PDB 现算，并按 priors/murj_priors.json 里的引用与规则逐条核对文献断言。复现不了的条目会以 FLAG / FAIL 显式列出，而不是被静默接受。' })),
    div({ class: 'actions' },
      el('button', {
        class: 'btn', text: '重跑证据链',
        onclick: () => ctx.navigate('run'),
      }),
      el('a', { class: 'btn', href: '/files/evidence/report.md?download=1', text: '下载报告' }))));

  root.appendChild(div({ class: 'grid cols-4', style: { marginBottom: '12px' } },
    card({ class: 'tight' }, stat({ label: '核对条目', value: status.n ?? rows.length, unit: '条', small: true })),
    card({ class: 'tight' }, stat({
      label: 'PASS', value: byStatus.PASS ?? 0, unit: '条', small: true, tone: 'good',
      foot: `通过率 ${budget.summary ? fmtPct(budget.summary.pass_rate) : '—'}`,
    })),
    card({ class: 'tight' }, stat({ label: 'FLAG', value: byStatus.FLAG ?? 0, unit: '条', small: true, tone: 'warn', foot: '需要人工判断' })),
    card({ class: 'tight' }, stat({ label: 'FAIL / MISMATCH', value: (byStatus.FAIL ?? 0) + (byStatus.MISMATCH ?? 0), unit: '条', small: true, tone: 'bad' }))));

  /* ---------- tabs ---------- */
  const TABS = [
    { id: 'report', label: '证据链报告' },
    { id: 'priors', label: '先验逐条核对', n: rows.length },
    { id: 'census', label: '结构普查与原子预算' },
    { id: 'interface', label: '界面测量' },
    { id: 'anchors', label: '锚点几何' },
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
    const md = markdownNode(report.markdown || '_（缺少 evidence/report.md，请先运行 run_evidence.py）_');
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
      card({ class: 'tight', title: '目录' }, toc),
      card(null, md));
  }

  /* ---------- priors ---------- */
  function renderPriors() {
    const counts = {};
    for (const row of rows) counts[row.status] = (counts[row.status] || 0) + 1;
    let filter = '';
    const table = new DataTable({
      columns: [
        { key: 'kind', label: '类型', width: '110px' },
        { key: 'target', label: '目标', mono: true, width: '110px' },
        { key: 'verify_type', label: '核对方式', mono: true, width: '140px' },
        { key: 'expectation', label: '期望', width: '70px' },
        { key: 'status', label: '结果', width: '86px', render: (r) => statusChip(r.status) },
        { key: 'observed', label: '观测', wrap: true, render: (r) => span({ class: 'mono-sm', text: r.observed || '—' }) },
        { key: 'locator', label: '出处', wrap: true, render: (r) => span({ class: 'tiny muted', text: r.locator || '—' }) },
      ],
      rows,
      keyField: 'target',
      pageSize: 24,
      onRowClick: showPrior,
    });
    const chips = div({ class: 'split', style: { marginBottom: '10px' } },
      el('button', {
        class: 'btn sm' + (filter === '' ? ' primary' : ''), text: `全部 ${rows.length}`,
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
      div({ class: 'card-note', text: '点击任意一行查看原文引用、核规则与观测细节。FLAG 表示结构里能测到但与论文说法不完全一致；FAIL 表示核规则不成立；MISMATCH 表示先验文件与结构编号体系不一致。' }));
  }

  function showPrior(row) {
    openDrawer({
      title: `${row.kind} · ${row.target}`,
      sub: `${row.verify_type} · 期望 ${row.expectation}`,
      body: div({ style: { display: 'grid', gap: '12px' } },
        div({ class: 'split' }, statusChip(row.status), row.value ? chip(`值 ${row.value}`, 'mono') : null),
        card({ class: 'flat', title: '观测' }, div({ class: 'mono-sm', style: { whiteSpace: 'pre-wrap' }, text: row.observed || '—' })),
        card({ class: 'flat', title: '文献断言' },
          div({ class: 'tiny', text: row.claim || '—' }),
          row.quote ? div({ class: 'md', style: { marginTop: '6px' } }, el('blockquote', { text: row.quote })) : null,
          div({ class: 'tiny muted', style: { marginTop: '4px' }, text: row.locator || '' })),
        card({ class: 'flat', title: '说明' }, div({ class: 'tiny', text: row.detail || '—' }))),
    });
  }

  /* ---------- census ---------- */
  function renderCensus() {
    const censusRows = budget.census || census.census || [];
    const table = new DataTable({
      columns: [
        { key: 'tag', label: '复合物', width: '80px' },
        { key: 'murj_residues', label: 'MurJ 残基', num: true },
        { key: 'murj_atoms', label: 'MurJ 原子', num: true },
        { key: 'murj_range', label: '编号范围', mono: true },
        { key: 'murj_gaps', label: '缺失区间', mono: true },
        { key: 'construct_residues', label: '构建体残基（已剔除）', num: true },
        { key: 'peptide_residues', label: '肽残基', num: true },
        { key: 'peptide_atoms', label: '肽原子', num: true },
        { key: 'peptide_range', label: '肽编号', mono: true },
      ],
      rows: censusRows,
      pageSize: 5,
    });

    const iface = budget.interface || {};
    const precomputed = budget.budget?.per_complex || {};
    const sm = budget.budget?.small_molecule || {};
    const ratio = budget.budget?.median_ratio;

    const budgetHeadCells = ['复合物', '界面 MurJ 残基', '界面肽残基', '占建模肽', '原子对',
      '肽侧接触重原子', '肽侧埋藏 Å²', '不对称度']
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
      card({ title: '输入结构里到底有什么（构建体残基已剔除）' },
        table.el,
        div({ class: 'card-note', style: { marginTop: '8px' }, text: 'PDB 中的 BRIL 融合构建体残基（负编号与 1–4）不是 MurJ，已在所有分析中剔除；MurJ 实际范围恒为 5–500（496 残基）。' })),
      card({ title: '原子预算：肽 vs 小分子' },
        div({ class: 'grid cols-3', style: { gap: '10px', marginBottom: '10px' } },
          stat({ label: '肽接触重原子（三结构最大）', value: fmtInt(Math.max(...Object.values(precomputed).map((p) => p.peptide_heavy_atoms_in_contact || 0))), unit: '个', small: true }),
          stat({ label: '小分子骨架重原子', value: `${sm.heavy_atoms_min ?? '—'}–${sm.heavy_atoms_max ?? '—'}`, unit: '个', small: true, foot: `中位数 ${sm.heavy_atoms_median ?? '—'}` }),
          stat({ label: '原子数比值', value: ratio ? `1 / ${fmt(ratio, 0)}` : '—', small: true, tone: 'accent' })),
        budgetTable,
        div({ class: 'card-note', style: { marginTop: '8px' }, text: '结论不是“把肽缩小”，而是“用更少的原子复刻同样的相互作用向量”——这正是把设计目标写成药效团点间距离约束、而不是形状匹配的定量理由。' })));
  }

  /* ---------- interface ---------- */
  function renderInterface() {
    const murjRows = iface.murj || [];
    const pepRows = iface.peptide || [];
    const murjTable = new DataTable({
      columns: [
        { key: 'complex', label: '复合物', width: '84px' },
        { key: 'resid', label: '残基号', num: true, width: '70px' },
        { key: 'one', label: '残基', width: '56px' },
        { key: 'n_atom_contacts', label: '接触原子', num: true },
        { key: 'n_peptide_residues', label: '肽残基', num: true },
        { key: 'dmin_A', label: 'dmin Å', num: true, render: (r) => fmt(r.dmin_A, 2) },
        { key: 'dSASA_A2', label: 'dSASA Å²', num: true, render: (r) => fmt(r.dSASA_A2, 1) },
        { key: 'convergence', label: '收敛度', num: true, render: (r) => fmt(r.convergence, 2) },
      ],
      rows: murjRows,
      pageSize: 20,
      initialSort: { key: 'n_atom_contacts', dir: 'desc' },
    });
    const pepTable = new DataTable({
      columns: [
        { key: 'complex', label: '复合物', width: '84px' },
        { key: 'pep_resid', label: '肽残基号', num: true, width: '80px' },
        { key: 'one', label: '残基', width: '56px' },
        { key: 'n_atom_contacts', label: '接触原子', num: true },
        { key: 'n_murj_residues', label: 'MurJ 残基', num: true },
        { key: 'dmin_A', label: 'dmin Å', num: true, render: (r) => fmt(r.dmin_A, 2) },
        { key: 'dSASA_A2', label: 'dSASA Å²', num: true, render: (r) => fmt(r.dSASA_A2, 1) },
        { key: 'heavy_atoms', label: '重原子', num: true },
      ],
      rows: pepRows,
      pageSize: 20,
      initialSort: { key: 'n_atom_contacts', dir: 'desc' },
    });
    return div({ style: { display: 'grid', gap: '12px' } },
      card({ title: 'MurJ 侧：每个界面残基的接触与埋藏（三结构逐行）' }, murjTable.el),
      card({ title: '肽侧：每个进入界面的肽残基' }, pepTable.el),
      div({ class: 'card-note', text: 'dSASA = 单体构象 SASA − 复合物 SASA，按残基求和（Shrake–Rupley，92 点/原子）。收敛度 = 三种肽中接触该残基的比例。' }));
  }

  /* ---------- anchors ---------- */
  function renderAnchors() {
    const anchors = iface.anchors || [];
    const table = new DataTable({
      columns: [
        { key: 'complex', label: '复合物', width: '90px' },
        { key: 'a', label: '锚点 A', mono: true, width: '80px' },
        { key: 'b', label: '锚点 B', mono: true, width: '80px' },
        { key: 'sidechain_centroid_A', label: '侧链质心距离 Å', num: true, render: (r) => fmt(r.sidechain_centroid_A, 2) },
        { key: 'min_atom_A', label: '最近原子距离 Å', num: true, render: (r) => fmt(r.min_atom_A, 2) },
      ],
      rows: anchors,
      pageSize: 20,
      initialSort: { key: 'min_atom_A', dir: 'desc' },
    });
    return div({ style: { display: 'grid', gap: '12px' } },
      card({ title: '正电锚点之间的几何（骨架轴长的来源）' },
        table.el,
        div({ class: 'card-note', style: { marginTop: '8px' }, text: 'K46 / R53 / R270 / K368 是四个功能必需的碱性残基；它们两两之间的距离决定了小分子需要跨多长，这也是药效团轴长 9–18 Å 的来源。' })),
      card({ title: '盐桥（来自报告第 4 节）' }, saltBridgeTable()));
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
    const head = ['MurJ 原子', '肽原子', '距离 Å', '复合物']
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
