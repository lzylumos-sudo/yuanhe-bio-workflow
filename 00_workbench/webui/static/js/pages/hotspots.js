/* pages/hotspots.js -- MurJ interface hotspots and the per-structure contacts. */

import { api } from '../api.js';
import {
  el, div, span, card, stat, chip, barList, DataTable, kvList, legend,
  fmt, fmtInt, loadingPane, statusChip,
} from '../ui.js';

export async function mount(root, ctx) {
  root.appendChild(loadingPane('正在汇总界面残基…'));
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
      el('h1', { text: '界面热点' }),
      div({ class: 'sub', text: '三种独立进化的噬菌体肽在 MurJ 上的接触面。score = 0.30·跨肽收敛度 + 0.25·埋藏 + 0.15·接触密度 + 0.20·功能必需 + 0.10·耐药证据——全部可复算，没有任何数字是手写的。' })),
    div({ class: 'actions' },
      el('a', { class: 'btn', href: '/files/results/interface_hotspots.csv?download=1', text: '下载 CSV' }),
      el('button', { class: 'btn', text: '在 3D 中查看', onclick: () => ctx.navigate('structure', { ligand: '0' }) }))));

  root.appendChild(div({ class: 'grid cols-4', style: { marginBottom: '12px' } },
    card({ class: 'tight' }, stat({ label: '界面残基', value: rows.length, unit: '个', small: true })),
    card({ class: 'tight' }, stat({ label: '三种肽共同接触', value: nCommon, unit: '个', small: true, tone: 'accent' })),
    card({ class: 'tight' }, stat({ label: '有遗传学证据', value: nEvidence, unit: '个', small: true, tone: 'warn' })),
    card({ class: 'tight' }, stat({
      label: '最大 dSASA', unit: 'Å²', small: true,
      value: fmt(Math.max(...rows.map((r) => Number(r.dSASA_A2) || 0)), 0),
    }))));

  /* ---- chart + common residues ---- */
  const chartCard = card({ title: '热点排序（前 14）' },
    barList(rows.slice(0, 14).map((r) => ({
      label: r.label, value: r.score, title: r.evidence || '几何证据',
      color: r.has_evidence ? '#ff8a5b' : '#3f7fb5',
    })), { format: (v) => fmt(v, 3) }),
    div({ class: 'hr' }),
    card({ class: 'flat tight' },
      div({ class: 'tiny muted', style: { marginBottom: '6px' } },
        '三种肽共同接触（几何收敛 = 1.0，n = ', String(nCommon), '）'),
      div({ class: 'chips' },
        ...rows.filter((r) => Number(r.convergence) >= 1).map((r) => span({ class: 'resid', text: `${r.resname}${r.murj_resid}` })))),
    div({ class: 'card-note', style: { marginTop: '8px' }, text: '注意：这里的“收敛”是几何收敛（三结构都接触），不是序列保守性；两者数值不同，混用会出错。' }));

  /* ---- detail panel ---- */
  const detail = div({ class: 'card' },
    div({ class: 'card-head' }, el('h3', { text: '残基详情' })),
    div({ class: 'card-note', text: '点击表格中的任意残基，查看它在三个结构里的接触原子数、最近距离与埋藏面积。' }));

  function showDetail(row) {
    const perComplex = interfaceRows.filter((r) => Number(r.resid) === Number(row.murj_resid));
    detail.replaceChildren(
      div({ class: 'card-head' },
        el('h3', { text: `${row.resname}${row.murj_resid}` }),
        span({ class: 'right chip mono', text: `score ${fmt(row.score, 3)}` })),
      div({ class: 'grid cols-2', style: { gap: '8px', marginBottom: '10px' } },
        stat({ label: '收敛度', value: fmt(row.convergence, 2), small: true, foot: '三结构中出现比例' }),
        stat({ label: '埋藏面积', value: fmt(row.dSASA_A2, 0), unit: 'Å²', small: true })),
      kvList([
        ['复合物', (row.complexes || '').split(';').filter(Boolean).join(' · ') || '—'],
        ['接触原子对数', fmtInt(row.n_contacts)],
        ['接触肽残基数', fmtInt(row.n_peptide_residues)],
        ['最近重原子距离', row.dmin_A ? `${fmt(row.dmin_A, 2)} Å` : '—'],
      ]),
      row.evidence ? div({ class: 'card flat tight', style: { marginTop: '10px' } },
        div({ class: 'tiny muted', text: '实验证据' }),
        div({ class: 'tiny', style: { marginTop: '2px' }, text: row.evidence })) : null,
      div({ class: 'hr' }),
      div({ class: 'tiny muted', style: { marginBottom: '4px' }, text: '逐结构测量（evidence/01_interface_murj.csv）' }),
      div({ class: 'table-wrap', style: { border: 'none' } },
        el('table', null,
          el('thead', null, el('tr', null,
            ...['复合物', '接触原子', '肽残基', 'dmin Å', 'dSASA Å²'].map((t, i) => el('th', { class: i ? 'num no-sort' : 'no-sort', text: t })))),
          el('tbody', null,
            ...(perComplex.length ? perComplex : [{ complex: '—' }]).map((r) => el('tr', null,
              el('td', { text: r.complex || '—' }),
              el('td', { class: 'num', text: fmtInt(r.n_atom_contacts) }),
              el('td', { class: 'num', text: fmtInt(r.n_peptide_residues) }),
              el('td', { class: 'num', text: fmt(r.dmin_A, 2) }),
              el('td', { class: 'num', text: fmt(r.dSASA_A2, 1) })))))),
      div({ class: 'split', style: { marginTop: '10px' } },
        el('button', {
          class: 'btn sm', text: '在口袋中定位 →',
          onclick: () => ctx.navigate('structure', { ligand: '0', focus: String(row.murj_resid) }),
        })));
    table.select(row.murj_resid);
  }

  /* ---- table ---- */
  const table = new DataTable({
    columns: [
      { key: 'murj_resid', label: '残基号', num: true, width: '70px' },
      { key: 'label', label: '残基', width: '60px' },
      { key: 'score', label: 'score', num: true, render: (r) => fmt(r.score, 3) },
      { key: 'convergence', label: '收敛度', num: true, render: (r) => fmt(r.convergence, 2) },
      { key: 'dSASA_A2', label: 'dSASA Å²', num: true, render: (r) => fmt(r.dSASA_A2, 1) },
      { key: 'n_contacts', label: '接触原子', num: true },
      { key: 'n_peptide_residues', label: '肽残基', num: true },
      { key: 'dmin_A', label: 'dmin Å', num: true, render: (r) => fmt(r.dmin_A, 2) },
      { key: 'complexes', label: '复合物', mono: true },
      {
        key: 'evidence', label: '实验证据', sortable: false, wrap: true,
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
  root.appendChild(card({ title: '全部界面残基', class: 'section' }, table.el));
  showDetail(rows[0]);

  return null;
}
