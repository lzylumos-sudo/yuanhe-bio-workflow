/* pages/pocket.js -- the buried cavity: volume, chemistry, sub-pocket anatomy. */

import { api } from '../api.js';
import {
  el, div, span, card, stat, chip, DataTable, kvList, meterRow,
  fmt, fmtInt, loadingPane, legend,
} from '../ui.js';

const CHEM = [
  { key: 'pos', label: '碱性 (K/R/H)', color: '#4f8ef7' },
  { key: 'neg', label: '酸性 (D/E)', color: '#e55353' },
  { key: 'polar', label: '极性 (S/T/N/Q)', color: '#38a169' },
  { key: 'aromatic', label: '芳香 (Y/W/F)', color: '#c99a3a' },
  { key: 'hydrophobic', label: '疏水', color: '#8b97a8' },
];

export async function mount(root, ctx) {
  root.appendChild(loadingPane('正在读取口袋体素划分…'));
  const data = await api.pocket();
  const pocket = data.pocket || {};
  const subs = pocket.subpockets || [];
  const total = Number(pocket.total_volume_A3) || 0;
  const deepest = subs.reduce((acc, s) => Math.max(acc, Number(s.mean_buriedness) || 0), 0);

  root.replaceChildren();
  root.appendChild(div({ class: 'page-head' },
    div(null,
      el('h1', { text: '埋藏腔体与亚口袋' }),
      div({ class: 'sub', text: '肽给出的是腔体的“负像”。0.8 Å 体素栅格化后，用 50 条射线计算每个体素的埋藏度（SiteMap 的轻量近似），再按最近的亚口袋质心划分化学区域。' })),
    div({ class: 'actions' },
      el('a', { class: 'btn', href: '/files/results/pocket.json?download=1', text: '下载 JSON' }),
      el('button', {
        class: 'btn primary', text: '在 3D 中查看腔体',
        onclick: () => ctx.navigate('structure', { ligand: '0', pocket: '1' }),
      }))));

  root.appendChild(div({ class: 'grid cols-4', style: { marginBottom: '12px' } },
    card({ class: 'tight' }, stat({ label: '腔体体积', value: fmtInt(total), unit: 'Å³', small: true, tone: 'accent' })),
    card({ class: 'tight' }, stat({ label: '亚口袋', value: subs.length, unit: '个', small: true })),
    card({ class: 'tight' }, stat({ label: '最深埋藏度', value: fmt(deepest, 2), small: true, tone: 'good', foot: 'S2 羧酸闸门' })),
    card({ class: 'tight' }, stat({
      label: '平均亚口袋体积', value: subs.length ? fmt(subs.reduce((a, s) => a + Number(s.volume_A3), 0) / subs.length, 0) : '—', unit: 'Å³', small: true,
    }))));

  /* ---- sub-pocket anatomy ---- */
  const cards = div({ class: 'grid cols-2', style: { gap: '12px' } },
    ...subs.map((sp) => {
      const chem = sp.chemistry || {};
      const chemTotal = CHEM.reduce((acc, c) => acc + (chem[c.key] || 0), 0) || 1;
      return card({ class: 'tight' },
        div({ class: 'split' },
          span({ class: 'chip mono', text: sp.subpocket }),
          span({ class: 'right tiny muted', style: { marginLeft: 'auto' }, text: `${fmtInt(sp.volume_A3)} Å³ · ${fmtInt(sp.n_voxels)} 体素` })),
        div({ style: { marginTop: '8px' } },
          div({ class: 'tiny muted', style: { marginBottom: '4px' }, text: '化学组成（按最近残基分类）' }),
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
          meterRow('平均埋藏度', sp.mean_buriedness, 1, 'var(--accent-2)'),
          div({ class: 'kv tiny' },
            el('div', { class: 'k', text: '质心' }),
            el('div', { class: 'v mono-sm', text: `(${(sp.centroid || []).map((v) => fmt(v, 1)).join(', ')})` }),
            el('div', { class: 'k', text: '范围 Å' }),
            el('div', { class: 'v mono-sm', text: (sp.extent_A || []).map((v) => fmt(v, 1)).join(' × ') }))),
        div({ class: 'chips', style: { marginTop: '10px' } },
          ...(sp.residues || []).map((r) => span({ class: 'resid', text: r }))));
    }));

  root.appendChild(card({ title: '四个亚口袋', class: 'section' }, cards));

  /* ---- table + method ---- */
  const table = new DataTable({
    columns: [
      { key: 'subpocket', label: '亚口袋', mono: true },
      { key: 'volume_A3', label: '体积 Å³', num: true, render: (r) => fmt(r.volume_A3, 1) },
      { key: 'n_voxels', label: '体素', num: true },
      { key: 'mean_buriedness', label: '平均埋藏度', num: true, render: (r) => fmt(r.mean_buriedness, 3) },
      { key: 'chemistry', label: '主导化学', render: (r) => dominant(r.chemistry) },
      { key: 'residues', label: '衬里残基', wrap: true, render: (r) => (r.residues || []).join(' ') },
    ],
    rows: subs,
    pageSize: 10,
    initialSort: { key: 'volume_A3', dir: 'desc' },
  });

  const method = card({ title: '这个数字怎么来的' },
    div({ class: 'tiny', style: { display: 'grid', gap: '8px' } },
      div(null, span({ class: 'muted', text: '1 · 体素化：' }), '在肽包络（±7 Å）内以 0.8 Å 栅格撒点，剔除与 MurJ 重原子距离 < 3.0 Å 的体素。'),
      div(null, span({ class: 'muted', text: '2 · 埋藏度：' }), '每个体素发 50 条射线（球面 Fabonacci 分布），统计 13 Å 内被受体原子（probe 1.8 Å）挡住的射线比例。'),
      div(null, span({ class: 'muted', text: '3 · 保留：' }), '只保留埋藏度 ≥ 0.5 的体素，其体积和就是 5,440 Å³。'),
      div(null, span({ class: 'muted', text: '4 · 化学标注：' }), '每个体素取最近的 MurJ 残基，按残基类型分类；亚口袋归属取最近的亚口袋残基质心。'),
      div({ class: 'card-note', text: '这是几何量而非实验量：它能说明“腔体足够大、足够深、化学是混合的”，但不能单独证明可成药。' })));

  root.appendChild(div({ class: 'section grid wide-left' },
    card({ title: '亚口袋数据' }, table.el),
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
