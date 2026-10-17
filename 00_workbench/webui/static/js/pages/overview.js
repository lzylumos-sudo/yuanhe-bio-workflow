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
  root.appendChild(loadingPane('正在读取管线结果…'));
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
      el('h1', { text: '设计管线概览' }),
      div({ class: 'sub', text: '从三个冷冻电镜结构（9NU4 / 9NU5 / 9NU8）出发，把三种噬菌体裂解肽趋同的抑制逻辑，翻译成可枚举、可解释、可合成的小分子骨架。所有数字都由本机从坐标文件现算。' })),
    div({ class: 'actions' },
      el('button', { class: 'btn', text: '3D 观察', onclick: () => ctx.navigate('structure', { ligand: '0' }) }),
      el('button', { class: 'btn primary', text: '运行管线', onclick: () => ctx.navigate('run') }))));

  /* ---- status strip ---- */
  const byStatus = ev.by_status || {};
  const statusChips = Object.entries(byStatus)
    .sort((a, b) => b[1] - a[1])
    .map(([key, n]) => span({ class: 'split' }, statusChip(key), span({ class: 'tiny muted', text: `×${n}` })));

  root.appendChild(div({ class: 'grid cols-auto', style: { marginBottom: '12px' } },
    card({ class: 'tight' },
      div({ class: 'split' },
        span({ class: 'tiny muted', text: '证据链（先验核对）' }),
        span({ class: 'right tiny muted', style: { marginLeft: 'auto' }, text: timeAgo(ev.generated || 0) })),
      div({ class: 'split', style: { marginTop: '6px' } }, ...statusChips)),
    card({ class: 'tight' },
      div({ class: 'tiny muted', text: '结构过滤（F1–F4）' }),
      div({ class: 'split', style: { marginTop: '6px' } },
        chip(`${filters.n_input ?? 0} 输入`, 'info'),
        span({ class: 'tiny muted', text: '→' }),
        chip(`${nKeep} 存活`, 'pass'),
        span({ class: 'tiny muted', text: `最大跨结构分差 ${filters.max_spread ?? '—'} Å·score` }))),
    card({ class: 'tight' },
      div({ class: 'tiny muted', text: '结构输入' }),
      div({ class: 'split', style: { marginTop: '6px' } },
        ...['9NU4', '9NU5', '9NU8'].map((id) => chip(id, 'mono info')),
        span({ class: 'tiny muted', text: '3.6–3.7 Å' }))),
    card({ class: 'tight' },
      div({ class: 'tiny muted', text: '结果文件生成时间' }),
      div({ class: 'mono-sm', style: { marginTop: '6px' }, text: fmtDate(Math.max(...Object.values(meta.files || { x: 0 }))) }))));

  /* ---- KPI row ---- */
  root.appendChild(div({ class: 'grid cols-4' },
    card(null, stat({
      label: '枚举骨架', value: fmtInt(nEnumerated), unit: '个',
      foot: `${axes.length} 条设计轴 · 性质包络内 ${fmtInt(fun.property_pass)} 个进入位姿打分`,
    })),
    card(null, stat({
      label: '三重达标', value: fmtInt(nDrug), unit: '个', tone: 'accent',
      foot: `无位阻 ${fmtInt(fun.no_clash)} · 肽模拟 ${fmtInt(fun.peptide_mimic)}`,
    })),
    card(null, stat({
      label: '结构过滤存活', value: fmtInt(nKeep), unit: '个', tone: 'good',
      foot: '锚点 + 跨亚口袋 + 三结构一致性（F4 是主要瓶颈）',
    })),
    card(null, stat({
      label: '最佳骨架总分', value: bestTotal === null ? '—' : fmt(bestTotal, 2), tone: 'warn',
      foot: best.length ? `${best[0].props?.MW} Da · cLogP ${best[0].props?.cLogP} · TPSA ${best[0].props?.TPSA}` : '',
    }))));

  /* ---- funnel + best ---- */
  const funnelCard = card({
    title: '筛选漏斗',
    right: span({ class: 'tiny muted', text: '同一轴上的确定性枚举' }),
  },
    funnel([
      { label: '骨架枚举', value: nEnumerated, hint: '长度约束精确求解' },
      { label: '性质包络通过', value: fun.property_pass ?? 0, hint: 'MW ≤ 500, cLogP/TPSA 包络' },
      { label: '完成 3D 位姿', value: fun.posed ?? 0, hint: '刚体放置 + 12 个旋转' },
      { label: '无位阻冲突', value: fun.no_clash ?? 0, hint: 'clash < 0.25' },
      { label: '肽模拟达标', value: fun.peptide_mimic ?? 0, hint: '与三肽空间重合度 ≥ 0.5' },
      { label: '三重达标', value: nDrug, hint: '位阻 + 类药 + 肽模拟', final: true },
    ]),
    div({ class: 'hr' }),
    div({ class: 'split tiny muted' },
      span({ text: `结构衍生过滤：${filters.n_input ?? 0} → ` }),
      span({ style: { color: 'var(--accent-2)' }, text: `${nKeep}` }),
      span({ text: `（锚点 ${filters.by_stage?.anchor ?? '—'} · 桥接 ${filters.by_stage?.bridge ?? '—'} · 三结构一致 ${filters.by_stage?.ensemble ?? '—'}）` })));

  const bestCard = card({
    title: '排名前十骨架',
    right: el('button', { class: 'btn sm ghost', text: '查看全部 →', onclick: () => ctx.navigate('scaffolds') }),
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
    title: '界面热点（前 10）',
    right: el('button', { class: 'btn sm ghost', text: '全部 52 个 →', onclick: () => ctx.navigate('hotspots') }),
  },
    barList(hot.slice(0, 10).map((h) => ({
      label: `${h.resname}${h.murj_resid}`,
      value: h.score,
      color: h.evidence ? '#ff8a5b' : '#3f7fb5',
      title: h.evidence || `${h.resname}${h.murj_resid}`,
    })), { format: (v) => fmt(v, 3) }),
    div({ class: 'hr' }),
    legend([
      { color: '#ff8a5b', label: '有遗传学/耐药证据' },
      { color: '#3f7fb5', label: '仅几何证据' },
    ]));

  const pharmCard = card({
    title: '药效团共识',
    right: el('button', { class: 'btn sm ghost', text: '详情 →', onclick: () => ctx.navigate('pharmacophore') }),
  },
    div({ class: 'grid cols-2', style: { gap: '8px' } },
      stat({ label: '共识特征', value: feats.length, unit: '个', small: true }),
      stat({ label: '设计轴', value: axes.length, unit: '条', small: true })),
    div({ class: 'hr' }),
    div({ class: 'chips' },
      ...Object.entries(typeCount).map(([type, n]) => span({
        class: 'chip', style: { borderColor: 'rgba(255,255,255,.12)' },
      }, span({ class: 'swatch', style: { background: FEATURE_COLOR[type] || '#8b97a8' } }), `${type} ×${n}`))),
    div({ class: 'card-note', style: { marginTop: '8px' }, text: '每种特征由三种噬菌体肽的功能原子在三结构共帧下聚类得到，weight 已按趋同度加权。' }));

  root.appendChild(div({ class: 'section grid cols-2' }, hotCard, pharmCard));

  /* ---- pocket + figures ---- */
  const subs = pocket.subpockets || [];
  const pocketCard = card({
    title: '埋藏腔体与亚口袋',
    right: el('button', { class: 'btn sm ghost', text: '口袋详情 →', onclick: () => ctx.navigate('pocket') }),
  },
    div({ class: 'split', style: { marginBottom: '8px' } },
      stat({ label: '腔体体积', value: fmtInt(pocket.total_volume_A3 || 0), unit: 'Å³', small: true }),
      stat({ label: '亚口袋', value: subs.length, unit: '个', small: true })),
    div({ class: 'bar-list' },
      ...subs.map((sp) => div({ class: 'row' },
        div({ class: 'k', text: sp.subpocket.replace(/^S\d_/, ''), title: sp.subpocket }),
        div({ class: 'track' },
          div({ class: 'fill', style: { width: `${(sp.volume_A3 / Math.max(...subs.map((x) => x.volume_A3), 1)) * 100}%`, background: SUBPOCKET_COLOR[sp.subpocket] || '#8b97a8' } })),
        div({ class: 'num', text: fmtInt(sp.volume_A3) })))),
    div({ class: 'card-note', style: { marginTop: '8px' }, text: '体积 = 体素 × 0.8³ Å³；S2 羧酸闸门埋藏度 0.94，是最深的位点。' }));

  const fileCard = card({ title: '产物清单' },
    div({ class: 'kv tiny' },
      ...Object.entries(meta.files || {}).map(([name, mtime]) => [
        el('div', { class: 'k mono-sm', style: { wordBreak: 'break-all' }, text: name }),
        div({ class: 'v tiny muted', text: mtime ? timeAgo(mtime) : '缺失' }),
      ]).flat()));

  root.appendChild(div({ class: 'section grid cols-2' }, pocketCard, fileCard));

  const figs = (meta.figures || []).filter((f) => f.name.endsWith('.png'));
  if (figs.length) {
    root.appendChild(card({ title: '管线示意图', class: 'section' },
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
