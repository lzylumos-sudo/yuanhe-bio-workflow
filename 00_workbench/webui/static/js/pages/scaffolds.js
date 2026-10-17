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
  root.appendChild(loadingPane('正在读取 200 个排名骨架…'));
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
      { key: 'axis', label: '设计轴', width: '76px' },
      { key: 'MW', label: 'MW', num: true, render: (r) => fmt(r.MW, 1) },
      { key: 'cLogP', label: 'cLogP', num: true, render: (r) => fmt(r.cLogP, 2) },
      { key: 'TPSA', label: 'TPSA', num: true, render: (r) => fmt(r.TPSA, 0) },
      { key: 'heavy_atoms', label: '重原子', num: true },
      { key: 'peptide_overlap', label: '肽模拟', num: true, render: (r) => fmt(r.peptide_overlap, 2) },
      { key: 'clash', label: '位阻', num: true, render: (r) => fmt(r.clash, 2) },
      {
        key: 'total', label: '总分', num: true, render: (r) => r.total,
      },
      {
        key: 'f_keep', label: '结构过滤', sortable: false,
        render: (r) => (r.f_keep === true ? chip('存活', 'pass') : r.f_keep === false ? chip('未过', 'weak') : chip('—')),
      },
      {
        key: 'n_contacts', label: '接触残基', sortable: false,
        render: (r) => (r.murj_contacts ? String(r.murj_contacts).split(/\s+/).length : 0),
      },
    ],
    rows,
    keyField: 'smiles',
    pageSize: 50,
    initialSort: { key: 'total', dir: 'desc' },
    empty: '没有满足当前条件的骨架',
    onRowClick: (row) => openDetail(row),
  });

  /* ---- toolbar ---- */
  const search = el('input', {
    type: 'search', class: 'search', placeholder: '搜索 SMILES / 设计轴 / 接触残基（如 K368）',
    oninput: (event) => { state.q = event.target.value.trim(); apply(); },
  });
  const axisSel = el('select', { onchange: (event) => { state.axis = event.target.value; apply(); } },
    el('option', { value: '' }, '全部设计轴'),
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
    label(tripleBox, '仅三重达标'),
    label(keepBox, '仅结构过滤存活'),
    label(mwInput, 'MW ≤'),
    div({ class: 'grow' }),
    countLabel,
    el('button', { class: 'btn sm', text: '导出当前视图 CSV', onclick: exportCsv }));

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
    countLabel.textContent = `${filtered.length} / ${rows.length} 条`;
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
    toast(`已导出 ${table.rows.length} 行`, 'ok');
  }

  /* ---- detail drawer ---- */
  function openDetail(row) {
    const poseIndex = row.pose_index === null || row.pose_index === undefined ? null : Number(row.pose_index);
    const contacts = String(row.murj_contacts || '').split(/\s+/).filter(Boolean);
    const f = row.f_scores || {};
    const filters = row.f_keep === undefined
      ? div({ class: 'card-note', text: '该骨架不在结构衍生过滤的输入集合内（仅三重达标子集被评估）。' })
      : div({ style: { display: 'grid', gap: '8px' } },
        div({ class: 'split' },
          chip(row.f_keep ? 'F1–F4 全部通过' : '未通过结构过滤', row.f_keep ? 'pass' : 'weak'),
          row.f_anchor_ids ? chip(`锚点 ${row.f_anchor_ids}`, 'info') : null,
          chip(`桥接 ${row.f_bridge ? '是' : '否'}`, row.f_bridge ? 'pass' : ''),
          chip(`跨结构分差 ${fmt(row.f_ens_spread, 2)}`, Math.abs(row.f_ens_spread ?? 99) <= 5 ? 'pass' : 'weak')),
        div({ class: 'kv tiny' },
          el('div', { class: 'k', text: '三结构得分' }),
          el('div', { class: 'v mono-sm', text: Object.entries(f).map(([k, v]) => `${k.replace('Sgl', 'Sgl-')} ${fmt(v, 2)}`).join('   ') }),
          ...(row.f_ens_rmsd !== undefined ? [
            el('div', { class: 'k', text: '重放 RMSD' }),
            el('div', { class: 'v mono-sm', text: `${fmt(row.f_ens_rmsd, 2)} Å` }),
          ] : [])));

    const body = div({ style: { display: 'grid', gap: '14px' } },
      div({ class: 'smiles-box', text: row.smiles }),
      div({ class: 'split' },
        chip(`排名 #${row.rank}`, 'info'),
        chip(`轴 ${row.axis}`, 'mono'),
        chip(`家族 ${row.family}`, 'mono'),
        row.triple ? chip('三重达标', 'pass') : null),
      div({ class: 'grid cols-2', style: { gap: '10px' } },
        div({ class: 'card flat' },
          div({ class: 'card-head' }, el('h3', { text: '性质包络' })),
          kvList([
            ['分子量 MW', `${fmt(row.MW, 1)} Da`],
            ['cLogP', fmt(row.cLogP, 2)],
            ['TPSA', `${fmt(row.TPSA, 1)} Å²`],
            ['氢键供体/受体', `${row.HBD} / ${row.HBA}`],
            ['可旋转键', row.RotB],
            ['环数', row.Rings],
            ['重原子', row.heavy_atoms],
          ])),
        div({ class: 'card flat' },
          div({ class: 'card-head' }, el('h3', { text: '打分分解' })),
          div({ style: { display: 'grid', gap: '6px' } },
            meterRow('总分', row.total, 16, 'var(--accent)'),
            meterRow('位姿得分', row.pose_score, 12, 'var(--accent-2)'),
            meterRow('肽模拟', row.peptide_overlap, 1, '#a78bfa'),
            meterRow('埋藏度', row.buried, 1, '#38bdf8'),
            meterRow('刚性', row.rigidity, 1, '#38a169'),
            meterRow('盐桥', row.salt, 3, '#4f8ef7'),
            meterRow('氢键', row.hb, 12, '#c99a3a'),
            meterRow('疏水', row.hyd, 25, '#8b97a8'),
            meterRow('位阻惩罚', row.clash, 1, 'var(--danger)'),
            meterRow('性质惩罚', row.property_penalty, 1, 'var(--warn)')))),

      div({ class: 'card flat' },
        div({ class: 'card-head' }, el('h3', { text: '结构衍生过滤（F1–F4）' })),
        filters),

      div({ class: 'card flat' },
        div({ class: 'card-head' },
          el('h3', { text: '接触的 MurJ 残基' }),
          span({ class: 'right tiny muted', text: `${contacts.length} 个，4.2 Å 内` })),
        div({ class: 'chips' }, ...contacts.map((c) => span({ class: 'resid', text: c })))),

      div({ class: 'card flat' },
        div({ class: 'card-head' }, el('h3', { text: '构建方式' })),
        kvList([
          ['头帽', row.head], ['环系统', row.rings], ['连接子', row.links], ['尾帽', row.tail],
          ['设计轴长', `${fmt(row.axis_length_A, 2)} Å`],
        ])));

    const footer = [
      el('button', {
        class: 'btn primary',
        text: poseIndex === null ? '该骨架无 3D 位姿' : '在 3D 中查看位姿',
        disabled: poseIndex === null,
        onclick: () => {
          if (poseIndex === null) return;
          close();
          ctx.navigate('structure', { ligand: String(poseIndex) });
        },
      }),
      el('button', {
        class: 'btn', text: '复制 SMILES',
        onclick: async () => {
          try {
            await navigator.clipboard.writeText(row.smiles);
            toast('SMILES 已复制', 'ok', 1500);
          } catch (err) {
            toast('浏览器拒绝了剪贴板访问，请手动选择复制', 'err');
          }
        },
      }),
      el('button', {
        class: 'btn ghost', text: '下载该行 JSON',
        onclick: () => downloadText(`scaffold-${row.rank}.json`, JSON.stringify(row, null, 2), 'application/json'),
      }),
    ];

    const { close } = openDrawer({
      title: `骨架 #${row.rank}`,
      sub: `${row.axis} · ${row.family} · 总分 ${fmt(row.total, 3)}`,
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
      el('h1', { text: '骨架库' }),
      div({ class: 'sub', text: `按总分排序的前 ${rows.length} 个骨架，其中 ${nTriple} 个通过三重达标（位阻 + 类药 + 肽模拟）；全部 ${fsum.n_input ?? 538} 个三重达标里有 ${nKeepAll ?? 32} 个通过结构衍生过滤，本页覆盖 ${nKeepHere} 个。位姿是刚体放置而非对接——最终需下游 MD/MM-PBSA 精修。` })),
    div({ class: 'actions' },
      el('a', { class: 'btn', href: '/files/results/scaffolds_ranked.csv?download=1', text: '下载 CSVs' }),
      el('a', { class: 'btn', href: '/files/results/top_scaffold_poses.sdf?download=1', text: '下载 SDF' }))));

  root.appendChild(div({ class: 'grid cols-4', style: { marginBottom: '12px' } },
    card({ class: 'tight' }, stat({ label: '本页骨架', value: rows.length, unit: '个', small: true, foot: '按总分排名' })),
    card({ class: 'tight' }, stat({ label: '其中三重达标', value: nTriple, unit: '个', small: true, tone: 'accent', foot: '位阻 + 类药 + 肽模拟' })),
    card({ class: 'tight' }, stat({
      label: '结构过滤存活', value: nKeepAll ?? '—', unit: '个', small: true, tone: 'good',
      foot: `全部三重达标中；本页 ${nKeepHere} 个`,
    })),
    card({ class: 'tight' }, stat({ label: '有 3D 位姿', value: nPoses, unit: '个', small: true, foot: 'top_scaffold_poses.sdf' }))));

  root.appendChild(card({ title: '排名与筛选' }, toolbar, table.el));

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
