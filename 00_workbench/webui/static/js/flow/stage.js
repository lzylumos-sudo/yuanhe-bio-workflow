/* flow/stage.js -- one shared layout for all 7 pages:
   left  = 3D molecule view (+ 2D card), right = screening-count funnel + molecule list. */
import { get } from '../api.js';
import { el, div, span } from '../ui.js';
import { funnel, viewerPanel, moleculeTable, molCard } from './widgets.js';

const STATUS = { ready: '已接入', demo: '布局预览', pending: '等待数据' };
const nf = (v) => (v === null || v === undefined ? '—' : Number(v).toLocaleString('en-US'));

const card = (title, sub, ...kids) => div({ class: 'card' },
  el('h3', null, title, sub ? span({ class: 'sub', text: sub }) : null), ...kids);

function head(stage, ctx) {
  return div({ class: 'head' },
    div({ class: 'titles' },
      div({ class: 'eyebrow', text: `STEP ${stage.n} / 7` }),
      el('h1', null, stage.title, ' ', span({ class: `badge ${stage.status}`, text: STATUS[stage.status] })),
      el('p', { text: stage.subtitle })),
    div({ class: 'nav' },
      el('button', { class: 'btn', text: '← 上一步', disabled: stage.n === 1, onclick: () => ctx.go(stage.n - 1) }),
      el('button', { class: 'btn primary', text: stage.n === 7 ? '回到第 1 步' : '下一步 →', onclick: () => ctx.go(stage.n === 7 ? 1 : stage.n + 1) })));
}

function banner(stage) {
  if (stage.status === 'demo') {
    return div({ class: 'banner', text: '布局预览：此阶段尚未接入模型，下面的分子与数字取自第 2 步的真实骨架，仅用于检查版面，不是该阶段的计算结果。' });
  }
  if (stage.status === 'pending') {
    return div({ class: 'banner pending', text: '等待数据：模型接入后，把结果文件放到下方路径，本页会自动显示分子和计数。右上角“布局预览”可先看版面。' });
  }
  return null;
}

function filesCard(stage) {
  if (!stage.files?.length) return null;
  return card('数据接口', stage.contract || '',
    div({ class: 'files' }, ...stage.files.map((f) => div(null,
      span({ class: f.exists ? 'ok' : 'miss', text: f.exists ? '●' : '○' }),
      el('code', { text: f.path }), span({ class: 'miss', text: f.role })))),
    el('p', { class: 'note', text: '约定：分子表至少含 id、smiles 两列，可选 keep 列；漏斗文件为 {"funnel":[{"label":"…","count":N}]}；3D 位姿 SDF 的标题行写分子 id。' }));
}

/* ---------------- step 1: complexes ---------------- */
async function mountComplex(root, ctx) {
  const { stage } = ctx;
  const tip = document.getElementById('tip');
  const panel = await viewerPanel({ ligandHint: false, tip });
  const pep = panel.scene.peptides;
  const shown = Object.fromEntries(pep.map((p) => [p.tag, true]));
  const list = div({ class: 'cx' });
  const draw = () => {
    list.replaceChildren(...pep.map((p) => el('button', {
      class: ['cx-item', shown[p.tag] && 'on'],
      onclick: () => {
        shown[p.tag] = !shown[p.tag];
        panel.viewer.setOptions({ peptideFilter: { ...shown } });
        draw();
      },
    },
    el('span', { class: 'sw', style: { background: p.color } }),
    div(null, el('b', { text: p.tag }), ' ', el('small', { text: `${p.pdbid} · ${p.n_residues} 残基` })),
    span({ class: shown[p.tag] ? 'yes' : 'no', text: shown[p.tag] ? '显示' : '隐藏' }),
    div({ class: 'seq', text: p.seq }))));
  };
  draw();
  const sub = panel.scene.subpockets || [];
  const chips = div({ class: 'chips' }, ...sub.map((s) => span({ class: 'chip' },
    el('i', { style: { display: 'inline-block', width: '8px', height: '8px', borderRadius: '50%', background: s.color, marginRight: '5px' } }),
    s.name.split('_')[0], el('b', { text: `${Math.round(s.volume || 0)} Å³` }))));

  root.append(
    head(stage, ctx),
    div({ class: 'grid' },
      div({ class: 'col' }, div({ class: 'card flush' }, panel.root)),
      div({ class: 'col' },
        card('筛选数量变化', '界面 → 口袋 → 药效团', funnel(stage.funnel), chips),
        card('复合物', '点击切换显示', list),
        card('导入文件', '结构来自 PDB 坐标文件', div({ class: 'files' },
          ...pep.map((p) => div(null, span({ class: 'ok', text: '●' }), el('code', { text: `${p.pdbid.toLowerCase()}.pdb` }), span({ class: 'miss', text: p.tag })))),
        el('p', { class: 'note', text: '三个复合物按共享 CA 叠合到 9NU8 坐标系；绿/橙/黄为三条裂解肽，彩色点云为口袋。' })))));
  return () => panel.destroy();
}

/* ---------------- step 7 summary ---------------- */
function pipelineFunnel(wf) {
  const [, s2] = wf.stages;
  const f2 = s2.funnel;
  const rows = [
    { label: '枚举骨架', count: f2[0]?.count, note: '第 2 步' },
    { label: '类药骨架', count: f2[2]?.count },
    { label: '结构过滤存活', count: f2[3]?.count },
  ];
  for (const s of wf.stages.slice(2)) {
    rows.push({ label: `${s.n}. ${s.title}`, count: s.out, note: s.status === 'pending' ? '等待数据' : s.status === 'demo' ? '布局预览' : null });
  }
  return rows;
}

/* ---------------- steps 2-7 ---------------- */
async function mountMolecules(root, ctx) {
  const { stage, demo } = ctx;
  const tip = document.getElementById('tip');
  const q = demo ? '?demo=1' : '';
  const [panel, mols] = await Promise.all([
    viewerPanel({ tip }),
    get(`/api/stage/${stage.key}/molecules${q}`),
  ]);
  const mcard = molCard();
  const caption = div({ class: 'note' });
  let token = 0;

  async function pick(row) {
    mcard.show(row);
    if (!row) { panel.setLigand(null); return; }
    const t = ++token;
    caption.textContent = '';
    const sep = q ? '&' : '?';
    const r = await get(`/api/stage/${stage.key}/pose${q}${sep}id=${encodeURIComponent(row.id)}`);
    if (t !== token) return;
    panel.setLigand(r.ligand, row.id);
    caption.textContent = r.ligand ? '' : (r.reason || '');
  }

  const rightTop = stage.n === 7
    ? card('全流程漏斗', '各阶段保留的分子数', funnel(pipelineFunnel(ctx.wf)))
    : card('筛选数量变化', stage.status === 'demo' ? '版面示意' : null, funnel(stage.funnel),
      stage.extras?.independent
        ? div({ class: 'chips' }, ...stage.extras.independent.map((c) => span({ class: 'chip' }, c.label, el('b', { text: nf(c.count) }))))
        : null);

  let listBody;
  let table = null;
  if (mols.n) {
    table = moleculeTable(mols, { onSelect: pick });
    listBody = table.root;
  } else {
    listBody = div({ class: 'empty', text: '暂无分子。接入后在此显示列表，选中一行即可在左侧 3D 视图中查看。' });
  }

  const grid = div({ class: 'grid' },
    div({ class: 'col' },
      div({ class: 'card flush' }, panel.root),
      card('分子', stage.status === 'pending' ? null : '2D 结构与属性', mcard.root, caption)),
    div({ class: 'col' },
      banner(stage),
      rightTop,
      card('分子列表', `${nf(mols.n)} 个 · 点击行查看`, listBody),
      filesCard(stage)));
  root.append(head(stage, ctx), grid);

  if (table) {
    const rows = table.rows();
    const first = rows.find((r) => r.keep) || rows[0];
    table.select(first);
  } else {
    mcard.show(null);
  }
  return () => panel.destroy();
}

export async function mountStage(root, ctx) {
  if (ctx.stage.key === 'complex') return mountComplex(root, ctx);
  return mountMolecules(root, ctx);
}
