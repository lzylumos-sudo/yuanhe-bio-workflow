/* pages/structure.js -- interactive 3D view of the complexes, the cavity and the
   designed scaffolds (custom Canvas renderer, no WebGL, no external libs). */

import { api } from '../api.js';
import {
  el, div, span, card, chip, kvList, fmt, fmtInt, loadingPane, toast,
} from '../ui.js';
import { StructureViewer, SUB_COLORS } from '../viewer.js';

const SUB_LABEL = {
  S1_basic_bowl: 'S1 碱性碗',
  S2_carboxylate_gate: 'S2 羧酸闸门',
  S3_TM2_TM7_groove: 'S3 TM2/TM7 沟',
  S4_TM7_floor: 'S4 TM7 底板',
};

export async function mount(root, ctx) {
  const params = ctx.params || {};
  let ref = params.ref && ['SglM', 'SglPP7', 'SglCJ3'].includes(params.ref) ? params.ref : 'SglCJ3';
  let ligandIndex = params.ligand !== undefined ? params.ligand : '0';
  if (ligandIndex === 'none') ligandIndex = null;

  root.appendChild(loadingPane('正在把三个复合物叠到同一坐标系…'));
  const [poses, pocketData] = await Promise.all([
    api.poses(),
    api.pocket().catch(() => ({ pocket: null })),
  ]);
  const top = poses.top || [];
  const pocket = pocketData.pocket || {};

  /* ---------- layout ---------- */
  const viewerHost = div({ class: 'viewer-shell viewer-panel' });
  const overlay = div({ class: 'viewer-overlay' });
  viewerHost.appendChild(overlay);

  const tooltip = document.getElementById('tooltip');
  const infoCard = card({ class: 'tight' }, div({ class: 'card-note', text: '把鼠标移到残基或配体原子上查看标签；点击可选中。' }));

  const controls = div({ class: 'viewer-controls' });

  root.replaceChildren();
  root.appendChild(div({ class: 'page-head' },
    div(null,
      el('h1', { text: '3D 结构' }),
      div({ class: 'sub', text: '三种噬菌体肽通过共享 CA 做 Kabsch 叠合后放进同一个参考坐标系；骨架位姿由算法放置在该坐标系中，可直接与肽的空间占据比较。' })),
    div({ class: 'actions' },
      el('button', { class: 'btn', text: '复位视角', onclick: () => viewer.fit() }),
      el('a', { class: 'btn', href: '/files/results/diverse_scaffolds.sdf?download=1', text: '下载位姿 SDF' }))));

  root.appendChild(div({ class: 'viewer-layout' },
    card({ class: 'tight' }, controls),
    div(null, viewerHost,
      div({ class: 'card-note', style: { marginTop: '6px' } },
        '拖动旋转 · 滚轮缩放 · Shift/右键拖动平移 · 双击复位。',
        span({ id: 'pocket-status', class: 'muted' })))));

  /* ---------- viewer ---------- */
  const viewer = new StructureViewer(viewerHost, {
    onHover: (item, event) => {
      if (!item) { tooltip.hidden = true; return; }
      tooltip.hidden = false;
      tooltip.textContent = item.label;
      tooltip.style.left = `${(event?.clientX ?? 0) + 14}px`;
      tooltip.style.top = `${(event?.clientY ?? 0) + 12}px`;
    },
    onPick: (item) => {
      if (item.kind === 'residue' || item.kind === 'atom') {
        const resid = item.resid || Number(String(item.label).match(/(\d+)/)?.[1]);
        const hot = scene?.hotspots?.find((h) => Number(h.resid) === Number(resid));
        infoCard.replaceChildren(
          div({ class: 'card-head' },
            el('h3', { text: `已选中 ${item.label}` }),
            hot ? span({ class: 'right chip mono', text: `热点 ${fmt(hot.score, 3)}` }) : null),
          hot ? kvList([
            ['dSASA', `${fmt(hot.dSASA, 1)} Å²`],
            ['证据', hot.evidence || '几何证据'],
          ]) : div({ class: 'card-note', text: '该残基不在前 24 个热点内。' }),
          div({ class: 'tiny muted', style: { marginTop: '6px' }, text: '提示：热点球体已按 score 着色并标注为橙色（有遗传学证据）。' }));
      } else if (item.kind === 'hotspot' && item.hotspot) {
        const h = item.hotspot;
        infoCard.replaceChildren(
          div({ class: 'card-head' }, el('h3', { text: `热点 ${h.one}${h.resid}` })),
          kvList([
            ['score', fmt(h.score, 3)],
            ['dSASA', `${fmt(h.dSASA, 1)} Å²`],
            ['证据', h.evidence || '几何证据'],
          ]));
      }
    },
  });

  /* ---------- scene loading ---------- */
  let scene = null;
  let pocketPoints = null;
  let pocketTimer = null;

  async function loadScene() {
    overlay.replaceChildren(div({ class: 'line', text: '载入结构中…' }));
    scene = await api.scene({ ref, ligand: ligandIndex === null ? undefined : ligandIndex, radius: 13 });
    viewer.setScene(scene);
    paintOverlay();
    updateInfo();
    if (viewer.options.showPocket) ensurePocket();
  }

  function paintOverlay() {
    if (!scene) return;
    const structure = scene.structure || {};
    const ligand = scene.ligand;
    overlay.replaceChildren(
      div({ class: 'line', text: `${scene.ref_pdbid} · ${ref} · ${fmt(structure.resolution, 2)} Å` }),
      div({ class: 'line', text: `参考坐标系 ${scene.ref}（MurJ 5–500）` }),
      ligand ? div({ class: 'line', text: `配体 ${ligand.smiles}` }) : div({ class: 'line', text: '无配体（仅显示腔体）' }),
    );
  }

  function updateInfo() {
    if (!scene) return;
    const ligand = scene.ligand;
    const contacts = ligandContacts(scene);
    const infoChildren = [];
    if (ligand) {
      const poseRow = top[Number(ligandIndex)];
      infoChildren.push(
        div({ class: 'card-head' },
          el('h3', { text: '当前位姿' }),
          poseRow ? span({ class: 'right chip mono', text: poseRow.legend.split('score=')[1] || '' }) : null),
        div({ class: 'smiles-box', style: { fontSize: '11.5px' }, text: ligand.smiles }),
        kvList([
          ['重原子', ligand.atoms.length],
          ['接触残基', `${contacts.length} 个（4.2 Å 内）`],
        ]),
        div({ class: 'chips', style: { marginTop: '6px' } },
          ...contacts.slice(0, 24).map((c) => span({ class: 'resid', text: c }))),
        div({ class: 'split', style: { marginTop: '8px' } },
          el('button', {
            class: 'btn sm', text: '骨架详情 →',
            onclick: () => ctx.navigate('scaffolds', { focus: ligand.smiles }),
          })));
    } else {
      infoChildren.push(div({ class: 'card-note', text: '未选择配体：视图现在只显示受体、肽路径、药效团和（可选的）腔体点云。' }));
    }
    infoChildren.push(div({ class: 'hr' }));
    infoChildren.push(div({ class: 'tiny muted', text: '三种肽的序列（PDB 中建模部分）' }));
    infoChildren.push(div({ style: { display: 'grid', gap: '3px', marginTop: '4px' } },
      ...(scene.peptides || []).map((p) => div({ class: 'split' },
        span({ class: 'swatch', style: { background: p.color, width: '8px', height: '8px', borderRadius: '2px' } }),
        span({ class: 'tiny muted', text: p.pdbid }),
        span({ class: 'mono-sm right', style: { marginLeft: 'auto', fontSize: '10px' }, text: p.seq })))));
    infoCard.replaceChildren(...infoChildren);
  }

  /* ---------- pocket cloud ---------- */
  function setPocketStatus(text, busy = false) {
    const node = document.getElementById('pocket-status');
    if (!node) return;
    node.replaceChildren(busy ? span({ class: 'spinner', style: { display: 'inline-block', marginRight: '6px' } }) : span(),
      span({ text }));
  }

  async function ensurePocket() {
    if (pocketPoints) { viewer.setPocket(pocketPoints); return; }
    try {
      const data = await api.pocketPoints(0.8);
      if (data.status === 'ready') {
        pocketPoints = data;
        viewer.setPocket(data);
        setPocketStatus(`腔体点云：${fmtInt(data.n_voxels)} 体素 / ${fmtInt(data.volume)} Å³`);
        return;
      }
      setPocketStatus('首次计算腔体体素（约 15 s）…', true);
      clearTimeout(pocketTimer);
      pocketTimer = setTimeout(ensurePocket, 2500);
    } catch (err) {
      setPocketStatus(`腔体点云计算失败：${err.message}`, false);
    }
  }

  /* ---------- controls ---------- */
  function paintControls() {
    const o = viewer.options;
    controls.replaceChildren(
      div({ class: 'ctrl-group' },
        div({ class: 'cap', text: '参考结构' }),
        div({ class: 'seg' },
          ...['SglCJ3', 'SglM', 'SglPP7'].map((tag) => el('button', {
            class: tag === ref ? 'active' : '',
            text: { SglCJ3: '9NU8', SglM: '9NU4', SglPP7: '9NU5' }[tag],
            title: tag,
            onclick: async () => {
              if (ref === tag) return;
              ref = tag;
              paintControls();
              await loadScene();
            },
          }))),
        div({ class: 'tiny muted', text: '叠合由共享 Cα 的 Kabsch 变换完成，肽与药效团跟随变换。' })),

      div({ class: 'ctrl-group' },
        div({ class: 'cap', text: '配体位姿（前 25）' }),
        el('select', {
          onchange: async (event) => {
            ligandIndex = event.target.value === 'none' ? null : event.target.value;
            paintOverlay();
            await loadScene();
          },
        },
          el('option', { value: 'none', selected: ligandIndex === null }, '无配体（只看口袋）'),
          ...top.map((m, idx) => el('option', {
            value: String(idx), selected: String(idx) === String(ligandIndex),
          }, `#${idx + 1} ${m.smiles.length > 34 ? `${m.smiles.slice(0, 34)}…` : m.smiles}`)))),

      div({ class: 'ctrl-group' },
        div({ class: 'cap', text: '显示' }),
        toggle('受体轨迹', o.showTrace, (v) => viewer.setOptions({ showTrace: v })),
        toggle('衬里侧链', o.showSidechains, (v) => viewer.setOptions({ showSidechains: v })),
        toggle('配体位姿', o.showLigand, (v) => viewer.setOptions({ showLigand: v })),
        toggle('药效团特征', o.showFeatures, (v) => viewer.setOptions({ showFeatures: v })),
        toggle('热点残基', o.showHotspots, (v) => viewer.setOptions({ showHotspots: v })),
        toggle('标签', o.showLabels, (v) => viewer.setOptions({ showLabels: v })),
        toggle('自转', o.spin, (v) => viewer.setOptions({ spin: v }))),

      div({ class: 'ctrl-group' },
        div({ class: 'cap', text: '噬菌体肽' }),
        ...(scene?.peptides || []).map((p) => toggle(`${p.pdbid} ${p.tag}`, o.peptideFilter[p.tag] !== false, (v) => {
          o.peptideFilter[p.tag] = v;
          viewer.setOptions({ showPeptides: true });
        }))),

      div({ class: 'ctrl-group' },
        div({ class: 'cap', text: '埋藏腔体点云' }),
        toggle('显示腔体', o.showPocket, (v) => {
          viewer.setOptions({ showPocket: v });
          if (v) ensurePocket();
          else setPocketStatus('已隐藏');
        }),
        div({ class: 'seg' },
          ...[['chem', '化学'], ['sub', '亚口袋']].map(([mode, label]) => el('button', {
            class: o.pocketMode === mode ? 'active' : '',
            text: label,
            onclick: () => {
              viewer.setOptions({ pocketMode: mode });
              paintControls();
            },
          })),)),

      div({ class: 'ctrl-group' },
        div({ class: 'cap', text: '图例' }),
        div({ class: 'legend', style: { display: 'grid', gap: '4px' } },
          ...Object.entries(SUB_LABEL).map(([key, label]) => div({ class: 'item' },
            span({ class: 'swatch', style: { background: SUB_COLORS[key] } }), span({ class: 'tiny', text: label })))),
        div({ class: 'legend', style: { display: 'grid', gap: '4px', marginTop: '6px' } },
          ...(scene?.peptides || []).map((p) => div({ class: 'item' },
            span({ class: 'swatch', style: { background: p.color, borderRadius: '50%' } }), span({ class: 'tiny', text: `${p.tag}（${p.pdbid}）` }))),
          div({ class: 'item' }, span({ class: 'swatch', style: { background: '#ff8a5b', borderRadius: '50%' } }), span({ class: 'tiny', text: '热点（有遗传学证据）' })))),

      infoCard);
  }

  function toggle(labelText, checked, onchange) {
    return el('label', { class: 'check' },
      el('input', { type: 'checkbox', checked, onchange: (event) => onchange(event.target.checked) }),
      labelText);
  }

  await loadScene();
  paintControls();

  if (params.focus) {
    const rid = Number(params.focus);
    const res = (scene.receptor.residues || []).find((r) => Number(r.resid) === rid);
    const hot = (scene.hotspots || []).find((h) => Number(h.resid) === rid);
    if (res) {
      infoCard.replaceChildren(
        div({ class: 'card-head' },
          el('h3', { text: `定位残基 ${res.one}${res.resid}` }),
          hot ? span({ class: 'right chip mono', text: `热点 ${fmt(hot.score, 3)}` }) : null),
        kvList([
          ['亚口袋', res.sub ? SUB_LABEL[res.sub] || res.sub : '—'],
          ['dSASA', hot ? `${fmt(hot.dSASA, 1)} Å²` : '—'],
          ['证据', hot?.evidence || '几何证据'],
        ]),
        div({ class: 'tiny muted', style: { marginTop: '6px' }, text: '该残基已在视图中以热点球体标出；旋转视图可查看它与配体的相对位置。' }),
        div({ class: 'split', style: { marginTop: '8px' } },
          el('button', { class: 'btn sm', text: '打开界面热点表 →', onclick: () => ctx.navigate('hotspots') })));
    }
  }

  if (params.pocket === '1') {
    viewer.setOptions({ showPocket: true });
    ensurePocket();
  }

  return () => {
    clearTimeout(pocketTimer);
    viewer.destroy();
    tooltip.hidden = true;
  };
}

/* residues of the receptor within 4.2 A of the ligand */
function ligandContacts(scene) {
  const ligand = scene.ligand;
  if (!ligand) return [];
  const set = new Map();
  for (const atom of scene.receptor.atoms || []) {
    for (const la of ligand.atoms) {
      const d = Math.hypot(atom.p[0] - la.p[0], atom.p[1] - la.p[1], atom.p[2] - la.p[2]);
      if (d <= 4.2) {
        set.set(atom.resid, `${atom.one}${atom.resid}`);
        break;
      }
    }
  }
  return [...set.values()];
}
