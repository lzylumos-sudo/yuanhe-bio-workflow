/* pages/structure.js -- interactive 3D view of the complexes, the cavity and the
   designed scaffolds (custom Canvas renderer, no WebGL, no external libs). */

import { api } from '../api.js';
import {
  el, div, span, card, chip, kvList, fmt, fmtInt, loadingPane, toast,
} from '../ui.js';
import { StructureViewer, SUB_COLORS } from '../viewer.js';

const SUB_LABEL = {
  S1_basic_bowl: 'S1 Basic bowl',
  S2_carboxylate_gate: 'S2 Carboxylate gate',
  S3_TM2_TM7_groove: 'S3 TM2/TM7 Groove',
  S4_TM7_floor: 'S4 TM7 Floor',
};

export async function mount(root, ctx) {
  const params = ctx.params || {};
  let ref = params.ref && ['SglM', 'SglPP7', 'SglCJ3'].includes(params.ref) ? params.ref : 'SglCJ3';
  let ligandIndex = params.ligand !== undefined ? params.ligand : '0';
  if (ligandIndex === 'none') ligandIndex = null;

  root.appendChild(loadingPane('Aligning three complexes into one coordinate frame...'));
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
  const infoCard = card({ class: 'tight' }, div({ class: 'card-note', text: 'Hover over a residue or ligand atom for its label; click to select.' }));

  const controls = div({ class: 'viewer-controls' });

  root.replaceChildren();
  root.appendChild(div({ class: 'page-head' },
    div(null,
      el('h1', { text: '3D Structure' }),
      div({ class: 'sub', text: 'Three phage peptides use shared CA for Kabsch alignment into a common reference frame. Scaffold poses occupy this frame for direct spatial comparison.' })),
    div({ class: 'actions' },
      el('button', { class: 'btn', text: 'Reset view', onclick: () => viewer.fit() }),
      el('a', { class: 'btn', href: '/files/results/diverse_scaffolds.sdf?download=1', text: 'Download pose SDF' }))));

  root.appendChild(div({ class: 'viewer-layout' },
    card({ class: 'tight' }, controls),
    div(null, viewerHost,
      div({ class: 'card-note', style: { marginTop: '6px' } },
        'Drag to rotate · Scroll to zoom · Shift/Right-drag to pan · Double-click to reset.',
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
            el('h3', { text: `Selected ${item.label}` }),
            hot ? span({ class: 'right chip mono', text: `Hotspot ${fmt(hot.score, 3)}` }) : null),
          hot ? kvList([
            ['dSASA', `${fmt(hot.dSASA, 1)} Å²`],
            ['Evidence', hot.evidence || 'Geometric evidence'],
          ]) : div({ class: 'card-note', text: 'This residue is not among the top 24 hotspots.' }),
          div({ class: 'tiny muted', style: { marginTop: '6px' }, text: 'Hotspot spheres are colored by score and marked orange when genetic evidence is available.' }));
      } else if (item.kind === 'hotspot' && item.hotspot) {
        const h = item.hotspot;
        infoCard.replaceChildren(
          div({ class: 'card-head' }, el('h3', { text: `Hotspot ${h.one}${h.resid}` })),
          kvList([
            ['score', fmt(h.score, 3)],
            ['dSASA', `${fmt(h.dSASA, 1)} Å²`],
            ['Evidence', h.evidence || 'Geometric evidence'],
          ]));
      }
    },
  });

  /* ---------- scene loading ---------- */
  let scene = null;
  let pocketPoints = null;
  let pocketTimer = null;

  async function loadScene() {
    overlay.replaceChildren(div({ class: 'line', text: 'Loading structures...' }));
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
      div({ class: 'line', text: `Reference coordinate frame ${scene.ref}(MurJ 5–500)` }),
      ligand ? div({ class: 'line', text: `Ligand ${ligand.smiles}` }) : div({ class: 'line', text: 'No ligand (cavity only)' }),
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
          el('h3', { text: 'Current pose' }),
          poseRow ? span({ class: 'right chip mono', text: poseRow.legend.split('score=')[1] || '' }) : null),
        div({ class: 'smiles-box', style: { fontSize: '11.5px' }, text: ligand.smiles }),
        kvList([
          ['Heavy atoms', ligand.atoms.length],
          ['Contact residues', `${contacts.length} items (4.2 Å within)`],
        ]),
        div({ class: 'chips', style: { marginTop: '6px' } },
          ...contacts.slice(0, 24).map((c) => span({ class: 'resid', text: c }))),
        div({ class: 'split', style: { marginTop: '8px' } },
          el('button', {
            class: 'btn sm', text: 'Scaffold details →',
            onclick: () => ctx.navigate('scaffolds', { focus: ligand.smiles }),
          })));
    } else {
      infoChildren.push(div({ class: 'card-note', text: 'No ligand selected: showing receptor, peptide traces, pharmacophores and optional cavity points.' }));
    }
    infoChildren.push(div({ class: 'hr' }));
    infoChildren.push(div({ class: 'tiny muted', text: 'Sequences of three peptides (PDB modeled portion)' }));
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
        setPocketStatus(`Cavity point cloud:${fmtInt(data.n_voxels)} Voxels / ${fmtInt(data.volume)} Å³`);
        return;
      }
      setPocketStatus('Computing cavity voxels for the first time (approximately 15 s)…', true);
      clearTimeout(pocketTimer);
      pocketTimer = setTimeout(ensurePocket, 2500);
    } catch (err) {
      setPocketStatus(`Cavity calculation failed:${err.message}`, false);
    }
  }

  /* ---------- controls ---------- */
  function paintControls() {
    const o = viewer.options;
    controls.replaceChildren(
      div({ class: 'ctrl-group' },
        div({ class: 'cap', text: 'Reference structure' }),
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
        div({ class: 'tiny muted', text: 'Alignment uses shared Cα associated Kabsch transformation; peptides and pharmacophores follow the same transform.' })),

      div({ class: 'ctrl-group' },
        div({ class: 'cap', text: 'Ligand poses (top 25)' }),
        el('select', {
          onchange: async (event) => {
            ligandIndex = event.target.value === 'none' ? null : event.target.value;
            paintOverlay();
            await loadScene();
          },
        },
          el('option', { value: 'none', selected: ligandIndex === null }, 'No ligand (pocket only)'),
          ...top.map((m, idx) => el('option', {
            value: String(idx), selected: String(idx) === String(ligandIndex),
          }, `#${idx + 1} ${m.smiles.length > 34 ? `${m.smiles.slice(0, 34)}…` : m.smiles}`)))),

      div({ class: 'ctrl-group' },
        div({ class: 'cap', text: 'Show' }),
        toggle('Receptor trace', o.showTrace, (v) => viewer.setOptions({ showTrace: v })),
        toggle('Lining side chains', o.showSidechains, (v) => viewer.setOptions({ showSidechains: v })),
        toggle('Ligand pose', o.showLigand, (v) => viewer.setOptions({ showLigand: v })),
        toggle('Pharmacophore features', o.showFeatures, (v) => viewer.setOptions({ showFeatures: v })),
        toggle('Hotspot residues', o.showHotspots, (v) => viewer.setOptions({ showHotspots: v })),
        toggle('Labels', o.showLabels, (v) => viewer.setOptions({ showLabels: v })),
        toggle('Auto-rotate', o.spin, (v) => viewer.setOptions({ spin: v }))),

      div({ class: 'ctrl-group' },
        div({ class: 'cap', text: 'Phage peptide' }),
        ...(scene?.peptides || []).map((p) => toggle(`${p.pdbid} ${p.tag}`, o.peptideFilter[p.tag] !== false, (v) => {
          o.peptideFilter[p.tag] = v;
          viewer.setOptions({ showPeptides: true });
        }))),

      div({ class: 'ctrl-group' },
        div({ class: 'cap', text: 'Buried cavity points' }),
        toggle('Show cavity', o.showPocket, (v) => {
          viewer.setOptions({ showPocket: v });
          if (v) ensurePocket();
          else setPocketStatus('Hidden');
        }),
        div({ class: 'seg' },
          ...[['chem', 'Chemistry'], ['sub', 'Subpocket']].map(([mode, label]) => el('button', {
            class: o.pocketMode === mode ? 'active' : '',
            text: label,
            onclick: () => {
              viewer.setOptions({ pocketMode: mode });
              paintControls();
            },
          })),)),

      div({ class: 'ctrl-group' },
        div({ class: 'cap', text: 'Legend' }),
        div({ class: 'legend', style: { display: 'grid', gap: '4px' } },
          ...Object.entries(SUB_LABEL).map(([key, label]) => div({ class: 'item' },
            span({ class: 'swatch', style: { background: SUB_COLORS[key] } }), span({ class: 'tiny', text: label })))),
        div({ class: 'legend', style: { display: 'grid', gap: '4px', marginTop: '6px' } },
          ...(scene?.peptides || []).map((p) => div({ class: 'item' },
            span({ class: 'swatch', style: { background: p.color, borderRadius: '50%' } }), span({ class: 'tiny', text: `${p.tag}(${p.pdbid})` }))),
          div({ class: 'item' }, span({ class: 'swatch', style: { background: '#ff8a5b', borderRadius: '50%' } }), span({ class: 'tiny', text: 'Hotspot with genetic evidence' })))),

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
          el('h3', { text: `Locate residue ${res.one}${res.resid}` }),
          hot ? span({ class: 'right chip mono', text: `Hotspot ${fmt(hot.score, 3)}` }) : null),
        kvList([
          ['Subpocket', res.sub ? SUB_LABEL[res.sub] || res.sub : '—'],
          ['dSASA', hot ? `${fmt(hot.dSASA, 1)} Å²` : '—'],
          ['Evidence', hot?.evidence || 'Geometric evidence'],
        ]),
        div({ class: 'tiny muted', style: { marginTop: '6px' }, text: 'This residue is marked by a hotspot sphere; rotate to inspect its position relative to the ligand.' }),
        div({ class: 'split', style: { marginTop: '8px' } },
          el('button', { class: 'btn sm', text: 'Open hotspot table →', onclick: () => ctx.navigate('hotspots') })));
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
