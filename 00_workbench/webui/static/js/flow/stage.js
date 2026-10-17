/* flow/stage.js -- one shared layout for all 7 pages:
   left  = 3D molecule view (+ 2D card), right = screening-count funnel + molecule list. */
import { get } from '../api.js';
import { el, div, span } from '../ui.js';
import { funnel, viewerPanel, moleculeTable, molCard } from './widgets.js';

const STATUS = { ready: 'Connected', demo: 'Layout preview', pending: 'Waiting for data' };
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
      el('button', { class: 'btn', text: '← Previous step', disabled: stage.n === 1, onclick: () => ctx.go(stage.n - 1) }),
      el('button', { class: 'btn primary', text: stage.n === 7 ? 'Return to step 1 step' : 'Next step →', onclick: () => ctx.go(stage.n === 7 ? 1 : stage.n + 1) })));
}

function banner(stage) {
  if (stage.status === 'demo') {
    return div({ class: 'banner', text: 'Layout preview: this stage has no model connection yet. Molecules and counts below come from step 2 and illustrate layout only; they are not outputs of this stage.' });
  }
  if (stage.status === 'pending') {
    return div({ class: 'banner pending', text: 'Waiting for data: after connecting a model, place results at the paths below. Use Layout preview to inspect the page.' });
  }
  return null;
}

function filesCard(stage) {
  if (!stage.files?.length) return null;
  return card('Data interface', stage.contract || '',
    div({ class: 'files' }, ...stage.files.map((f) => div(null,
      span({ class: f.exists ? 'ok' : 'miss', text: f.exists ? '●' : '○' }),
      el('code', { text: f.path }), span({ class: 'miss', text: f.role })))),
    el('p', { class: 'note', text: 'Contract: molecule tables must include id, smiles two columns; optional keep columns; funnel file format: {"funnel":[{"label":"…","count":N}]}; 3D Pose SDF title line contains molecule id. ' }));
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
    div(null, el('b', { text: p.tag }), ' ', el('small', { text: `${p.pdbid} · ${p.n_residues} Residue` })),
    span({ class: shown[p.tag] ? 'yes' : 'no', text: shown[p.tag] ? 'Show' : 'Hide' }),
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
        card('Screening counts', 'Interface → Pocket → Pharmacophore', funnel(stage.funnel), chips),
        card('Complex', 'Click to toggle', list),
        card('Import file', 'Structures from PDB coordinate files', div({ class: 'files' },
          ...pep.map((p) => div(null, span({ class: 'ok', text: '●' }), el('code', { text: `${p.pdbid.toLowerCase()}.pdb` }), span({ class: 'miss', text: p.tag })))),
        el('p', { class: 'note', text: 'Three complexes aligned through shared CA onto the 9NU8 coordinate frame; green,/orange/and yellow denote the three peptides; colored point clouds show the pocket.' })))));
  return () => panel.destroy();
}

/* ---------------- step 7 summary ---------------- */
function pipelineFunnel(wf) {
  const [, s2] = wf.stages;
  const f2 = s2.funnel;
  const rows = [
    { label: 'Enumerated scaffolds', count: f2[0]?.count, note: 'Number 2 step' },
    { label: 'Drug-like scaffolds', count: f2[2]?.count },
    { label: 'Structure filters passed', count: f2[3]?.count },
  ];
  for (const s of wf.stages.slice(2)) {
    rows.push({ label: `${s.n}. ${s.title}`, count: s.out, note: s.status === 'pending' ? 'Waiting for data' : s.status === 'demo' ? 'Layout preview' : null });
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
    ? card('Workflow funnel', 'Molecules retained at each stage', funnel(pipelineFunnel(ctx.wf)))
    : card('Screening counts', stage.status === 'demo' ? 'Layout illustration' : null, funnel(stage.funnel),
      stage.extras?.independent
        ? div({ class: 'chips' }, ...stage.extras.independent.map((c) => span({ class: 'chip' }, c.label, el('b', { text: nf(c.count) }))))
        : null);

  let listBody;
  let table = null;
  if (mols.n) {
    table = moleculeTable(mols, { onSelect: pick });
    listBody = table.root;
  } else {
    listBody = div({ class: 'empty', text: 'No molecules yet. Select a row to view its structure in the left-hand 3D viewer.' });
  }

  const grid = div({ class: 'grid' },
    div({ class: 'col' },
      div({ class: 'card flush' }, panel.root),
      card('Molecule', stage.status === 'pending' ? null : '2D Structure and properties', mcard.root, caption)),
    div({ class: 'col' },
      banner(stage),
      rightTop,
      card('Molecule list', `${nf(mols.n)} items · Click a row to view`, listBody),
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
