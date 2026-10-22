// Controls write the same frozen backend configuration used by real tasks.
export function screeningControls(parent, cfg, el) {
  const engine = cfg.research?.predocking?.engine;
  if (!engine) return;
  const section = el('fieldset'); section.className = 'panel';
  section.append(el('legend', 'Pre-Docking Screening'));
  section.append(el('p', 'Predictions support research decisions. Missing results require review or stop the stage. Computed synthesis routes still require experimental validation.'));
  const input = (path, title, type, values) => {
    const label = el('label', title);
    const control = el(type === 'select' ? 'select' : 'input');
    control.dataset.screening = path;
    const value = path.split('.').reduce((v,k) => v?.[k], engine);
    if (type === 'select') {
      for (const item of values) { const option = el('option', item); option.value = item; control.append(option); }
      control.value = value;
    } else {
      control.type = type;
      if (type === 'checkbox') control.checked = value; else control.value = value ?? '';
      if (type === 'number') { control.step = 'any'; control.required = true; }
    }
    label.append(control); section.append(label);
  };
  input('mode', 'Screening mode', 'select', ['balanced', 'discovery', 'delivery']);
  input('keep', 'Maximum molecules retained', 'number');
  input('missing_policy', 'Missing model results', 'select', ['review', 'abort']);
  input('structure_enabled', 'Enable structural alerts (disable only for reviewed research policies)', 'checkbox');
  input('catalog_actions.PAINS', 'PAINS alert action', 'select', ['review', 'warn', 'exclude']);
  input('catalog_actions.BRENK', 'Brenk alert action', 'select', ['review', 'warn', 'exclude']);
  input('toxicity.python', 'ADMET-AI environment Python', 'text');
  input('toxicity.device', 'AI inference device', 'select', ['cpu', 'auto', 'cuda']);
  input('toxicity.batch_size', 'AI batch size', 'number');
  input('toxicity.cpu_threads', 'AI CPU threads', 'number');
  for (const endpoint of ['AMES', 'hERG', 'DILI', 'LD50_Zhu']) {
    const label = el('label', 'Predict '+endpoint);
    const control = el('input'); control.type = 'checkbox'; control.dataset.endpoint = endpoint;
    control.checked = engine.toxicity.endpoints.includes(endpoint); label.append(control); section.append(label);
    input('toxicity.thresholds.'+endpoint, endpoint+' review threshold'+(endpoint==='LD50_Zhu'?' / log10(1/(mol/kg))':' / probability'), 'number');
  }
  input('synthesis.sa_max', 'SA maximum / 1–10, higher is harder', 'number');
  input('synthesis.sc_max', 'SCScore maximum / 1–5, higher is harder', 'number');
  input('synthesis.scscore_path', 'Official SCScore 1024bool JSON.GZ weights', 'text');
  input('synthesis.scscore_sha256', 'SCScore weight SHA256', 'text');
  input('retro_max_candidates', 'Maximum priority molecules for retrosynthesis', 'number');
  section.append(el('p', 'Route models, search time/depth and per-rule overrides remain editable in the pre-docking JSON field. A complete validated computed route is required.'));
  parent.append(section);
}

export function readScreening(parent, cfg) {
  const engine = cfg.research?.predocking?.engine;
  if (!engine) return;
  for (const control of parent.querySelectorAll('[data-screening]')) {
    const keys = control.dataset.screening.split('.');
    const target = keys.slice(0,-1).reduce((v,k)=>v[k], engine);
    target[keys.at(-1)] = control.type === 'checkbox' ? control.checked : control.type === 'number' ? Number(control.value) : control.value;
  }
  engine.toxicity.endpoints = [...parent.querySelectorAll('[data-endpoint]:checked')].map(x=>x.dataset.endpoint);
  cfg.predocking = structuredClone(cfg.research.predocking);
}


// Values are rendered with textContent.
export function screeningDetail(panel, detail, el) {
  panel.replaceChildren();
  const r = detail.record;
  panel.append(el('h3', 'Candidate screening'), el('p', detail.note));
  const list = (title, pairs) => {
    const box = el('section'); box.append(el('h4', title));
    for (const [name, value] of pairs) box.append(el('p', name+': '+(value ?? 'Not evaluated')));
    panel.append(box);
  };
  list('Decision and identity', [['Molecule', r.molecule_id || r.compound_id || r.id], ['Final decision', r.decision], ['Structure filter', r.filter_status], ['Reason', r.reason], ['Original structure', r.original_smiles], ['Standardized structure', r.canonical_smiles], ['Rank', r.ranking?.rank]]);
  const table = el('table'), head = el('tr');
  for (const name of ['Endpoint', 'Prediction', 'Unit', 'Review threshold', 'Status']) head.append(el('th',name));
  table.append(head);
  for (const p of r.toxicity_predictions || []) {
    const tr = el('tr');
    for (const value of [p.endpoint, p.raw_prediction ?? '—', p.prediction_unit, p.threshold, p.prediction_status+' / '+p.risk_level]) tr.append(el('td', String(value ?? '—')));
    table.append(tr);
  }
  panel.append(el('h4','AI toxicity predictions'),table);
  list('Synthesis assessment', [['SA score (1–10)',r.synthesis?.sa_score], ['SCScore (1–5)',r.synthesis?.sc_score], ['Synthesis assessment',r.synthesis?.synthesis_difficulty], ['Route assessment',r.availability?.status], ['Evidence labels',(r.status_labels || []).join('; ')]]);
  const full = el('details'); full.append(el('summary','All alerts, model provenance, routes and planning building blocks'));
  const pre = el('pre',JSON.stringify(r,null,2));pre.style.cssText='white-space:pre-wrap;overflow-wrap:anywhere';full.append(pre);panel.append(full);
}
