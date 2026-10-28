// Read-only until a validated submission succeeds. Never edit frozen policies.
export async function openEvidenceRevision({api, runId, stage, reload, el}) {
  const context = await api(`/api/project/run/${runId}/predocking/${stage}/revision`);
  const dialog = el('dialog'); dialog.id = 'evidenceRevisionDialog'; dialog.className = 'tool-dialog';
  const form = el('form'), title = el('h2', 'Revise pre-docking evidence'); form.className = 'tool-body';
  const error = el('p'); error.id = 'evidenceRevisionError'; error.setAttribute('role', 'alert');
  const input = (key, text, multiline = false, parent = form) => {
    const label = el('label', text), control = el(multiline ? 'textarea' : 'input');
    control.name = key; control.required = true;
    if (multiline) {control.rows = 12; control.style.width = '100%';}
    label.append(control); parent.append(label); return control;
  };
  form.append(title, el('p', context.policy), el('p',
    'Opening or cancelling leaves results unchanged. Save preserves old results and marks the affected stages for rerun; it does not qualify any molecule.'));
  const reviewer = input('reviewer', 'Reviewer');
  const reason = input('reason', 'Reason for this revision');
  const references = el('details'), summary = el('summary', 'Current alerts, routes and existing review records');
  references.append(summary, el('pre', JSON.stringify({candidates: context.candidates, existing_reviews: context.existing_reviews}, null, 2)));
  form.append(references);
  const field = stage === 'risk' ? 'risk_reviews' : 'route_reviews';
  const editor = el('details'); editor.append(el('summary', 'Add an individual evidence review'));
  const fields = el('div'); fields.className = 'form-grid'; editor.append(fields); form.append(editor);
  function select(name, text) {
    const label = el('label', text), control = el('select'); control.name = name;
    label.append(control); fields.append(label); return control;
  }
  const candidate = select('candidate', 'Candidate'); candidate.append(new Option('Select a candidate', ''));
  for (const item of context.candidates) candidate.append(new Option(item.compound_id, item.compound_id));
  const target = select('review_target', stage === 'risk' ? 'Structural alert' : 'Synthesis route');
  const targetStatus = el('p'); editor.append(targetStatus);
  candidate.onchange = () => {
    target.replaceChildren(new Option('Select evidence to review', ''));
    const item = context.candidates.find(c => c.compound_id === candidate.value);
    const evidence = stage === 'risk' ? item?.alerts : item?.routes;
    for (const [index, record] of (evidence || []).entries()) {
      const permitted = stage === 'risk' ? record.action === 'review' : record.reviewable;
      const option = new Option(stage === 'risk' ? `${record.rule_id} (${record.action})`
        : `Route ${index + 1}: ${record.problems.join(', ') || record.status}`, stage === 'risk' ? record.rule_id : record.route_sha256);
      option.disabled = !permitted; target.append(option);
    }
    targetStatus.textContent = item?.reason || '';
  };
  candidate.onchange();
  const date = input('reviewed_at', 'Review date', false, fields); date.type = 'date'; date.required = false;
  const notes = input('review_notes', 'Evidence and chemical reasoning', true, fields); notes.rows = 3; notes.required = false;
  const confirmed = el('input'); confirmed.type = 'checkbox'; confirmed.name = 'chemistry_confirmed';
  if (stage === 'availability') {
    const label = el('label', 'I have checked the target identity, stereochemistry and reaction selectivity.');
    label.prepend(confirmed); editor.append(label);
  }
  const advanced = el('details'); advanced.append(el('summary', 'Advanced review records (JSON)'));
  const reviews = input(field, 'New evidence review records (JSON array)', true, advanced); reviews.value = '[]'; reviews.rows = 5;
  const queued = el('p', 'No new reviews added. Saving can rerun unchanged evidence after a UI repair.');
  const add = el('button', 'Add review record'); add.type = 'button'; editor.append(add); form.append(queued, advanced);
  add.onclick = () => {
    error.textContent = '';
    try {
      if (!candidate.value || !target.value || !reviewer.value.trim() || !date.value || !notes.value.trim())
        throw Error('Select a candidate and reviewable evidence, then enter reviewer, date and chemical reasoning.');
      if (stage === 'availability' && !confirmed.checked) throw Error('Chemistry confirmation is required for a route review.');
      const records = JSON.parse(reviews.value);
      if (!Array.isArray(records)) throw Error('Evidence records must be a JSON array');
      records.push({compound_id: candidate.value, reviewer: reviewer.value, reason: notes.value, reviewed_at: date.value,
        ...(stage === 'risk' ? {rule_id: target.value} : {route_sha256: target.value, identity_stereo_selectivity_confirmed: true})});
      reviews.value = JSON.stringify(records, null, 2); queued.textContent = `${records.length} review record(s) ready to submit.`;
      target.value = ''; notes.value = ''; confirmed.checked = false;
    } catch (e) {error.textContent = e.message;}
  };
  form.append(el('p', 'Hard exclusions and incomplete synthesis routes remain blocked. Each review applies only to its recorded molecule and evidence.'));
  const save = el('button', 'Save evidence revision'), cancel = el('button', 'Cancel');
  save.type = 'submit'; cancel.type = 'button';
  cancel.onclick = () => dialog.close();
  dialog.onclose = () => dialog.remove();
  form.append(error, save, cancel); dialog.append(form); document.body.append(dialog); dialog.showModal();
  form.onsubmit = async event => {
    event.preventDefault(); error.textContent = '';
    try {
      const records = JSON.parse(reviews.value);
      if (!Array.isArray(records)) throw Error('Evidence records must be a JSON array');
      save.disabled = true; cancel.disabled = true;
      await api(`/api/project/run/${runId}/predocking/${stage}/revise`, {
        revision_token: context.revision_token, reviewer: reviewer.value, reason: reason.value,
        risk_reviews: [], route_reviews: [], [field]: records,
      });
      dialog.close(); await reload(runId);
    } catch (e) {error.textContent = e.message;}
    finally {save.disabled = false; cancel.disabled = false;}
  };
  dialog.oncancel = event => {if (save.disabled) event.preventDefault();};
}
