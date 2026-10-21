/* flow/main.js -- hash router + stepper for the 7-step workflow UI. */
import { get } from '../api.js';
import { el } from '../ui.js';
import { mountStage } from './stage.js';

const page = document.getElementById('page');
const stepper = document.getElementById('stepper');
const toggle = document.getElementById('demo-toggle');

const state = { demo: localStorage.getItem('flow.demo') === '1', wf: null, cleanup: null, token: 0 };
toggle.checked = state.demo;

const stageNo = () => {
  const m = location.hash.match(/#\/(\d)/);
  const n = m ? Number(m[1]) : 1;
  return Math.min(7, Math.max(1, n));
};

export const go = (n) => { location.hash = `#/${n}`; };

async function loadWorkflow() {
  state.wf = await get(`/api/workflow?demo=${state.demo ? 1 : 0}`);
  return state.wf;
}

function drawStepper(active) {
  stepper.replaceChildren();
  for (const s of state.wf.stages) {
    const reached = s.n <= active;
    const done = s.status !== 'pending' && s.n < active;
    const cnt = s.out === null || s.out === undefined ? '—' : s.out.toLocaleString('en-US');
    const b = el('button', {
      class: ['step', reached && 'reached', done && 'done', s.n === active && 'active'],
      title: s.subtitle,
      'aria-current': s.n === active ? 'step' : null,
      onclick: () => go(s.n),
    },
    el('span', { class: 'dot', text: s.n }),
    el('span', { class: 'name', text: s.short }),
    el('span', { class: 'cnt', text: `${cnt} ${s.unit}` }));
    stepper.appendChild(b);
  }
}

async function render() {
  const token = ++state.token;
  const n = stageNo();
  if (state.cleanup) { try { state.cleanup(); } catch (e) { /* ignore */ } state.cleanup = null; }
  page.replaceChildren(el('div', { class: 'loading', text: 'Loading...' }));
  try {
    await loadWorkflow();
    if (token !== state.token) return;
    drawStepper(n);
    page.replaceChildren();
    const stage = state.wf.stages[n - 1];
    document.title = `${stage.n}/7 ${stage.title} · Yuanhe Bio`;
    const cleanup = await mountStage(page, { stage, wf: state.wf, demo: state.demo, go });
    if (token !== state.token) { cleanup?.(); return; }
    state.cleanup = cleanup;
    window.scrollTo(0, 0);
  } catch (err) {
    page.replaceChildren(el('div', { class: 'err', text: `Load failed:${err.message}` }));
  }
}

toggle.addEventListener('change', () => {
  state.demo = toggle.checked;
  localStorage.setItem('flow.demo', state.demo ? '1' : '0');
  render();
});
window.addEventListener('hashchange', render);
if (!location.hash) location.hash = '#/1'; else render();

