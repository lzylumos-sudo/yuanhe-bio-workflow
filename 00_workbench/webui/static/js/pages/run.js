/* pages/run.js -- launch the pipeline entry points and watch stdout live. */

import { api, invalidate } from '../api.js';
import {
  el, div, span, card, chip, kvList, fmt, toast, downloadText, formatNow,
} from '../ui.js';

const STATUS_LABEL = {
  queued: 'Queued', running: 'Running', done: 'Completed', failed: 'Failed', cancelled: 'Cancelled',
};

export async function mount(root, ctx) {
  const data = await api.jobs();
  const kinds = data.kinds || [];
  let jobs = data.jobs || [];
  let selected = ctx.params.job || (jobs[0] && jobs[0].id) || null;
  let cursor = 0;
  let logLines = [];
  let timer = null;
  let stopped = false;

  root.replaceChildren();
  root.appendChild(div({ class: 'page-head' },
    div(null,
      el('h1', { text: 'Run pipeline' }),
      div({ class: 'sub', text: 'Rerun the four bundled entry scripts. Commands are allowlisted, arguments validated and subprocess logs streamed. Results are written to results/ and evidence/, ; refresh other pages for updated values.' })),
    div({ class: 'actions' },
      el('button', { class: 'btn', text: 'Refresh status', onclick: () => reloadJobs() }))));

  /* ---------- launcher cards ---------- */
  const launchCards = div({ class: 'grid cols-2', style: { gap: '12px' } });
  for (const kind of kinds) {
    const inputs = {};
    const fields = [];
    for (const [name, spec] of Object.entries(kind.args || {})) {
      const input = el('input', {
        type: spec.type === 'path' ? 'text' : 'number',
        value: String(spec.default ?? ''),
        step: spec.type === 'float' ? '0.5' : '1',
        style: { width: '110px' },
      });
      inputs[name] = input;
      fields.push(el('label', { class: 'check' }, span({ class: 'tiny muted', text: name }), input));
    }
    const btn = el('button', {
      class: kind.kind === 'pipeline' ? 'btn primary' : 'btn',
      text: 'Run',
      onclick: async () => {
        btn.disabled = true;
        try {
          const params = {};
          for (const [name, input] of Object.entries(inputs)) params[name] = input.value;
          const res = await api.startJob(kind.kind, params);
          toast(`Started:${kind.label}`, 'ok');
          selected = res.id;
          cursor = 0;
          logLines = [];
          invalidate('/api/jobs');
          await reloadJobs();
          paintConsole();
          startPolling();
        } catch (err) {
          toast(err.message, 'err', 5000);
        } finally {
          btn.disabled = false;
        }
      },
    });
    launchCards.appendChild(card({ class: 'tight' },
      div({ class: 'split' },
        span({ class: 'tiny muted mono', text: kind.script }),
        span({ class: 'right tiny muted', style: { marginLeft: 'auto' }, text: kind.kind })),
      div({ style: { fontWeight: 600, margin: '2px 0 4px' }, text: kind.label }),
      div({ class: 'card-note', text: kind.docs }),
      div({ class: 'split', style: { marginTop: '8px' } }, ...fields,
        span({ style: { marginLeft: 'auto' } }, btn))));
  }
  root.appendChild(launchCards);

  /* ---------- history + console ---------- */
  const history = div({ style: { display: 'grid', gap: '6px' } });
  const consoleNode = div({ class: 'console' });
  const consoleTitle = div({ class: 'card-head' });
  const consoleTools = div({ class: 'split' });
  const consoleCard = card({}, consoleTitle, consoleTools, consoleNode);

  function paintHistory() {
    history.replaceChildren();
    if (!jobs.length) {
      history.appendChild(div({ class: 'empty-state', text: 'No run records yet.' }));
      return;
    }
    for (const job of jobs) {
      const running = job.status === 'running' || job.status === 'queued';
      history.appendChild(div({
        style: {
          padding: '7px 9px', border: '1px solid', borderRadius: '8px', cursor: 'pointer',
          borderColor: job.id === selected ? 'rgba(56,189,248,.4)' : 'var(--line)',
          background: job.id === selected ? 'rgba(56,189,248,.08)' : 'transparent',
        },
        onclick: () => { selected = job.id; cursor = 0; logLines = []; paintHistory(); paintConsole(); startPolling(); },
      },
        div({ class: 'split' },
          span({ class: `status-chip ${job.status}`, text: STATUS_LABEL[job.status] || job.status }),
          span({ class: 'tiny', text: job.label }),
          span({ class: 'right tiny muted mono', style: { marginLeft: 'auto' }, text: job.id })),
        div({ class: 'tiny muted mono', style: { marginTop: '2px' }, text: `exit ${job.exit ?? '—'} · ${job.total_lines} rows` })));
    }
  }

  function paintConsole() {
    const job = jobs.find((j) => j.id === selected);
    consoleTitle.replaceChildren(
      el('h3', { text: job ? `${job.label} · ${job.id}` : 'Log' }),
      job ? span({ class: `right status-chip ${job.status}`, text: STATUS_LABEL[job.status] || job.status }) : null);
    consoleTools.replaceChildren(
      job ? span({ class: 'tiny muted mono', text: (job.cmd || []).slice(1).join(' ') || job.kind }) : span({ class: 'tiny muted', text: 'Select a run on the left' }),
      div({ class: 'grow' }),
      job && (job.status === 'running' || job.status === 'queued')
        ? el('button', { class: 'btn sm danger', text: 'Stop', onclick: () => cancel(job.id) })
        : null,
      el('button', {
        class: 'btn sm ghost', text: 'Export log',
        onclick: () => downloadText(`job-${job?.id || 'log'}-${formatNow()}.txt`, logLines.join('\n'), 'text/plain;charset=utf-8'),
      }),
      el('button', { class: 'btn sm ghost', text: 'Clear', onclick: () => { logLines = []; consoleNode.replaceChildren(); } }));
    consoleNode.replaceChildren(...logLines.map((line) => {
      const cls = /\[Completed\]/.test(line.text) ? 'done' : /\[Failed\]|Traceback|Error/.test(line.text) ? 'fail' : /^\$/.test(line.text) ? 'cmd' : '';
      return div({ class: `ln ${cls}` },
        line.t ? span({ class: 't', text: new Date(line.t * 1000).toLocaleTimeString('en') }) : null,
        line.text);
    }));
    consoleNode.scrollTop = consoleNode.scrollHeight;
  }

  async function cancel(id) {
    try {
      await api.cancelJob(id);
      toast('Stop requested', 'warn');
    } catch (err) {
      toast(err.message, 'err');
    }
  }

  async function reloadJobs() {
    try {
      invalidate('/api/jobs');
      const fresh = await api.jobs();
      jobs = fresh.jobs || [];
      if (!selected) selected = jobs[0]?.id || null;
      paintHistory();
      paintConsole();
    } catch (err) {
      toast(`Could not read job status:${err.message}`, 'err');
    }
  }

  async function pollOnce() {
    if (!selected) return;
    try {
      const job = await api.job(selected, cursor);
      if (job.lines?.length) {
        logLines = logLines.concat(job.lines);
        cursor = job.cursor;
        paintConsole();
      }
      const idx = jobs.findIndex((j) => j.id === job.id);
      if (idx >= 0) {
        const wasRunning = jobs[idx].status === 'running' || jobs[idx].status === 'queued';
        jobs[idx] = { ...jobs[idx], ...job, lines: undefined };
        if (wasRunning && job.status !== 'running' && job.status !== 'queued') {
          clearInterval(timer);
          paintHistory();
          paintConsole();
          if (job.status === 'done') {
            invalidate();
            toast(`${job.label} Completed; result files updated. Refresh other pages for new values.`, 'ok', 6000);
          } else {
            toast(`${job.label} ${STATUS_LABEL[job.status] || job.status}`, 'err', 5000);
          }
          return;
        }
        paintHistory();
      }
    } catch (err) {
      /* transient */
    }
  }

  function startPolling() {
    clearInterval(timer);
    timer = setInterval(() => {
      if (stopped) return;
      pollOnce();
    }, 800);
  }

  paintHistory();
  paintConsole();
  startPolling();
  reloadJobs();

  root.appendChild(div({ class: 'section grid wide-left' },
    card({ title: 'Run history' }, history),
    consoleCard));

  return () => {
    stopped = true;
    clearInterval(timer);
  };
}
