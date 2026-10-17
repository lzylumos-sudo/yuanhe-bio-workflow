/* pages/run.js -- launch the pipeline entry points and watch stdout live. */

import { api, invalidate } from '../api.js';
import {
  el, div, span, card, chip, kvList, fmt, toast, downloadText, formatNow,
} from '../ui.js';

const STATUS_LABEL = {
  queued: '排队中', running: '运行中', done: '完成', failed: '失败', cancelled: '已取消',
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
      el('h1', { text: '运行管线' }),
      div({ class: 'sub', text: '在这里重跑仓库自带的四个入口脚本。命令被白名单限定，参数经过校验；日志从子进程逐行流式返回。结果写回 results/ 与 evidence/，其它页面刷新后即可看到新数字。' })),
    div({ class: 'actions' },
      el('button', { class: 'btn', text: '刷新状态', onclick: () => reloadJobs() }))));

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
      text: '运行',
      onclick: async () => {
        btn.disabled = true;
        try {
          const params = {};
          for (const [name, input] of Object.entries(inputs)) params[name] = input.value;
          const res = await api.startJob(kind.kind, params);
          toast(`已启动：${kind.label}`, 'ok');
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
      history.appendChild(div({ class: 'empty-state', text: '还没有运行记录。' }));
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
        div({ class: 'tiny muted mono', style: { marginTop: '2px' }, text: `exit ${job.exit ?? '—'} · ${job.total_lines} 行` })));
    }
  }

  function paintConsole() {
    const job = jobs.find((j) => j.id === selected);
    consoleTitle.replaceChildren(
      el('h3', { text: job ? `${job.label} · ${job.id}` : '日志' }),
      job ? span({ class: `right status-chip ${job.status}`, text: STATUS_LABEL[job.status] || job.status }) : null);
    consoleTools.replaceChildren(
      job ? span({ class: 'tiny muted mono', text: (job.cmd || []).slice(1).join(' ') || job.kind }) : span({ class: 'tiny muted', text: '选择左侧任意一次运行' }),
      div({ class: 'grow' }),
      job && (job.status === 'running' || job.status === 'queued')
        ? el('button', { class: 'btn sm danger', text: '中止', onclick: () => cancel(job.id) })
        : null,
      el('button', {
        class: 'btn sm ghost', text: '导出日志',
        onclick: () => downloadText(`job-${job?.id || 'log'}-${formatNow()}.txt`, logLines.join('\n'), 'text/plain;charset=utf-8'),
      }),
      el('button', { class: 'btn sm ghost', text: '清屏', onclick: () => { logLines = []; consoleNode.replaceChildren(); } }));
    consoleNode.replaceChildren(...logLines.map((line) => {
      const cls = /\[完成\]/.test(line.text) ? 'done' : /\[失败\]|Traceback|Error/.test(line.text) ? 'fail' : /^\$/.test(line.text) ? 'cmd' : '';
      return div({ class: `ln ${cls}` },
        line.t ? span({ class: 't', text: new Date(line.t * 1000).toLocaleTimeString('zh-CN') }) : null,
        line.text);
    }));
    consoleNode.scrollTop = consoleNode.scrollHeight;
  }

  async function cancel(id) {
    try {
      await api.cancelJob(id);
      toast('已请求中止', 'warn');
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
      toast(`读取任务状态失败：${err.message}`, 'err');
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
            toast(`${job.label} 完成，结果文件已更新（其它页面刷新即可看到新数字）`, 'ok', 6000);
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
    card({ title: '运行历史' }, history),
    consoleCard));

  return () => {
    stopped = true;
    clearInterval(timer);
  };
}
