/* app.js -- shell: hash router, navigation, theme, global status. */

import { api, invalidate } from './api.js';
import { el, div, span, toast } from './ui.js';

const NAV = [
  { group: '设计与结果', items: [
    { id: 'overview', label: '概览', file: 'overview.js', icon: 'M4 13h5v7H4zM10 4h4v16h-4zM15 9h5v11h-5z' },
    { id: 'scaffolds', label: '骨架库', file: 'scaffolds.js', icon: 'M5 6h14M5 12h14M5 18h9' },
    { id: 'hotspots', label: '界面热点', file: 'hotspots.js', icon: 'M12 3v18M3 12h18' },
    { id: 'pharmacophore', label: '药效团', file: 'pharmacophore.js', icon: 'M6 18a3 3 0 1 1 0-6 3 3 0 0 1 0 6zM18 9a3 3 0 1 1 0-6 3 3 0 0 1 0 6zM12 17a3 3 0 1 1 0-6 3 3 0 0 1 0 6z' },
    { id: 'pocket', label: '口袋', file: 'pocket.js', icon: 'M4 8c3-4 13-4 16 0 3 4-2 10-8 10S1 12 4 8z' },
  ] },
  { group: '结构与证据', items: [
    { id: 'structure', label: '3D 结构', file: 'structure.js', icon: 'M12 3l8 4.5v9L12 21 4 16.5v-9z' },
    { id: 'evidence', label: '证据链', file: 'evidence.js', icon: 'M5 4h11l3 3v13H5zM8 11h8M8 15h5' },
  ] },
  { group: '运行与文档', items: [
    { id: 'run', label: '运行管线', file: 'run.js', icon: 'M7 5l12 7-12 7z' },
    { id: 'docs', label: '文档', file: 'docs.js', icon: 'M6 3h9l4 4v14H6zM9 9h6M9 13h6M9 17h4' },
  ] },
];

const state = {
  meta: null,
  page: null,
  pageModule: null,
  route: 'overview',
  params: {},
  jobs: [],
  theme: localStorage.getItem('wedge-theme') || 'dark',
};

const view = document.getElementById('view');
const navEl = document.getElementById('nav');
const crumbEl = document.getElementById('crumb');

/* ---------------- theme ---------------- */
function applyTheme() {
  document.documentElement.dataset.theme = state.theme;
  localStorage.setItem('wedge-theme', state.theme);
}

/* ---------------- navigation ---------------- */
function buildNav() {
  navEl.replaceChildren();
  for (const group of NAV) {
    navEl.appendChild(div({ class: 'nav-group', text: group.group }));
    for (const item of group.items) {
      const a = el('a', { class: 'nav-item', href: `#/${item.id}`, dataset: { id: item.id } },
        icon(item.icon), span({ text: item.label }));
      navEl.appendChild(a);
    }
  }
}

function icon(path, extra = '') {
  const node = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  node.setAttribute('viewBox', '0 0 24 24');
  node.setAttribute('fill', 'none');
  node.innerHTML = `<path d="${path}" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round"/>${extra}`;
  return node;
}

function parseHash() {
  const raw = (location.hash || '#/overview').replace(/^#\/?/, '');
  const [route, query] = raw.split('?');
  const params = {};
  for (const [key, value] of new URLSearchParams(query || '')) params[key] = value;
  return { route: route || 'overview', params };
}

function navigate(route, params = {}) {
  const query = new URLSearchParams(params).toString();
  location.hash = `#/${route}${query ? `?${query}` : ''}`;
}

function highlightNav(route) {
  const known = NAV.flatMap((g) => g.items);
  const current = known.find((item) => item.id === route) || known[0];
  for (const node of navEl.querySelectorAll('.nav-item')) {
    node.classList.toggle('active', node.dataset.id === current.id);
  }
  crumbEl.replaceChildren(
    span({ class: 'crumb-root', text: 'WEDGE-Design' }),
    span({ class: 'sep', text: '/' }),
    el('b', { text: current.label }),
  );
  document.title = `${current.label} · WEDGE-Design 控制台`;
  return current;
}

/* ---------------- page mounting ---------------- */
async function mountPage() {
  const { route, params } = parseHash();
  state.route = route;
  state.params = params;
  const current = highlightNav(route);

  if (state.pageModule?.unmount) {
    try { state.pageModule.unmount(); } catch (err) { console.warn(err); }
  }
  state.pageModule = null;
  document.getElementById('sidebar')?.classList.remove('open');
  view.replaceChildren(div({ class: 'loading-pane' }, span({ class: 'spinner' }), `正在载入「${current.label}」…`));

  try {
    const module = await import(`./pages/${current.file}`);
    const ctx = {
      params,
      navigate,
      toast,
      meta: state.meta,
      refresh: async (options = {}) => {
        invalidate();
        state.meta = await api.meta().catch(() => state.meta);
        paintStatus();
        await mountPage();
        if (!options.silent) toast('已重新读取结果文件', 'ok', 1800);
      },
    };
    view.replaceChildren();
    const cleanup = await module.mount(view, ctx);
    state.pageModule = { ...module, unmount: cleanup };
    view.scrollTop = 0;
  } catch (err) {
    console.error(err);
    view.replaceChildren(
      div({ class: 'card' },
        el('h2', { text: '页面载入失败', style: { fontSize: '15px', marginBottom: '8px', color: 'var(--danger)' } }),
        div({ class: 'mono-sm', text: String(err && err.message ? err.message : err) }),
        div({ class: 'card-note', style: { marginTop: '8px' }, text: '检查终端里的服务端日志；如果是接口错误，README 中的“排错”一节列出了常见原因。' }),
        div({ style: { marginTop: '12px' } }, el('button', { class: 'btn', text: '重试', onclick: () => mountPage() })),
      ),
    );
  }
  updateFootStatus();
}

/* ---------------- global status ---------------- */
async function paintStatus() {
  const pill = document.getElementById('pill-data');
  const foot = document.getElementById('foot-meta');
  try {
    state.meta = state.meta || await api.meta();
    const meta = state.meta;
    const ok = meta.available.results && meta.available.evidence;
    pill.textContent = ok ? '数据就绪' : '数据不完整';
    pill.className = `pill ${ok ? 'ok' : 'warn'}`;
    pill.title = 'results/ 与 evidence/ 由 run_pipeline.py / run_evidence.py 生成';
    foot.textContent = meta.available.filters ? `538→32 结构过滤已就绪` : '结构过滤未运行';
  } catch (err) {
    pill.textContent = '接口异常';
    pill.className = 'pill warn';
    foot.textContent = String(err.message || err);
  }
}

async function pollJobs() {
  try {
    invalidate('/api/jobs');
    const data = await api.jobs();
    state.jobs = data.jobs || [];
    const running = state.jobs.filter((j) => j.status === 'running' || j.status === 'queued');
    updateFootStatus(running);
  } catch (err) {
    /* status polling is best effort */
  }
}

function updateFootStatus(runningOverride) {
  const node = document.getElementById('foot-status');
  if (!node) return;
  const running = runningOverride || state.jobs.filter((j) => j.status === 'running' || j.status === 'queued');
  const dot = span({ class: 'dot' });
  let text = '就绪';
  if (running && running.length) {
    dot.className = 'dot dot-busy';
    text = `运行中：${running[0].label}`;
  } else {
    dot.className = 'dot dot-ok';
  }
  node.replaceChildren(dot, span({ text }));
}

/* ---------------- boot ---------------- */
function onHashChange() { mountPage(); }

async function boot() {
  applyTheme();
  buildNav();
  document.getElementById('btn-theme').addEventListener('click', () => {
    state.theme = state.theme === 'dark' ? 'light' : 'dark';
    applyTheme();
    state.pageModule?.onThemeChange?.();
  });
  document.getElementById('btn-refresh').addEventListener('click', () => {
    invalidate();
    mountPage();
    toast('已重新读取结果文件', 'ok', 1600);
  });
  document.getElementById('nav-toggle').addEventListener('click', () => {
    document.getElementById('sidebar').classList.toggle('open');
  });
  window.addEventListener('hashchange', onHashChange);
  window.addEventListener('resize', () => state.pageModule?.onResize?.());

  await paintStatus();
  await mountPage();
  pollJobs();
  setInterval(pollJobs, 4000);
}

boot();
