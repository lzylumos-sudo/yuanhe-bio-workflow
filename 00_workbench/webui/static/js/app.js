/* app.js -- shell: hash router, navigation, theme, global status. */

import { api, invalidate } from './api.js';
import { el, div, span, toast } from './ui.js';

const NAV = [
  { group: 'Design and results', items: [
    { id: 'overview', label: 'Overview', file: 'overview.js', icon: 'M4 13h5v7H4zM10 4h4v16h-4zM15 9h5v11h-5z' },
    { id: 'scaffolds', label: 'Scaffold library', file: 'scaffolds.js', icon: 'M5 6h14M5 12h14M5 18h9' },
    { id: 'hotspots', label: 'Interface hotspots', file: 'hotspots.js', icon: 'M12 3v18M3 12h18' },
    { id: 'pharmacophore', label: 'Pharmacophore', file: 'pharmacophore.js', icon: 'M6 18a3 3 0 1 1 0-6 3 3 0 0 1 0 6zM18 9a3 3 0 1 1 0-6 3 3 0 0 1 0 6zM12 17a3 3 0 1 1 0-6 3 3 0 0 1 0 6z' },
    { id: 'pocket', label: 'Pocket', file: 'pocket.js', icon: 'M4 8c3-4 13-4 16 0 3 4-2 10-8 10S1 12 4 8z' },
  ] },
  { group: 'Structures and evidence', items: [
    { id: 'structure', label: '3D Structure', file: 'structure.js', icon: 'M12 3l8 4.5v9L12 21 4 16.5v-9z' },
    { id: 'evidence', label: 'Evidence chain', file: 'evidence.js', icon: 'M5 4h11l3 3v13H5zM8 11h8M8 15h5' },
  ] },
  { group: 'Execution and documentation', items: [
    { id: 'run', label: 'Run pipeline', file: 'run.js', icon: 'M7 5l12 7-12 7z' },
    { id: 'docs', label: 'Documentation', file: 'docs.js', icon: 'M6 3h9l4 4v14H6zM9 9h6M9 13h6M9 17h4' },
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
  document.title = `${current.label} · WEDGE-Design Console`;
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
  view.replaceChildren(div({ class: 'loading-pane' }, span({ class: 'spinner' }), `Loading「${current.label}」…`));

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
        if (!options.silent) toast('Result files reloaded', 'ok', 1800);
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
        el('h2', { text: 'Page load failed', style: { fontSize: '15px', marginBottom: '8px', color: 'var(--danger)' } }),
        div({ class: 'mono-sm', text: String(err && err.message ? err.message : err) }),
        div({ class: 'card-note', style: { marginTop: '8px' }, text: 'Check server logs; for endpoint errors, see the troubleshooting section in README for common causes.' }),
        div({ style: { marginTop: '12px' } }, el('button', { class: 'btn', text: 'Retry', onclick: () => mountPage() })),
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
    pill.textContent = ok ? 'Data ready' : 'Incomplete data';
    pill.className = `pill ${ok ? 'ok' : 'warn'}`;
    pill.title = 'results/ and evidence/ by run_pipeline.py / run_evidence.py Generate';
    foot.textContent = meta.available.filters ? `538→32 Structure filters ready` : 'Structure filters not run';
  } catch (err) {
    pill.textContent = 'Endpoint error';
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
  let text = 'Ready';
  if (running && running.length) {
    dot.className = 'dot dot-busy';
    text = `Running:${running[0].label}`;
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
    toast('Result files reloaded', 'ok', 1600);
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
