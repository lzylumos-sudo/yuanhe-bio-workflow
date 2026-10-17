/* api.js -- thin fetch layer over the stdlib JSON API, with an in-memory cache
   that can be dropped when the pipeline has produced new files. */

const cache = new Map();

export async function get(path, { force = false } = {}) {
  if (!force && cache.has(path)) return cache.get(path);
  const res = await fetch(path, { headers: { Accept: 'application/json' } });
  const text = await res.text();
  let data;
  try {
    data = text ? JSON.parse(text) : null;
  } catch (err) {
    throw new Error(`Endpoint returned non- JSON content (${path}): ${text.slice(0, 160)}`);
  }
  if (!res.ok) {
    const detail = data && data.error ? data.error : `${res.status} ${res.statusText}`;
    const error = new Error(detail);
    error.status = res.status;
    error.payload = data;
    throw error;
  }
  cache.set(path, data);
  return data;
}

export async function post(path, body) {
  const res = await fetch(path, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json', Accept: 'application/json' },
    body: JSON.stringify(body || {}),
  });
  const text = await res.text();
  const data = text ? JSON.parse(text) : null;
  if (!res.ok) {
    const error = new Error((data && data.error) || `${res.status} ${res.statusText}`);
    error.status = res.status;
    error.payload = data;
    throw error;
  }
  return data;
}

export function invalidate(prefix = '') {
  for (const key of [...cache.keys()]) {
    if (!prefix || key.startsWith(prefix)) cache.delete(key);
  }
}

export function qs(params) {
  const parts = [];
  for (const [key, value] of Object.entries(params || {})) {
    if (value === undefined || value === null || value === '') continue;
    parts.push(`${encodeURIComponent(key)}=${encodeURIComponent(value)}`);
  }
  return parts.length ? `?${parts.join('&')}` : '';
}

/* ---- endpoint helpers ------------------------------------------------ */
export const api = {
  meta: () => get('/api/meta'),
  overview: () => get('/api/overview'),
  hotspots: () => get('/api/hotspots'),
  pharmacophore: () => get('/api/pharmacophore'),
  pocket: () => get('/api/pocket'),
  pocketPoints: (step = 0.8) => get(`/api/pocket/points${qs({ step })}`),
  axes: () => get('/api/axes'),
  scaffolds: () => get('/api/scaffolds'),
  filters: () => get('/api/filters'),
  poses: () => get('/api/poses'),
  scene: (params) => get(`/api/scene${qs(params)}`),
  evidenceReport: () => get('/api/evidence/report'),
  evidencePriors: () => get('/api/evidence/priors'),
  evidenceCensus: () => get('/api/evidence/census'),
  evidenceInterface: () => get('/api/evidence/interface'),
  docs: () => get('/api/docs'),
  doc: (name) => get(`/api/docs/${encodeURIComponent(name)}`),
  jobs: () => get('/api/jobs'),
  job: (id, since = 0) => get(`/api/jobs/${id}?since=${since}`),
  startJob: (kind, params) => post('/api/jobs', { kind, params }),
  cancelJob: (id) => post(`/api/jobs/${id}/cancel`, {}),
  health: () => get('/api/health'),
};
