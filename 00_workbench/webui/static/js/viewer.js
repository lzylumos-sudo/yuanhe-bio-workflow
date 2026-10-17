/* viewer.js -- a small interactive 3D structure viewer for the WEDGE-Design
   console, written from scratch on a Canvas 2D context (no WebGL, no external
   libraries, works offline).

   Geometry comes from /api/scene: C-alpha trace, pocket-lining side chains,
   the three phage peptides (superposed into one frame), the pharmacophore
   features, validated hotspot residues, and the designed ligand pose.

   Rendering: painter's algorithm over depth-sorted primitives; the buried
   cavity point cloud is drawn as colour-bucketed batches behind everything. */

const TAU = Math.PI * 2;

/* ---------------- small matrix helpers (row-major 3x3) ---------------- */
function matMul(a, b) {
  const out = new Array(9);
  for (let r = 0; r < 3; r += 1) {
    for (let c = 0; c < 3; c += 1) {
      out[r * 3 + c] = a[r * 3] * b[c] + a[r * 3 + 1] * b[3 + c] + a[r * 3 + 2] * b[6 + c];
    }
  }
  return out;
}
const matRotX = (t) => [1, 0, 0, 0, Math.cos(t), -Math.sin(t), 0, Math.sin(t), Math.cos(t)];
const matRotY = (t) => [Math.cos(t), 0, Math.sin(t), 0, 1, 0, -Math.sin(t), 0, Math.cos(t)];
const matRotZ = (t) => [Math.cos(t), -Math.sin(t), 0, Math.sin(t), Math.cos(t), 0, 0, 0, 1];

function normalize(v) {
  const n = Math.hypot(v[0], v[1], v[2]) || 1;
  return [v[0] / n, v[1] / n, v[2] / n];
}
function cross(a, b) {
  return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
}
function dot(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

function hexToRgb(hex) {
  const clean = hex.replace('#', '');
  const full = clean.length === 3 ? clean.split('').map((c) => c + c).join('') : clean;
  return [parseInt(full.slice(0, 2), 16), parseInt(full.slice(2, 4), 16), parseInt(full.slice(4, 6), 16)];
}
function shade(hex, factor) {
  const [r, g, b] = hexToRgb(hex);
  const f = (v) => Math.max(0, Math.min(255, Math.round(v * factor)));
  return `rgb(${f(r)},${f(g)},${f(b)})`;
}
function withAlpha(hex, alpha) {
  const [r, g, b] = hexToRgb(hex);
  return `rgba(${r},${g},${b},${alpha})`;
}

/* ---------------- viewer ---------------- */
export class StructureViewer {
  constructor(container, opts = {}) {
    this.container = container;
    this.opts = Object.assign({
      onHover: null, onPick: null,
      background: null,
    }, opts);

    this.canvas = document.createElement('canvas');
    this.ctx = this.canvas.getContext('2d');
    container.appendChild(this.canvas);

    this.scene = null;
    this.pocket = null;
    this.view = {
      rot: [1, 0, 0, 0, 1, 0, 0, 0, 1],
      zoom: 0.9,
      panX: 0,
      panY: 0,
      center: [0, 0, 0],
      baseScale: 1,
    };
    this.options = {
      showPocket: true,
      pocketMode: 'chem',          // 'chem' | 'sub'
      showTrace: true,
      showSidechains: true,
      showPeptides: true,
      peptideFilter: {},           // tag -> bool
      showLigand: true,
      showFeatures: true,
      showHotspots: true,
      showLabels: true,
      spin: false,
      ligandRadius: 0.55,
    };
    this.hover = null;
    this._pickables = [];
    this._pocketBuckets = [];
    this._depthRange = [0, 1];
    this._frame = null;
    this._spinTimer = null;
    this._bindEvents();
    this._observeSize();
    this._resize();
  }

  /* ---------- lifecycle ---------- */
  destroy() {
    this._resizeObserver?.disconnect();
    this._stopSpin();
    this.canvas.remove();
  }

  _observeSize() {
    this._resizeObserver = new ResizeObserver(() => { this._resize(); });
    this._resizeObserver.observe(this.container);
  }

  _resize() {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const w = this.container.clientWidth || 640;
    const h = this.container.clientHeight || 420;
    this.canvas.width = Math.round(w * dpr);
    this.canvas.height = Math.round(h * dpr);
    this.canvas.style.width = `${w}px`;
    this.canvas.style.height = `${h}px`;
    this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    this.width = w;
    this.height = h;
    this.requestRender();
  }

  /* ---------- data ---------- */
  setScene(scene) {
    this.scene = scene;
    const ca = scene.receptor?.ca || [];
    const ligand = scene.ligand;
    let sum = [0, 0, 0];
    for (const p of ca) { sum = [sum[0] + p[0], sum[1] + p[1], sum[2] + p[2]]; }
    const n = ca.length || 1;
    const proteinCenter = [sum[0] / n, sum[1] / n, sum[2] / n];
    this.proteinCenter = proteinCenter;

    // initial orientation: put the ligand (or cavity) between camera and protein
    const target = ligand ? meanOf(ligand.atoms.map((a) => a.p)) : (scene.center || proteinCenter);
    const view = normalize([target[0] - proteinCenter[0], target[1] - proteinCenter[1], target[2] - proteinCenter[2]]);
    let up = [0, 1, 0];
    if (Math.abs(dot(up, view)) > 0.95) up = [0, 0, 1];
    const xAxis = normalize(cross(up, view));
    const yAxis = cross(view, xAxis);
    this.view.rot = [xAxis[0], xAxis[1], xAxis[2], yAxis[0], yAxis[1], yAxis[2], view[0], view[1], view[2]];
    this.view.center = proteinCenter;
    this.view.panX = 0;
    this.view.panY = 0;
    this.view.zoom = 1;
    this._pocketBuckets = [];
    this.pocket = null;
    this.fit();
    this.requestRender();
  }

  setPocket(pocket) {
    this.pocket = pocket;
    this._buildPocketBuckets();
    this.requestRender();
  }

  setOptions(patch) {
    Object.assign(this.options, patch);
    if (patch.pocketMode && this.pocket) this._buildPocketBuckets();
    if (patch.spin !== undefined) {
      if (patch.spin) this._startSpin(); else this._stopSpin();
    }
    this.requestRender();
  }

  _buildPocketBuckets() {
    const pocket = this.pocket;
    if (!pocket || !pocket.xyz || !pocket.xyz.length) { this._pocketBuckets = []; return; }
    const chemColors = ['#4f8ef7', '#e55353', '#38a169', '#c99a3a', '#8b97a8'];
    const subColors = ['#4f8ef7', '#e55353', '#38a169', '#c99a3a'];
    const useSub = this.options.pocketMode === 'sub' && Array.isArray(pocket.sub) && pocket.sub[0]?.length;
    const bucketOf = (i) => (useSub ? (pocket.sub[0][i] % 4) : (pocket.chem[i] ?? 4));
    const palette = useSub ? subColors : chemColors;
    const buckets = palette.map(() => []);
    for (let i = 0; i < pocket.xyz.length / 3; i += 1) {
      buckets[bucketOf(i)].push(pocket.xyz[i * 3], pocket.xyz[i * 3 + 1], pocket.xyz[i * 3 + 2]);
    }
    this._pocketBuckets = buckets.map((coords, idx) => ({ color: palette[idx], coords: Float32Array.from(coords) }));
  }

  /* ---------- camera ---------- */
  fit() {
    const pts = [];
    const scene = this.scene;
    if (!scene) return;
    for (const p of scene.receptor?.ca || []) pts.push(p);
    for (const pep of scene.peptides || []) for (const p of pep.ca || []) pts.push(p);
    if (scene.ligand) for (const a of scene.ligand.atoms) pts.push(a.p);
    for (const f of scene.features || []) pts.push(f.p);
    if (!pts.length) return;
    let min = [Infinity, Infinity, Infinity];
    let max = [-Infinity, -Infinity, -Infinity];
    for (const p of pts) {
      for (let k = 0; k < 3; k += 1) { min[k] = Math.min(min[k], p[k]); max[k] = Math.max(max[k], p[k]); }
    }
    this.view.center = [(min[0] + max[0]) / 2, (min[1] + max[1]) / 2, (min[2] + max[2]) / 2];
    let radius = 1;
    for (const p of pts) {
      radius = Math.max(radius, Math.hypot(p[0] - this.view.center[0], p[1] - this.view.center[1], p[2] - this.view.center[2]));
    }
    this._radius = radius;
    this.view.zoom = 1;
    this.view.panX = 0;
    this.view.panY = 0;
    this.requestRender();
  }

  zoomBy(factor) {
    this.view.zoom = Math.max(0.25, Math.min(9, this.view.zoom * factor));
    this.requestRender();
  }

  /* ---------- projection ---------- */
  _scale() {
    const radius = this._radius || 40;
    return (Math.min(this.width, this.height) / (2.35 * radius)) * this.view.zoom;
  }

  _project(p) {
    const c = this.view.center;
    const r = this.view.rot;
    const x0 = p[0] - c[0];
    const y0 = p[1] - c[1];
    const z0 = p[2] - c[2];
    const x = r[0] * x0 + r[1] * y0 + r[2] * z0;
    const y = r[3] * x0 + r[4] * y0 + r[5] * z0;
    const z = r[6] * x0 + r[7] * y0 + r[8] * z0;
    const s = this._scale();
    return [this.width / 2 + this.view.panX + x * s, this.height / 2 + this.view.panY - y * s, z];
  }

  /* ---------- interaction ---------- */
  _bindEvents() {
    const canvas = this.canvas;
    let dragging = null;
    let last = { x: 0, y: 0 };
    let moved = 0;

    const onDown = (event) => {
      if (event.button !== 0 && event.button !== 1 && event.button !== 2) return;
      dragging = { pan: event.shiftKey || event.button === 1 || event.button === 2, button: event.button };
      last = { x: event.clientX, y: event.clientY };
      moved = 0;
      canvas.classList.add('dragging');
      canvas.setPointerCapture?.(event.pointerId);
      event.preventDefault();
    };
    const onMove = (event) => {
      if (!dragging) { this._hoverAt(event); return; }
      const dx = event.clientX - last.x;
      const dy = event.clientY - last.y;
      last = { x: event.clientX, y: event.clientY };
      moved += Math.abs(dx) + Math.abs(dy);
      if (dragging.pan) {
        this.view.panX += dx;
        this.view.panY += dy;
      } else {
        const k = 0.0075;
        const rot = matMul(matMul(matRotX(dy * k), matRotY(dx * k)), this.view.rot);
        this.view.rot = rot;
      }
      this.requestRender();
    };
    const onUp = (event) => {
      canvas.classList.remove('dragging');
      if (dragging && moved < 5) {
        const pick = this._pickAt(event);
        if (pick && this.opts.onPick) this.opts.onPick(pick);
      }
      dragging = null;
      canvas.releasePointerCapture?.(event.pointerId);
    };
    canvas.addEventListener('pointerdown', onDown);
    canvas.addEventListener('pointermove', onMove);
    canvas.addEventListener('pointerup', onUp);
    canvas.addEventListener('pointerleave', () => {
      dragging = null;
      this.hover = null;
      this.opts.onHover?.(null);
      this.requestRender();
    });
    canvas.addEventListener('wheel', (event) => {
      event.preventDefault();
      this.zoomBy(Math.exp(-event.deltaY * 0.0016));
    }, { passive: false });
    canvas.addEventListener('dblclick', () => this.fit());
    canvas.addEventListener('contextmenu', (event) => event.preventDefault());
  }

  _pickAt(event) {
    const rect = this.canvas.getBoundingClientRect();
    const x = event.clientX - rect.left;
    const y = event.clientY - rect.top;
    let best = null;
    let bestScore = Infinity;
    const radius = this._radius || 40;
    for (const item of this._pickables) {
      const d = Math.hypot(item.x - x, item.y - y);
      const limit = Math.max(6, item.r || 6);
      if (d > limit) continue;
      const score = d - (item.z / (2 * radius)) * 10;   // prefer the closest atom
      if (score < bestScore) { bestScore = score; best = item; }
    }
    return best;
  }

  _hoverAt(event) {
    const item = this._pickAt(event);
    const same = (item && this.hover && item.label === this.hover.label)
      || (!item && !this.hover);
    if (same) return;
    this.hover = item;
    this.opts.onHover?.(item, event);
    this.requestRender();
  }

  /* ---------- render ---------- */
  requestRender() {
    if (this._frame) return;
    this._frame = requestAnimationFrame(() => {
      this._frame = null;
      this.render();
    });
  }

  render() {
    const ctx = this.ctx;
    const W = this.width;
    const H = this.height;
    ctx.clearRect(0, 0, W, H);
    const gradient = ctx.createRadialGradient(W * 0.5, H * 0.42, 20, W * 0.5, H * 0.5, Math.max(W, H) * 0.75);
    gradient.addColorStop(0, '#0e1a26');
    gradient.addColorStop(1, '#060a0e');
    ctx.fillStyle = gradient;
    ctx.fillRect(0, 0, W, H);

    const scene = this.scene;
    if (!scene) return;

    this._pickables = [];
    const s = this._scale();
    const items = [];

    // depth normalisation for shading
    let minZ = Infinity;
    let maxZ = -Infinity;
    const track = (z) => { if (z < minZ) minZ = z; if (z > maxZ) maxZ = z; };
    const radius = this._radius || 40;
    track(-radius); track(radius);
    const spanZ = Math.max(1e-6, maxZ - minZ);
    const bright = (z) => 0.62 + 0.38 * ((z - minZ) / spanZ);

    // ---- pocket point cloud (behind everything)
    if (this.options.showPocket && this._pocketBuckets.length) {
      for (const bucket of this._pocketBuckets) {
        const path = new Path2D();
        const coords = bucket.coords;
        for (let i = 0; i < coords.length; i += 3) {
          const [x, y] = this._project([coords[i], coords[i + 1], coords[i + 2]]);
          if (x < -4 || y < -4 || x > W + 4 || y > H + 4) continue;
          path.rect(x - 1.3, y - 1.3, 2.6, 2.6);
        }
        ctx.fillStyle = withAlpha(bucket.color, 0.30);
        ctx.fill(path);
      }
    }

    // ---- receptor C-alpha trace
    if (this.options.showTrace) {
      const residues = scene.receptor?.residues || [];
      const ca = scene.receptor?.ca || [];
      for (let i = 0; i < ca.length; i += 1) {
        const [x, y, z] = this._project(ca[i]);
        const res = residues[i];
        const base = res?.sub ? SUB_COLORS[res.sub] || '#5f7d99' : '#4c6479';
        items.push({ kind: 'trace', x, y, z, r: Math.max(2.0, 2.0 * s), color: shade(base, bright(z)), label: res ? `${res.one}${res.resid}` : null, pick: false });
        if (res) this._pickables.push({ x, y, z, r: 7, label: `MurJ ${res.one}${res.resid}`, kind: 'residue', resid: res.resid });
      }
    }

    // ---- pocket-lining side chains
    if (this.options.showSidechains && scene.receptor?.atoms?.length) {
      const atoms = scene.receptor.atoms;
      for (const [a, b] of scene.receptor.bonds) {
        const p1 = this._project(atoms[a].p);
        const p2 = this._project(atoms[b].p);
        const z = (p1[2] + p2[2]) / 2;
        items.push({
          kind: 'stick', x: p1[0], y: p1[1], z: z - 0.01, x2: p2[0], y2: p2[1],
          width: Math.max(1.4, 0.22 * s), color: shade('#8fa3bf', bright(z) * 0.95),
        });
      }
      for (const atom of atoms) {
        const [x, y, z] = this._project(atom.p);
        items.push({
          kind: 'sphere', x, y, z, r: Math.max(1.3, 0.3 * s),
          color: shade(atom.rgb, bright(z)), label: `${atom.one}${atom.resid} ${atom.a}`,
          atom,
        });
        this._pickables.push({ x, y, z, r: 7, label: `${atom.one}${atom.resid} ${atom.a} (${atom.e})`, kind: 'atom' });
      }
    }

    // ---- peptides
    if (this.options.showPeptides) {
      for (const pep of scene.peptides || []) {
        if (this.options.peptideFilter[pep.tag] === false) continue;
        const base = pep.color;
        // CA ribbon
        for (let i = 0; i < pep.ca.length; i += 1) {
          const [x, y, z] = this._project(pep.ca[i]);
          items.push({ kind: 'trace', x, y, z, r: Math.max(2.2, 2.1 * s), color: shade(base, bright(z)) });
        }
        if (pep.atoms?.length) {
          for (const [a, b] of pep.bonds || []) {
            const p1 = this._project(pep.atoms[a].p);
            const p2 = this._project(pep.atoms[b].p);
            const z = (p1[2] + p2[2]) / 2;
            items.push({
              kind: 'stick', x: p1[0], y: p1[1], z: z - 0.01, x2: p2[0], y2: p2[1],
              width: Math.max(1.2, 0.2 * s), color: shade(base, bright(z) * 0.9),
            });
          }
        }
      }
    }

    // ---- ligand
    if (this.options.showLigand && scene.ligand) {
      const atoms = scene.ligand.atoms;
      for (const [a, b] of scene.ligand.bonds || []) {
        const p1 = this._project(atoms[a].p);
        const p2 = this._project(atoms[b].p);
        const z = (p1[2] + p2[2]) / 2;
        items.push({
          kind: 'stick', x: p1[0], y: p1[1], z: z - 0.02, x2: p2[0], y2: p2[1],
          width: Math.max(1.8, 0.16 * s), color: shade('#cfd8e3', bright(z)),
        });
      }
      for (const atom of atoms) {
        const [x, y, z] = this._project(atom.p);
        items.push({
          kind: 'sphere', x, y, z, r: Math.max(2.4, this.options.ligandRadius * s),
          color: shade(atom.rgb, bright(z) * 1.08), ligand: true,
        });
        this._pickables.push({ x, y, z, r: 8, label: `Ligand ${atom.e}`, kind: 'ligand' });
      }
    }

    // ---- hotspots
    if (this.options.showHotspots) {
      for (const h of scene.hotspots || []) {
        const [x, y, z] = this._project(h.p);
        const score = Number(h.score) || 0;
        items.push({
          kind: 'sphere', x, y, z, r: Math.max(2.0, 0.5 * s + score * 2.6),
          color: shade(h.evidence ? '#ff8a5b' : '#5b7fa6', bright(z) * (0.85 + score * 0.4)),
        });
        this._pickables.push({
          x, y, z, r: 9, label: `${h.one}${h.resid} Hotspot ${h.score}`, kind: 'hotspot', hotspot: h,
        });
      }
    }

    // ---- features
    if (this.options.showFeatures) {
      for (const f of scene.features || []) {
        const [x, y, z] = this._project(f.p);
        items.push({ kind: 'ring', x, y, z, r: Math.max(6, 1.1 * s), color: f.rgb, label: f.fid });
      }
    }

    items.sort((a, b) => a.z - b.z);
    const hovered = this.hover;
    for (const item of items) {
      if (item.kind === 'trace' || item.kind === 'sphere') {
        ctx.beginPath();
        ctx.arc(item.x, item.y, item.r, 0, TAU);
        ctx.fillStyle = item.color;
        ctx.fill();
        if (item.kind === 'sphere' && item.r > 3) {
          ctx.beginPath();
          ctx.arc(item.x - item.r * 0.3, item.y - item.r * 0.3, item.r * 0.42, 0, TAU);
          ctx.fillStyle = 'rgba(255,255,255,.20)';
          ctx.fill();
        }
        if (hovered && hovered.kind === 'hotspot' && Math.abs(hovered.x - item.x) < 1
            && Math.abs(hovered.y - item.y) < 1) {
          ctx.beginPath();
          ctx.arc(item.x, item.y, item.r + 4, 0, TAU);
          ctx.strokeStyle = '#ffffff';
          ctx.lineWidth = 1.5;
          ctx.stroke();
        }
      } else if (item.kind === 'stick') {
        ctx.beginPath();
        ctx.moveTo(item.x, item.y);
        ctx.lineTo(item.x2, item.y2);
        ctx.strokeStyle = item.color;
        ctx.lineWidth = item.width;
        ctx.lineCap = 'round';
        ctx.stroke();
      } else if (item.kind === 'ring') {
        ctx.beginPath();
        ctx.arc(item.x, item.y, item.r, 0, TAU);
        ctx.fillStyle = withAlpha(item.color, 0.16);
        ctx.fill();
        ctx.lineWidth = 1.4;
        ctx.strokeStyle = withAlpha(item.color, 0.9);
        ctx.stroke();
      }
    }

    // ---- labels on top
    if (this.options.showLabels) {
      ctx.font = '600 10px ui-monospace, Consolas, monospace';
      ctx.textAlign = 'center';
      for (const f of scene.features || []) {
        const [x, y, z] = this._project(f.p);
        ctx.fillStyle = 'rgba(6,10,14,.72)';
        const text = f.fid;
        const w = ctx.measureText(text).width + 8;
        ctx.fillRect(x - w / 2, y - 17, w, 13);
        ctx.fillStyle = withAlpha(f.rgb, 0.35 + 0.65 * ((z - minZ) / spanZ));
        ctx.fillText(text, x, y - 7);
      }
      ctx.textAlign = 'left';
    }

    // ---- hovered residue name in the corner (screen reader friendly too)
    if (hovered) {
      ctx.font = '600 11px ui-monospace, Consolas, monospace';
      ctx.fillStyle = 'rgba(219,228,238,.92)';
      ctx.fillText(hovered.label || '', 12, H - 12);
    }

    if (this.options.spin) {
      this.view.rot = matMul(matRotY(0.006), this.view.rot);
      this.requestRender();
    }
  }

  _startSpin() {
    this._stopSpin();
    this.options.spin = true;
    this.requestRender();
  }

  _stopSpin() {
    this.options.spin = false;
    if (this._spinTimer) clearInterval(this._spinTimer);
    this._spinTimer = null;
  }
}

const SUB_COLORS = {
  S1_basic_bowl: '#4f8ef7',
  S2_carboxylate_gate: '#e55353',
  S3_TM2_TM7_groove: '#38a169',
  S4_TM7_floor: '#c99a3a',
};

function meanOf(points) {
  const sum = [0, 0, 0];
  for (const p of points) { sum[0] += p[0]; sum[1] += p[1]; sum[2] += p[2]; }
  const n = points.length || 1;
  return [sum[0] / n, sum[1] / n, sum[2] / n];
}

export { SUB_COLORS, shade };
