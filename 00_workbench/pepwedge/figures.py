# -*- coding: utf-8 -*-
"""pepwedge.figures -- publication-style figures without PyMOL.

Everything is projected onto the membrane plane (the PCA normal of MurJ), which
makes the periplasm-open cavity, the three peptide paths and the designed
scaffold readable in a single 2D panel.
"""
from __future__ import annotations

import os
from typing import Dict, List

import os as _os
import tempfile as _tempfile
_os.environ.setdefault('MPLCONFIGDIR', _os.path.join(_tempfile.gettempdir(), 'mpl-cfg-pepwedge'))

import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

for _f in (r'C:/Windows/Fonts/msyh.ttc', r'C:/Windows/Fonts/simhei.ttf'):
    try:
        fm.fontManager.addfont(_f)
    except Exception:
        pass
mpl.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
mpl.rcParams['axes.unicode_minus'] = False
import numpy as np

from .core import apply_rt, membrane_frame, kabsch

CHEM_COLORS = {'pos': '#2b6cb0', 'neg': '#c53030', 'polar': '#38a169',
               'aromatic': '#b7791f', 'hydrophobic': '#a0aec0'}
PEPTIDE_COLORS = {'SglM': '#ecc94b', 'SglPP7': '#dd6b20', 'SglCJ3': '#38a169'}


def _projector(cx):
    n, c0 = membrane_frame(cx)
    e1 = np.cross(n, [0.0, 0.0, 1.0])
    if np.linalg.norm(e1) < 1e-3:
        e1 = np.cross(n, [0.0, 1.0, 0.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    return (lambda X: np.stack([(np.asarray(X) - c0) @ e1, (np.asarray(X) - c0) @ e2], 1)), n


def _ca_trace(cx, chain='receptor'):
    d = cx.receptor if chain == 'receptor' else cx.ligand
    pts = [d[k].atoms['CA'] for k in sorted(d) if 'CA' in d[k].atoms]
    return np.array(pts)


def fig_hotspots(hotspots, path, top=18):
    rows = hotspots[:top][::-1]
    fig, ax = plt.subplots(figsize=(6.2, 5.4), dpi=200)
    labels = [f'{r["one"]}{r["resid"]}' for r in rows]
    vals = [r['score'] for r in rows]
    cols = ['#c53030' if r['evidence'] else '#2b6cb0' for r in rows]
    ax.barh(labels, vals, color=cols, edgecolor='white', height=0.72)
    for i, r in enumerate(rows):
        ax.text(vals[i] + 0.005, i, f'{r["dSASA"]:.0f} Å²', va='center', fontsize=6.5, color='#4a5568')
    ax.set_xlabel('genetics-anchored hotspot score', fontsize=8)
    ax.set_title('MurJ interface hotspots (consensus of 3 phage peptides)', fontsize=9)
    ax.tick_params(labelsize=7)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.text(0.99, 0.02, '#c53030: validated by mutagenesis / resistance',
            transform=ax.transAxes, ha='right', fontsize=6, color='#4a5568')
    fig.tight_layout()
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)


def fig_pocket_map(pk, cxs, features, path, ref='SglCJ3'):
    cx = cxs[ref]
    proj, n = _projector(cx)
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.6), dpi=200,
                             gridspec_kw=dict(width_ratios=[1.45, 1]))
    ax = axes[0]
    P = proj(pk.points)
    cols = [CHEM_COLORS[c] for c in pk.chemistry]
    idx = np.random.RandomState(0).choice(len(P), size=min(4500, len(P)), replace=False)
    ax.scatter(P[idx, 0], P[idx, 1], c=[cols[i] for i in idx], s=3, alpha=0.55, lw=0)
    for tag, c in PEPTIDE_COLORS.items():
        tr = proj(_ca_trace(cxs[tag], 'ligand'))
        if tag != ref:
            ref_cx = cxs[ref]
            shared = [r for r in ref_cx.receptor if 'CA' in ref_cx.receptor[r].atoms
                      and r in cxs[tag].receptor]
            mob = np.array([cxs[tag].receptor[r].atoms['CA'] for r in shared])
            tgt = np.array([ref_cx.receptor[r].atoms['CA'] for r in shared])
            rt = kabsch(mob, tgt)
            tr = proj(apply_rt(_ca_trace(cxs[tag], 'ligand'), rt))
        ax.plot(tr[:, 0], tr[:, 1], color=c, lw=2.4, alpha=0.9,
                label={'SglM': 'Sgl$^{M}$', 'SglPP7': 'Sgl$^{PP7}$', 'SglCJ3': 'Sgl$^{CJ3}$'}[tag])
    ax.plot(proj(_ca_trace(cx))[:, 0], proj(_ca_trace(cx))[:, 1], '.', ms=1.2, color='#2d3748', alpha=0.5)
    for f in features:
        if f.ftype not in ('anion', 'aromatic'):
            continue
        p = proj(f.pos.reshape(1, 3))[0]
        mk = 'o' if f.ftype == 'anion' else 's'
        ax.scatter([p[0]], [p[1]], s=46, marker=mk, facecolor='none',
                   edgecolor='#9b2c2c' if f.ftype == 'anion' else '#975a16', lw=1.4, zorder=5)
        ax.annotate(f.fid, p, xytext=(4, 4), textcoords='offset points', fontsize=6,
                    color='#1a202c')
    ax.set_aspect('equal')
    ax.set_title('A  MurJ periplasm-open pocket + consensus pharmacophore', fontsize=8.5, loc='left')
    ax.legend(fontsize=6.5, loc='lower left', frameon=False)
    ax.set_xticks([]); ax.set_yticks([])

    ax2 = axes[1]
    from .pocket import SUBPOCKETS
    names, vols, cols2 = [], [], []
    for sp in sorted({s for s in pk.subpocket}):
        sel = np.array([s == sp for s in pk.subpocket])
        names.append(sp.replace('_', ' '))
        vols.append(float(sel.sum() * pk.step ** 3))
        chem = [c for c, k in zip(pk.chemistry, sel) if k]
        dom = max(set(chem), key=chem.count)
        cols2.append(CHEM_COLORS[dom])
    ax2.barh(names, vols, color=cols2, edgecolor='white')
    ax2.set_xlabel('volume (Å³)', fontsize=8)
    ax2.set_title('B  sub-pocket partitioning', fontsize=8.5, loc='left')
    ax2.tick_params(labelsize=6.5)
    for s in ('top', 'right'):
        ax2.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)


def fig_scaffold(sc, cxs, rec, path, ref='SglCJ3', label=None):
    cx = cxs[ref]
    proj, n = _projector(cx)
    fig, ax = plt.subplots(figsize=(5.6, 4.6), dpi=200)
    tr = proj(_ca_trace(cx))
    ax.plot(tr[:, 0], tr[:, 1], '-', color='#cbd5e0', lw=1.0)
    lig = proj(sc.pose)
    el = [e.upper() for e in sc.elements]
    cc = ['#e53e3e' if e == 'O' else '#3182ce' if e == 'N' else '#f6ad55' if e == 'S' else '#4a5568'
          for e in el]
    ax.scatter(lig[:, 0], lig[:, 1], s=[46 if e != 'H' else 10 for e in el], c=cc,
               edgecolor='white', lw=0.4, zorder=6)
    d = np.linalg.norm(sc.pose[:, None, :] - rec.xyz[None, :, :], axis=2)
    for j in np.where((d < 3.6).any(0))[0]:
        p = proj(rec.xyz[j:j + 1])[0]
        ax.scatter([p[0]], [p[1]], s=90, facecolor='none', edgecolor='#2f855a', lw=1.1, zorder=4)
    ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(label or 'designed scaffold in the outward-facing cavity', fontsize=8.5, loc='left')
    fig.tight_layout()
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)


def fig_funnel(stats: Dict[str, int], path):
    fig, ax = plt.subplots(figsize=(5.2, 2.6), dpi=200)
    keys = list(stats)
    vals = [stats[k] for k in keys]
    ax.barh(keys[::-1], vals[::-1], color='#3182ce', edgecolor='white')
    for i, v in enumerate(vals[::-1]):
        ax.text(v + max(vals) * 0.01, i, f'{v}', va='center', fontsize=7)
    ax.set_xscale('log')
    ax.set_xlabel('compounds (log scale)', fontsize=8)
    ax.tick_params(labelsize=7)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.set_title('virtual screening funnel', fontsize=9, loc='left')
    fig.tight_layout()
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
