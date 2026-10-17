/* pages/docs.js -- project documentation browser (methodology, interfaces,
   evidence chain, work log). */

import { api } from '../api.js';
import {
  el, div, span, card, markdownNode, loadingPane, fmtBytes, fmtDate, chip,
} from '../ui.js';

export async function mount(root, ctx) {
  const data = await api.docs();
  const docs = (data.docs || []).slice().sort((a, b) => {
    if (a.folder !== b.folder) return a.folder === '.' ? 1 : -1;
    return a.name.localeCompare(b.name);
  });
  if (!docs.length) {
    root.replaceChildren(card(null, div({ class: 'empty-state', text: '没有找到 Markdown 文档。' })));
    return null;
  }
  let current = ctx.params.doc ? docs.find((d) => d.name === ctx.params.doc) || docs[0] : docs[0];

  const list = div({ style: { display: 'grid', gap: '2px' } });
  const content = card(null, loadingPane('正在载入文档…'));
  const crumbNote = span({ class: 'tiny muted' });

  root.replaceChildren();
  root.appendChild(div({ class: 'page-head' },
    div(null,
      el('h1', { text: '文档' }),
      div({ class: 'sub', text: '方法学主文档、下游对接协议、证据链说明与工作日志。全部为仓库内原始 Markdown，不做改写。' }))));

  function paintList() {
    list.replaceChildren();
    let lastFolder = null;
    for (const doc of docs) {
      if (doc.folder !== lastFolder) {
        lastFolder = doc.folder;
        list.appendChild(div({ class: 'nav-group', text: doc.folder === '.' ? '根目录' : doc.folder }));
      }
      const active = current && doc.name === current.name;
      list.appendChild(el('a', {
        href: `#/docs?doc=${encodeURIComponent(doc.name)}`,
        style: {
          display: 'block', padding: '5px 8px', borderRadius: '6px', fontSize: '12px',
          color: active ? '#bfe8ff' : 'var(--text-dim)',
          background: active ? 'rgba(56,189,248,.1)' : 'transparent',
          textDecoration: 'none',
        },
        onclick: (event) => {
          event.preventDefault();
          show(doc);
        },
      },
        div({ text: doc.title || doc.name }),
        div({ class: 'tiny muted mono', text: `${doc.name} · ${fmtBytes(doc.size)}` })));
    }
  }

  async function show(doc) {
    current = doc;
    paintList();
    content.replaceChildren(loadingPane(`正在载入 ${doc.name}…`));
    try {
      const payload = await api.doc(doc.name);
      content.replaceChildren(
        div({ class: 'split', style: { marginBottom: '8px' } },
          chip(doc.name, 'mono'),
          span({ class: 'tiny muted', text: `${fmtBytes(doc.size)} · ${fmtDate(doc.mtime)}` }),
          span({ class: 'right' },
            el('a', { class: 'tiny', href: `/files/${doc.folder === '.' ? '' : `${doc.folder}/`}${doc.name}?download=1`, text: '下载原文 ↘' }))),
        markdownNode(payload.markdown || ''));
      history.replaceState(null, '', `#/docs?doc=${encodeURIComponent(doc.name)}`);
    } catch (err) {
      content.replaceChildren(card(null, div({ class: 'empty-state', text: `载入失败：${err.message}` })));
    }
  }

  root.appendChild(div({ class: 'grid wide-right', style: { alignItems: 'start' } },
    card({ class: 'tight', title: '目录' }, list, div({ class: 'hr' }), crumbNote),
    content));

  paintList();
  await show(current);
  return null;
}
