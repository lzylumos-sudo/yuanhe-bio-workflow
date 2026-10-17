# -*- coding: utf-8 -*-
"""webui/test_api.py -- end-to-end checks of the web console's JSON API.

Starts the real HTTP server on an ephemeral port (in-process, no external
process) and asserts the shape and the semantic content of every endpoint.

    python webui/test_api.py
"""
from __future__ import annotations

import json
import os
import sys
import threading
import time
import urllib.error
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

try:                                    # the console may be cp936 on Windows
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except (AttributeError, ValueError):
    pass

from webui import server                                            # noqa: E402

BASE = None
PASSED = 0
FAILED = []


def check(name, condition, detail=''):
    global PASSED
    if condition:
        PASSED += 1
        print(f'  ok   {name}')
    else:
        FAILED.append(name)
        print(f'  FAIL {name} {detail}')


def fetch(path, method='GET', body=None):
    url = BASE + path
    data = None
    headers = {'Accept': 'application/json'}
    if body is not None:
        data = json.dumps(body).encode('utf-8')
        headers['Content-Type'] = 'application/json'
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(req, timeout=180) as resp:
            raw = resp.read()
            return resp.status, raw, resp.headers.get('Content-Type', '')
    except urllib.error.HTTPError as exc:
        return exc.code, exc.read(), exc.headers.get('Content-Type', '')


def jget(path):
    status, raw, ctype = fetch(path)
    assert status in (200, 202), f'{path} -> {status}: {raw[:200]}'
    return json.loads(raw.decode('utf-8'))


def run(do_job=False):
    global BASE
    httpd = server.Server(('127.0.0.1', 0), server.Handler)
    BASE = f'http://127.0.0.1:{httpd.server_address[1]}'
    threading.Thread(target=httpd.serve_forever, kwargs=dict(poll_interval=0.2),
                     daemon=True).start()
    print(f'[api] {BASE}')

    # ---- health / meta
    print('\n· 基础')
    meta = jget('/api/meta')
    check('meta 有工作目录与结构定义', meta.get('root') and len(meta.get('structures', {})) == 3)
    check('meta 报告 results 与 evidence 可用',
          meta['available']['results'] and meta['available']['evidence'])
    check('meta 列出 figures', len(meta.get('figures', [])) >= 4)

    status, raw, ctype = fetch('/')
    check('首页是 HTML', status == 200 and 'text/html' in ctype and b'<div id="app"' in raw)
    status, raw, ctype = fetch('/js/app.js')
    check('前端模块以 JS MIME 提供', status == 200 and 'javascript' in ctype)

    # ---- overview
    print('\n· 概览')
    ov = jget('/api/overview')
    fun = ov['funnel']
    check('漏斗数字与 README 一致',
          fun['enumerated'] == 20307 and fun['drug_like'] == 538,
          f'(enumerated={fun.get("enumerated")}, drug_like={fun.get("drug_like")})')
    check('最佳骨架是联苯羧酸骨架',
          ov['best'][0]['smiles'].startswith('OC(=O)c1ccc'), ov['best'][0]['smiles'])
    check('证据状态 24 条、PASS 22', ov['evidence']['n'] == 24
          and ov['evidence']['by_status'].get('PASS') == 22, str(ov['evidence']))
    check('结构过滤 538 -> 32', ov['filters']['n_input'] == 538
          and ov['filters']['survivors'] == 32)

    # ---- design data
    print('\n· 设计与结果')
    hot = jget('/api/hotspots')
    check('界面热点 52 个、三种肽共同接触 14 个',
          hot['n'] == 52 and hot['n_common'] == 14, f"(n={hot['n']}, common={hot['n_common']})")
    scores = [r['score'] for r in hot['rows']]
    check('热点已按 score 降序', scores == sorted(scores, reverse=True))

    ph = jget('/api/pharmacophore')
    check('12 个共识药效团特征', len(ph['features']) == 12, str(len(ph['features'])))
    check('每个特征都有 3D 坐标与类型',
          all(len(f['pos']) == 3 and f['type'] for f in ph['features']))
    check('特征对带距离', all(p['distance'] > 3 for p in ph['pairs']))

    axes = jget('/api/axes')['axes']
    # A7-H12 has no surviving enumeration (its distance/weight combination is
    # outside the linker budget) -- it must still be listed, with candidates = 0.
    check('设计轴带枚举统计（A7-H12 无候选也保留）',
          len(axes) == 15
          and sum(1 for a in axes if a['candidates'] > 0) == 14
          and any(a['candidates'] == 0 for a in axes), str(len(axes)))

    pocket = jget('/api/pocket')['pocket']
    check('腔体 5,440 Å³、4 个亚口袋',
          round(pocket['total_volume_A3']) == 5440 and len(pocket['subpockets']) == 4,
          f"{pocket['total_volume_A3']}, {len(pocket['subpockets'])}")

    sc = jget('/api/scaffolds')
    check('200 个排名骨架', len(sc['rows']) == 200, str(len(sc['rows'])))
    row = sc['rows'][0]
    check('骨架行含 SMILES/性质/打分/接触残基',
          row['smiles'] and row['MW'] and 'total' in row and row['murj_contacts'])
    check('骨架行与结构过滤按 SMILES 连接',
          any('f_keep' in r for r in sc['rows']))
    check('过滤判定是布尔值而不是字符串',
          all(isinstance(r['f_keep'], bool) for r in sc['rows'] if 'f_keep' in r)
          and sum(1 for r in sc['rows'] if r.get('f_keep') is True) == 14)
    check('排名第一的骨架并未通过结构过滤（诚实的负结果）',
          row['f_keep'] is False and row['f_bridge'] is True,
          f"keep={row['f_keep']!r} bridge={row['f_bridge']!r}")

    filt = jget('/api/filters')
    check('结构过滤 538 行、存活 32', len(filt['rows']) == 538
          and filt['summary']['survivors'] == 32)
    check('全部过滤行 keep 为布尔且 32 个为真',
          sum(1 for r in filt['rows'] if r['keep'] is True) == 32
          and sum(1 for r in filt['rows'] if r['bridge'] is True) == 462)

    poses = jget('/api/poses')
    check('25 个 Top 位姿且原子/键齐全',
          len(poses['top']) == 25 and poses['top'][0]['atoms']
          and poses['top'][0]['bonds'])
    check('位姿坐标在 9NU8 坐标系内（Z > 140）',
          all(a['p'][2] > 140 for a in poses['top'][0]['atoms']))

    # ---- evidence
    print('\n· 证据链')
    rep = jget('/api/evidence/report')
    check('报告含 24 条先验核对章节',
          '先验' in rep['markdown'] and '## 5' in rep['markdown'])
    pri = jget('/api/evidence/priors')
    check('先验 24 条、PASS 22 / FLAG 1 / WEAK 1',
          len(pri['rows']) == 24 and pri['status']['by_status'] ==
          {'PASS': 22, 'FLAG': 1, 'WEAK': 1}, str(pri['status']['by_status']))
    check('先验行带原文引用与出处',
          pri['rows'][0]['claim'] and pri['rows'][0]['locator'])
    cen = jget('/api/evidence/census')
    check('普查：MurJ 496 残基（构建体已剔除）',
          all(c['murj_residues'] == 496 for c in cen['atom_budget']['census']))
    check('原子预算：肽接触 121 vs 骨架中位 24',
          cen['atom_budget']['budget']['small_molecule']['heavy_atoms_median'] == 24)
    iface = jget('/api/evidence/interface')
    check('界面测量有双边数据', len(iface['murj']) > 50 and len(iface['peptide']) > 50)
    check('锚点表非空', len(iface['anchors']) >= 12)

    # ---- 3D scene
    print('\n· 3D 结构')
    scene = jget('/api/scene?ref=SglCJ3&ligand=0')
    check('受体重原子轨迹 496 个 Cα', len(scene['receptor']['ca']) == 496)
    check('三肽都在同一坐标系内', len(scene['peptides']) == 3
          and all(len(p['ca']) > 20 for p in scene['peptides']))
    check('配体位姿 24 个重原子', len(scene['ligand']['atoms']) == 24)
    check('药效团与热点随场景下发',
          len(scene['features']) == 12 and len(scene['hotspots']) >= 20)
    check('亚口袋配色可用于轨迹着色',
          len(scene.get('subpockets', [])) == 4
          and any(r.get('sub') for r in scene['receptor']['residues']))
    check('叠合后 SglM 肽落在 SglCJ3 口袋附近',
          _dist(scene['peptides'][0]['ca'][0], scene['center']) < 60)

    scene_m = jget('/api/scene?ref=SglM')
    check('换参考结构可用', len(scene_m['receptor']['ca']) == 496
          and scene_m['ref'] == 'SglM')

    pts = jget('/api/pocket/points')
    if pts['status'] == 'ready':
        check('腔体点云 10k+ 体素、体积 5,440 Å³',
              pts['n_voxels'] > 9000 and round(pts['volume']) == 5440)
        check('点云按化学分类', len(pts['chem']) == pts['n_voxels'])
    else:
        check('腔体点云正在计算（首次约 15 s）', pts['status'] == 'computing')

    # ---- docs & jobs
    print('\n· 文档与任务')
    docs = jget('/api/docs')['docs']
    check('列出 5 篇以上文档', len(docs) >= 5, str(len(docs)))
    check('文档有标题', all(d['title'] for d in docs))
    doc = jget('/api/docs/' + urllib.parse.quote('01_分子骨架设计算法.md'))
    check('中文文件名文档可读取', len(doc['markdown']) > 5000)

    jobs = jget('/api/jobs')
    check('任务接口返回四类可运行脚本', len(jobs['kinds']) == 4)
    check('运行历史结构完整', isinstance(jobs['jobs'], list))

    status, raw, _ = fetch('/api/jobs', 'POST', {'kind': 'nope'})
    check('未知任务类型被拒绝', status == 400)

    status, raw, _ = fetch('/api/jobs', 'POST',
                           {'kind': 'pipeline', 'params': {'out': '../escape'}})
    check('非法输出目录被拒绝', status == 400)

    status, raw, _ = fetch('/api/jobs', 'POST',
                           {'kind': 'filters', 'params': {'targets': -5}})
    check('越界参数被拒绝', status == 400)

    # ---- file serving / traversal
    print('\n· 静态文件')
    status, raw, ctype = fetch('/files/figures/fig4_funnel.png')
    check('图片可下载', status == 200 and 'image/png' in ctype)
    status, raw, _ = fetch('/files/priors/murj_priors.json')
    check('先验文件可读取', status == 200 and b'structures' in raw)
    status, raw, _ = fetch('/files/webui/data.py')
    check('webui 目录不可越权读取', status == 403, str(status))
    status, raw, _ = fetch('/files/../webui/data.py')
    check('路径穿越被拒绝', status in (403, 404), str(status))
    status, raw, ctype = fetch('/nonexistent-page')
    check('未知路径回落到 SPA 首页', status == 200 and 'text/html' in ctype)
    status, raw, _ = fetch('/js/missing-file.js')
    check('缺失的前端资源返回 404 而不是 HTML', status == 404, str(status))

    # ---- 7-step workflow UI
    print('\n· 工作流（7 步）')
    wf = jget('/api/workflow')
    st = wf['stages']
    check('工作流共 7 步', [s['n'] for s in st] == list(range(1, 8)))
    check('第 1/2 步已接入且数字来自文件',
          st[0]['status'] == 'ready' and st[1]['status'] == 'ready'
          and st[1]['funnel'][0]['count'] == 20307 and st[1]['out'] == 32)
    check('第 3-7 步无数据时为 pending 且计数为空',
          all(s['status'] == 'pending' and all(f['count'] is None for f in s['funnel']) for s in st[2:]))
    wd = jget('/api/workflow?demo=1')
    check('布局预览模式标记为 demo', all(s['status'] == 'demo' for s in wd['stages'][2:]))
    mols = jget('/api/stage/scaffold/molecules')
    check('第 2 步分子列表 538 条、32 条通过',
          mols['n'] == 538 and sum(1 for r in mols['rows'] if r['keep']) == 32)
    pose = jget('/api/stage/scaffold/pose?id=' + mols['rows'][0]['id'])
    check('第 2 步分子有 3D 位姿', pose['ligand'] and len(pose['ligand']['atoms']) > 5
          and len(pose['ligand']['bonds']) > 5)
    check('未接入阶段分子列表为空', jget('/api/stage/docking/molecules')['n'] == 0)
    status, raw, _ = fetch('/api/stage/nope/molecules')
    check('未知阶段返回 404', status == 404, str(status))
    m2 = jget('/api/mol2d?smiles=' + 'OC(%3DO)c1ccccc1')
    check('2D 结构图可生成', bool(m2.get('svg')) or 'RDKit' in (m2.get('reason') or ''))

    # ---- live job (opt-in: runs the repo's own smoke test through the API)
    if do_job:
        print('\n· 实时任务（冒烟测试）')
        status, raw, _ = fetch('/api/jobs', 'POST', {'kind': 'tests'})
        check('任务已受理', status == 201, f'{status}: {raw[:120]}')
        job_id = json.loads(raw.decode('utf-8'))['id']
        seen_running = False
        deadline = time.time() + 900
        final = None
        while time.time() < deadline:
            time.sleep(1.0)
            listing = jget('/api/jobs')['jobs']
            if any(j['id'] == job_id and j['status'] == 'running' for j in listing):
                seen_running = True
            current = jget(f'/api/jobs/{job_id}?since=0')
            if current['status'] in ('done', 'failed', 'cancelled'):
                final = current
                break
        check('任务列表在运行期间可读（不返回进程对象）',
              seen_running and all('proc' not in j for j in listing))
        check('任务执行完成且退出码 0',
              final is not None and final['status'] == 'done' and final['exit'] == 0,
              str(final and (final['status'], final['exit'])))
        check('日志已流式回传', final is not None and final['total_lines'] > 5,
              str(final and final['total_lines']))

    httpd.shutdown()


def _dist(a, b):
    return sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5


if __name__ == '__main__':
    import argparse

    ap = argparse.ArgumentParser(description='WEDGE-Design 控制台 API 测试')
    ap.add_argument('--job', action='store_true',
                    help='额外跑一次实时任务（tests/test_smoke.py，约 1 分钟）')
    args = ap.parse_args()
    try:
        run(do_job=args.job)
    except Exception as exc:                                       # noqa: BLE001
        import traceback
        traceback.print_exc()
        FAILED.append(f'crash: {exc}')
    print(f'\n{PASSED} passed, {len(FAILED)} failed')
    if FAILED:
        print('failed:', ', '.join(FAILED))
    raise SystemExit(1 if FAILED else 0)
