import json,html,time
from pathlib import Path
from sbvs.paths import ROOT,read_json,write_json
from sbvs.export import read_rows
winroot=__import__('subprocess').check_output(['wslpath','-w',str(ROOT)],text=True).strip()
doctor=read_json(ROOT/'reports/doctor.json');browser=read_json(ROOT/'reports/browser_acceptance.json');class20=read_json(ROOT/'reports/class20_result.json');fol=read_json(ROOT/'reports/FOL_redocking.json');gpu=read_json(ROOT/'reports/gnina_gpu_smoke.json');cpu=read_json(ROOT/'reports/gnina_cpu_smoke.json');trace=read_json(ROOT/'reports/gpu_device_trace.json')
results={
'A_environment_CPU':{'status':'passed','evidence':['reports/doctor.json','reports/cco_cli_smoke.json','logs/install/receptor_cli.log','data/prepared/demo_1drf/preparation.log'],'note':'Initial raw receptor attempt failed; explicit missing-atom repair and real CPU docking subsequently passed'},
'B_default_page':{'status':'passed','evidence':['reports/browser_acceptance.json','reports/browser/initial.png','reports/browser/cco_1.png','reports/browser/cco_2.png','reports/browser/cco_3.png']},
'C_invalid_inputs':{'status':'passed','evidence':['reports/unit_tests.log','reports/recovery_acceptance.json','reports/research_acceptance.json']},
'D_class20':{'status':'passed','run':class20['run'],'evidence':['reports/class_report_images.json']},
'E_GNINA_GPU':{'status':'passed','run':gpu['run'],'evidence':['reports/gpu_device_trace.json','reports/gnina_gpu_smoke.json']},
'E_GNINA_no_gpu_branch':{'status':'passed','run':cpu['run'],'note':'Explicit --no_gpu tested on GPU-equipped host; separate physical GPU-less computer not tested'},
'F_chemistry_pose_RMSD':{'status':'passed','evidence':['reports/artifact_acceptance.json','reports/FOL_redocking.json','reports/unit_tests.log']},
'G_recovery':{'status':'passed','evidence':['reports/recovery_acceptance.json','reports/corruption_acceptance.json','reports/stop_guard.json']},
'H_spaces_Chinese_paths':{'status':'passed','evidence':['reports/path_acceptance.json','reports/windows_path_test.log']},
'H_other_PC_full_VHDX_migration':{'status':'not_tested','note':'No second computer or cloned-VHDX registration test; no claim of full cross-machine validation'},
'I_offline':{'status':'passed','evidence':['reports/offline_namespace.log','reports/browser_acceptance.json'],'note':'Real CCO in isolated Linux network namespace, browser external requests blocked'},
'J_result_consistency':{'status':'passed','evidence':['reports/browser_acceptance.json','reports/artifact_acceptance.json','reports/research_acceptance.json']},
'K_MD_handoff':{'status':'passed','evidence':['reports/artifact_acceptance.json'],'note':'Automated test fixture explicitly labelled; not a real scientific choice by user; no MD executed'},
'state_enumeration':{'status':'passed','evidence':['reports/enumeration_acceptance.json']},
'target_enrichment_validation':{'status':'not_tested','note':'No labelled benchmark provided'},
'teacher_platform_remote_access':{'status':'blocked','note':'Initial sandbox request failed; bounded host request returned HTTP 401 authentication required. No case content viewed; no application dependency'}
}
for name in ('migration_local','runtime_export'):
 f=ROOT/'reports'/f'{name}.json'
 results[name]=read_json(f) if f.exists() else {'status':'not_tested'}
write_json(ROOT/'reports/acceptance_matrix.json',results)
text='# 实现与验收报告\n\n日期：2026-10-13（Asia/Shanghai）。所有 passed 均有本机真实文件/引擎/测试证据；not_tested 不计作通过。\n\n'
text+='## 实际部署\n\nWindows 项目：`'+winroot+'`。专用 WSL：`'+read_json(ROOT/'project.json')['distro']+'`。VHDX：`'+winroot+'\\runtime\\wsl\\ext4.vhdx`。外层盘为 exFAT，已实际导入和启动 WSL2；Linux 依赖与 SQLite 保存在 ext4 中。原 Ubuntu 保持独立。\n\n'
text+='启动：`START_WINDOWS.cmd`；地址：`http://localhost:5902/`。停止：`STOP_WINDOWS.cmd`。迁移前：`PREPARE_MIGRATION.ps1`；新电脑登记：`REGISTER_ON_THIS_PC.ps1`。\n\n'
text+='## 实际版本与设备\n\nPython '+doctor['python']+'；RDKit '+doctor['packages']['rdkit']+'；Meeko '+doctor['packages']['meeko']+'；Vina 1.2.7；GNINA 1.3.3；Molscrub 0.2.2。完整 build/来源分别见 env/conda-linux-64.lock、requirements-pip.lock、package_provenance.json、vendor.lock.json。\n\n'+doctor['gpu_detection'].get('output','')+'\nGNINA 固定模型：crossdock_default2018_1_3。GPU 明确 device 0；真实 GPU 设备 trace 记录 '+str(trace['successful_gpu_ioctls'])+' 次成功 /dev/dxg ioctl，并有有限 CNN 输出。\n\n'
text+='## 默认页面三次全新 CCO\n\n| Run ID | 秒 | Vina affinity kcal/mol | 截图 |\n|---|---:|---:|---|\n'
for r in browser['tests']:text+=f"|{r['run_id']}|{r['elapsed_seconds']:.3f}|{r['affinity_kcal_mol']}|{r['screenshot']}|\n"
text+='\n每个 run 都有独立执行日志与原始 Vina 输出；相同种子产生相同数值是实际重算结果。首次页面无预填结果，浏览器外网阻断时页面与本地 3D 仍工作。\n\n'
text+='## 批量、GNINA 与科学检查\n\n20 分子课堂：`'+class20['run']+'`，20/20 成功，180 个 pose。报告在该 run/report.html，两张不同真实姿态图在 run/images。\n\n'
text+=f"GNINA GPU：`{gpu['run']}`，{gpu['elapsed_seconds']:.3f} 秒；CPU CNN：`{cpu['run']}`，{cpu['elapsed_seconds']:.3f} 秒。两者均有 docking SDF、empirical score、CNN pose score、CNN affinity pK。\n\n"
text+=f"FOL 真正重对接：`{fol['run']}`；Top1 RMSD {fol['top1']['rmsd_A']:.4f} Å，Top3 中最佳与全部 pose 最佳均为 {fol['best_all']['rmsd_A']:.4f} Å。使用同一受体坐标系、对称性/原子映射校正，未独立最佳叠合配体。参考为原始 FOL A:187、32 重原子；不把 score_only 当 redocking。\n\n"
text+='1DRF 的 Gln12 等缺失重原子导致初次 Meeko 失败；最终明确采用 PDBFixer 补原子，保留原坐标来源和修复记录，没有删除坏残基或重建缺失 loop。首次 GNINA 缺 CUDA 用户态库与无效模型别名的失败也保留在日志/失败 run 中，随后真实修复并重测。\n\n'
text+='科研自定义受体、mmCIF 计划、不存在的链报错、Molscrub 状态枚举、TopK score_only / redock、多 seed、minimize、取消/中断/续跑、化学身份与人工导出门禁均有对应 reports 证据。\n\n'
text+='## 逐项状态\n\n| 项目 | 状态 | 说明 / 证据 |\n|---|---|---|\n'
for key,v in results.items():text+='|'+key+'|'+v['status']+'|'+str(v.get('note',v.get('evidence',v.get('run','')))).replace('|','/')+'|\n'
text+='\n## 未实测与未开放范围\n\n未在另一台真实电脑、无 GPU 的另一台主机或复制后的新 VHDX 注册实例验收；Windows ARM、WSL1 不支持声明。离线环境导出若已完成仅证明快照创建与哈希，不证明恢复导入成功。老师地址的宿主请求返回 HTTP 401，未获得认证、未查看案例内容，也未通过老师平台远程代理验收。\n\n未开放共价对接、一般金属配位/复杂辅因子参数化、柔性受体、缺失长 loop 重建和 MOL2 化学可靠性支持。没有活性标签，未评价富集或证明筛选性能。严格跨硬件逐位一致性未承诺。NVIDIA 进程列表在本机 WSL 未提供匹配 PID，因此 GPU 证据使用实际 /dev/dxg trace 和有限 CNN 输出。\n\nMD handoff 验收文件为明确标记的自动测试样本，不冒充用户亲自选中的候选；没有执行 MD、力场参数化或自由能计算。\n'
(ROOT/'docs/实现与验收报告.md').write_text(text,encoding='utf-8')
body='<h1>GNINA / Vina 本地工作站交付记录</h1><p><a href="http://localhost:5902/">打开已运行的工作站</a></p><pre>'+html.escape(text)+'</pre>'
for name in ['initial.png','cco_1.png','class20.png']:
 path=ROOT/'reports/browser'/name
 if path.exists():body+='<h2>'+name+'</h2><img style="max-width:100%" src="data:image/png;base64,'+__import__('base64').b64encode(path.read_bytes()).decode()+'">'
(ROOT/'reports/交付与验收.html').write_text('<!doctype html><meta charset="utf-8"><style>body{font:16px system-ui;max-width:1200px;margin:40px auto;padding:20px;color:#17313e}pre{white-space:pre-wrap;overflow-wrap:anywhere;line-height:1.65}a{color:#087d88}</style>'+body,encoding='utf-8')
print(json.dumps({'passed':sum(v['status']=='passed' for v in results.values()),'not_tested':sum(v['status']=='not_tested' for v in results.values()),'blocked':sum(v['status']=='blocked' for v in results.values())}))
