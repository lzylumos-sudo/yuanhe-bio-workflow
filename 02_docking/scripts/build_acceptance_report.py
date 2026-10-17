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
text='# Implementation and acceptance report\n\n Date:2026-10-13(Asia/Shanghai). All passed have local files,/engine/and test evidence;not_tested is not counted as passing.\n\n'
text+='## Actual deployment\n\nWindows Project:`'+winroot+'`. Dedicated WSL: `'+read_json(ROOT/'project.json')['distro']+'`. VHDX: `'+winroot+'\\runtime\\wsl\\ext4.vhdx`. External disk format: exFAT, Imported and started WSL2; Linux Dependencies and SQLite stored in ext4 . Original Ubuntu remains separate.\n\n'
text+='Start:`START_WINDOWS.cmd`; URL:`http://localhost:5902/`. Stop:`STOP_WINDOWS.cmd`. Before migration:`PREPARE_MIGRATION.ps1`; Register on new computer:`REGISTER_ON_THIS_PC.ps1`. \n\n'
text+='## Observed versions and devices\n\nPython '+doctor['python']+'; RDKit '+doctor['packages']['rdkit']+'; Meeko '+doctor['packages']['meeko']+'; Vina 1.2.7; GNINA 1.3.3; Molscrub 0.2.2. Complete build/sources listed in env/conda-linux-64.lock, requirements-pip.lock, package_provenance.json, vendor.lock.json. \n\n'+doctor['gpu_detection'].get('output','')+'\nGNINA Fixed model:crossdock_default2018_1_3. GPU Explicit device 0; Actual GPU Device trace Record '+str(trace['successful_gpu_ioctls'])+' successful runs /dev/dxg ioctl, with finite CNN outputs.\n\n'
text+='## Three fresh default-page CCO\n\n| Run ID | seconds | Vina affinity kcal/mol | screenshots |\n|---|---:|---:|---|\n'
for r in browser['tests']:text+=f"|{r['run_id']}|{r['elapsed_seconds']:.3f}|{r['affinity_kcal_mol']}|{r['screenshot']}|\n"
text+='\n Each run has an independent execution log and raw Vina output; identical seeded values came from actual recomputation. The first page had no prefilled results, and local assets worked with browser network access blocked. 3D remained functional.\n\n'
text+='## Batch,GNINA and scientific checks\n\n20 molecule classroom run:`'+class20['run']+'`, 20/20 succeeded,180 items pose. report stored in that run/report.html, Two distinct actual pose images in run/images. \n\n'
text+=f"GNINA GPU: `{gpu['run']}`, {gpu['elapsed_seconds']:.3f} seconds;CPU CNN: `{cpu['run']}`, {cpu['elapsed_seconds']:.3f} seconds. Both have docking SDF, empirical score, CNN pose score, CNN affinity pK. \n\n"
text+=f"FOL Actual redocking:`{fol['run']}`; Top1 RMSD {fol['top1']['rmsd_A']:.4f} Å, Top3 best pose and all pose best results were {fol['best_all']['rmsd_A']:.4f} Å. using the same receptor frame and symmetry-aware/atom mapping, without independent ligand alignment. Reference: original FOL A:187, 32 heavy atoms; do not treat score_only as redocking. \n\n"
text+='1DRF associated Gln12 and other missing heavy atoms caused the initial Meeko failure; the final repair explicitly used PDBFixer to add atoms, preserving coordinate provenance and repair records, without deleting residues or rebuilding missing loop. Initial GNINA missing CUDA user-space libraries and invalid model aliases remain documented in logs/Failed run and were subsequently repaired and retested.\n\n'
text+='Custom research receptors,mmCIF plans, unknown-chain errors,Molscrub state enumeration,TopK score_only / redock, multiple seed, minimize, cancellation,/interruption,/resume, chemical identity and manual export gates have corresponding reports evidence.\n\n'
text+='## Itemized status\n\n| Project | Status | Description / Evidence |\n|---|---|---|\n'
for key,v in results.items():text+='|'+key+'|'+v['status']+'|'+str(v.get('note',v.get('evidence',v.get('run','')))).replace('|','/')+'|\n'
text+='\n## Untested and unsupported scope\n\n Not accepted on another physical computer, a host without GPU or a copied new VHDX registration;Windows ARM, WSL1 is not claimed supported. Offline export verifies snapshot creation and hashes, not successful restoration. The instructor-host request returned HTTP 401, without authentication or case access; remote-proxy acceptance was not performed.\n\n Covalent docking, general metal coordination,/complex cofactor parameterization, flexible receptors, long missing loop reconstruction and MOL2 chemical reliability are unsupported. No activity labels or enrichment evaluation; no promise of cross-hardware bitwise identity.NVIDIA The local process list for WSL did not provide matching PID, therefore GPU evidence uses actual /dev/dxg trace and finite CNN outputs.\n\nMD handoff Acceptance files are explicitly labeled automated test samples, not user-selected candidates; no MD, force-field parameterization or free-energy calculation was performed.\n'
(ROOT/'docs/Implementation and acceptance report.md').write_text(text,encoding='utf-8')
body='<h1>GNINA / Vina Local workstation delivery record</h1><p><a href="http://localhost:5902/">Open running workstation</a></p><pre>'+html.escape(text)+'</pre>'
for name in ['initial.png','cco_1.png','class20.png']:
 path=ROOT/'reports/browser'/name
 if path.exists():body+='<h2>'+name+'</h2><img style="max-width:100%" src="data:image/png;base64,'+__import__('base64').b64encode(path.read_bytes()).decode()+'">'
(ROOT/'reports/Delivery and acceptance.html').write_text('<!doctype html><meta charset="utf-8"><style>body{font:16px system-ui;max-width:1200px;margin:40px auto;padding:20px;color:#17313e}pre{white-space:pre-wrap;overflow-wrap:anywhere;line-height:1.65}a{color:#087d88}</style>'+body,encoding='utf-8')
print(json.dumps({'passed':sum(v['status']=='passed' for v in results.values()),'not_tested':sum(v['status']=='not_tested' for v in results.values()),'blocked':sum(v['status']=='blocked' for v in results.values())}))
