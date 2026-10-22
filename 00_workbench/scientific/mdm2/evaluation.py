"""Unblind only after frozen label-free ranking. Full-library denominators."""
import collections
import datetime
import shutil
from pathlib import Path
from rdkit import Chem
from scientific import artifacts as a
from .common import now,write_csv,command,ROOT
from .scoring import score,rank


def metrics_for(ranked,labels,total):
    active={r['parent_id'] for r in labels if r['label']=='active'};result={}
    for k in (10,20,50):
        hits=len(active & {r['parent_id'] for r in ranked[:k]})
        result[str(k)]={'active_hits':hits,'retrieved':min(k,len(ranked)),'Recall':hits/len(active) if active else None,
                        'EF':(hits/k)/(len(active)/total) if active and total else None,
                        'EF_denominator':'fixed k slots; unfilled slots are misses; frozen full-library prevalence'}
    return {'input_parent_count':total,'active_count':len(active),'eligible_ranked_count':len(ranked),'top_k':result}


def export_sdf(path,rows):
    with Chem.SDWriter(str(path)) as writer:
        for row in rows:
            mol=Chem.SDMolSupplier(row['pose_sdf'],removeHs=False)[int(row['pose_index'])]
            if mol is None:raise ValueError('Cannot export pose '+row['parent_id'])
            for name in list(mol.GetPropNames(includePrivate=True,includeComputed=False)):mol.ClearProp(name)
            mol.SetProp('_Name',row['parent_id'])
            for key in ('parent_id','state_id','pose_id','rank','profile','arm','vina_affinity','A1_recall','A2_recall'):
                mol.SetProp(key,str(row[key]))
            writer.write(mol)


def resource_summary(path,quota):
    rows=a.read(path,[])
    def usage(row):
        return int(dict(line.split() for line in row['cpu_stat'].splitlines())['usage_usec'])
    intervals=[]
    for left,right in zip(rows,rows[1:]):
        dt=right['monotonic']-left['monotonic']
        if dt>0 and left['cpu_stat'] and right['cpu_stat']:intervals.append((usage(right)-usage(left))/1e6/dt)
    gpu=[]
    for r in rows:
        try:gpu.append(float(r['gpu']['stdout'].split(',')[0]))
        except (ValueError,KeyError):pass
    return {'scope':'whole cgroup, not just Vina','samples':len(rows),
            'mean_cpu_cores_used':sum(intervals)/len(intervals) if intervals else None,
            'mean_cpu_quota_utilization_percent':100*sum(intervals)/len(intervals)/quota if intervals else None,
            'peak_memory_GiB':max([int(r['memory_current'])/2**30 for r in rows if r['memory_current']] or [0]),
            'mean_gpu_utilization_percent':sum(gpu)/len(gpu) if gpu else None}


def figures(root,result,labels):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for arm in ('A0','A1','A2'):
        values=result['profiles']['ppi_benchmark'][arm]['top_k']
        for ax,key in zip(axes,('Recall','EF')):ax.plot([10,20,50],[values[str(k)][key] for k in (10,20,50)],marker='o',label=arm);ax.set_ylabel(key);ax.set_xlabel('Top-k');ax.grid(alpha=.2)
    axes[0].set_ylim(0,1);axes[1].legend();fig.suptitle('Retrospective MDM2 pilot: ppi_benchmark');fig.tight_layout();fig.savefig(root/'recovery_metrics.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(11,6))
    for ax,key in zip(axes.flat,['MW','cLogP','TPSA','HBD','HBA','RotB']):
        for label in ('active','unknown_background'):
            ax.hist([r['properties'][key] for r in labels if r['label']==label],bins=15,density=True,alpha=.5,label=label)
        ax.set_title(key)
    axes.flat[0].legend(fontsize=8);fig.tight_layout();fig.savefig(root/'property_matching.png',dpi=160);plt.close(fig)


def evaluate(root):
    root=Path(root);frozen_path=root/'inputs/preregistered_protocol.json';frozen=a.read(frozen_path)
    if a.digest(frozen_path)!=frozen_path.with_suffix('.sha256').read_text().strip():raise ValueError('Protocol changed')
    for rel,digest in frozen['files'].items():
        if a.digest(root/rel)!=digest:raise ValueError('Frozen input changed: '+rel)
    for rel,digest in frozen['source_files'].items():
        if a.digest(ROOT/rel)!=digest:raise ValueError('Frozen source changed: '+rel)
    scores=score(root);filters={r['parent_id']:r for r in a.read(root/'inputs/filters.json')}
    # Ranking uses anonymous chemistry only. Labels are opened below, after all orders.
    orders={};ledgers=[];combined=[]
    for profile in ('current_product','ppi_benchmark'):
        eligible=[];failures={}
        for row in scores:
            f=filters[row['parent_id']]
            stages=[('parent_standardization',True),('state_3d_pdbqt',row['preparation_status']=='ok'),
                    ('property_filter',True if profile=='current_product' else f['ppi_pass']),
                    ('risk_filter',f['risk']['decision']=='pass' if profile=='current_product' else True),
                    ('acquisition_filter',f['acquisition']['decision']=='pass' if profile=='current_product' else True),
                    ('docking',row['docking_status'] in ('ok','failed_pose_stereo_qc')),('pose_stereo_qc',row.get('selected_pose_stereo_pass',False)),('peptide_feature_filter',True)]
            first=None
            for stage,passed in stages:
                if not passed and first is None:first=stage
                ledgers.append({'run_id':root.name,'profile':profile,'arm':'shared','seed':20261023,'parent_id':row['parent_id'],
                    'stage':stage,'passed':passed,'retained_after_stage':first is None,'first_exclusion_stage':first,
                    'risk_decision':f['risk']['decision'],'risk_alerts':f['risk']['alert_ids'],
                    'acquisition_status':f['acquisition']['status'],'error':row.get('error'),
                    'research_docking_shared_across_profiles':True})
            if first is None:eligible.append({**row,'profile':profile,'run_id':root.name,'seed':20261023,
                                             'deliverable_eligible':f['deliverable_pass'],'acquisition_status':f['acquisition']['status']})
            failures[row['parent_id']]=first
        for arm in ('A0','A1','A2'):
            ranked=rank(eligible,arm);orders[profile,arm]=ranked
            byid={r['parent_id']:r for r in ranked}
            for row in scores:
                combined.append(byid.get(row['parent_id'],{**row,'profile':profile,'arm':arm,'run_id':root.name,
                    'seed':20261023,'rank':None,'first_exclusion_stage':failures[row['parent_id']]}))
            path=root/'profiles'/profile/arm;path.mkdir(parents=True,exist_ok=True)
            columns=list(dict.fromkeys(k for r in combined for k in r))
            write_csv(path/'scientific_top50.csv',ranked[:50],columns)
            export_sdf(path/'scientific_top50.sdf',ranked[:50])
            deliver=rank([r for r in eligible if r['deliverable_eligible']],arm)
            write_csv(path/'deliverable_top50.csv',deliver[:50],columns);export_sdf(path/'deliverable_top50.sdf',deliver[:50])
    a.atomic(root/'anonymous_rankings.json',combined)
    # Seal rankings before opening evaluation labels.
    a.atomic(root/'evaluation/ranking_seal.json',{'sha256':a.digest(root/'anonymous_rankings.json'),'sealed_at':now()})
    if a.digest(root/'evaluation/labels.json')!=frozen['label_seal_sha256']:raise ValueError('Evaluation labels changed after freeze')
    labels=a.read(root/'evaluation/labels.json');labelby={r['parent_id']:r for r in labels}
    if set(labelby)!={r['parent_id'] for r in scores}:raise ValueError('Evaluation universe differs from frozen library')
    annotated=[{**r,'label':labelby[r['parent_id']]['label'],'name':labelby[r['parent_id']]['name'],
                'source_id':labelby[r['parent_id']]['source_id'],'scaffold':labelby[r['parent_id']]['scaffold']} for r in combined]
    write_csv(root/'all_rankings.csv',annotated)
    ledgers=[{**r,'label':labelby[r['parent_id']]['label'],'name':labelby[r['parent_id']]['name']} for r in ledgers]
    write_csv(root/'molecule_stage_ledger.csv',ledgers)
    recovery=[r for r in annotated if r['name'] in {'Nutlin-3a','SAR405838','Brigimadlin'}]
    write_csv(root/'reference_recovery.csv',recovery)
    result={'run_id':root.name,'seed':20261023,'evaluation_population':'entire frozen parent library; failures retained',
            'total_parents':len(labels),'actives':sum(r['label']=='active' for r in labels),
            'unknown_background':sum(r['label']=='unknown_background' for r in labels),'experimental_negatives':0,
            'ROC_AUC':None,'PR_AUC':None,'AUC_reason':'unknown backgrounds are not experimentally confirmed negatives',
            'docked':sum(r['docking_status']=='ok' for r in scores),
            'raw_vina_completed':sum(r.get('raw_vina_status')=='ok' for r in scores),
            'pose_stereo_qc_failures':sum(r['docking_status']=='failed_pose_stereo_qc' for r in scores),
            'study_design_status':'pilot: chemistry/QC protocol amended before label unblinding; not a fully prospective preregistration',
            'profiles':{},'scaffold_recall':[],
            'limitations':['17/20 actives from one medicinal chemistry series','single docking seed; no independent replication',
                'residue contacts are not a validated chemical pharmacophore','retrospective recovery, not new molecule discovery']}
    for profile in ('current_product','ppi_benchmark'):
        result['profiles'][profile]={arm:metrics_for(orders[profile,arm],labels,len(labels)) for arm in ('A0','A1','A2')}
        result['profiles'][profile]['differences']={f'{arm}_minus_A0':{str(k):result['profiles'][profile][arm]['top_k'][str(k)]['Recall']-result['profiles'][profile]['A0']['top_k'][str(k)]['Recall'] for k in (10,20,50)} for arm in ('A1','A2')}
        active_groups=collections.defaultdict(list)
        for label in labels:
            if label['label']=='active':active_groups[label['scaffold']].append(label['parent_id'])
        for arm in ('A0','A1','A2'):
            top={r['parent_id'] for r in orders[profile,arm][:50]}
            for scaffold,ids in active_groups.items():
                result['scaffold_recall'].append({'profile':profile,'arm':arm,'scaffold':scaffold,'active_count':len(ids),'hits50':len(top&set(ids)),'recall50':len(top&set(ids))/len(ids)})
    result['first_exclusion_counts']={}
    for profile in ('current_product','ppi_benchmark'):
        subset=[r for r in annotated if r['profile']==profile and r['arm']=='A1' and r['label']=='active']
        result['first_exclusion_counts'][profile]=dict(collections.Counter(r.get('first_exclusion_stage') or 'ranked' for r in subset))
    # Root scientific list is explicitly the ppi research arm A1; not a product claim.
    for name in ('scientific_top50.csv','scientific_top50.sdf','deliverable_top50.csv','deliverable_top50.sdf'):
        shutil.copyfile(root/'profiles/ppi_benchmark/A1'/name,root/name)
    shutil.copyfile(root/'inputs/blinded_library.sdf',root/'blinded_library.sdf')
    manifest=a.read(root/'run_manifest.json');finished=datetime.datetime.now(datetime.timezone.utc)
    elapsed=(finished-datetime.datetime.fromisoformat(manifest['started_at'])).total_seconds()/60
    result['protocol_amendment']=manifest.get('protocol_amendment')
    result['earlier_protocol_amendment']=manifest.get('earlier_protocol_amendment')
    result['reused_docking_jobs']=a.read(root/'reuse_manifest.json',{}).get('imported_jobs',0)
    result['docking_summary']=a.read(root/'docking_summary.json',{})
    result['elapsed_minutes']=elapsed;result['within_180_minutes']=elapsed<=180
    complete=result['docked']==result['total_parents'];result['complete_docking']=complete
    research=result['profiles']['ppi_benchmark'];a1=research['A1']['top_k']['50']['active_hits']
    ref={r['name']:r for r in recovery if r['profile']=='ppi_benchmark' and r['arm']=='A1'}
    result['engineering_acceptance']={'full_frozen_library_completed':complete,'all_requested_references_top50':all(ref[n]['rank'] is not None and ref[n]['rank']<=50 for n in ref),
        'at_least_one_active_top50':a1>0,'A1_point_estimate_exceeds_A0':a1>research['A0']['top_k']['50']['active_hits'],
        'A1_point_estimate_exceeds_A2':a1>research['A2']['top_k']['50']['active_hits'],'generalizable_peptide_gain_proven':False,
        'current_product_success':result['profiles']['current_product']['A1']['top_k']['50']['active_hits']>0}
    environment=a.read(root/'environment.json')
    resource_root=Path(manifest.get('compute_resource_root',str(root)))
    result['resource_evidence_root']=str(resource_root)
    result['resources']={phase:resource_summary(resource_root/phase/'resource_samples.json',environment['effective_cpu_quota']) for phase in ('development','inputs')}
    result['throughput_preflight']=a.read(root/'throughput.json')
    a.atomic(root/'metrics.json',result)
    figures(root,result,labels)
    manifest.update(status='completed' if complete else 'incomplete',finished_at=finished.isoformat(),elapsed_minutes=elapsed,
                    frozen_protocol_sha256=a.digest(frozen_path),git_commit=frozen.get('git_commit'),results_alias='root scientific_top50 = ppi_benchmark/A1; deliverable separately qualified')
    a.atomic(root/'run_manifest.json',manifest)
    lines=['# MDM2–p53 肽引导虚拟筛选回收验证报告','',f'本次因解盲前技术质检协议修订，标记为pilot。实际耗时 **{elapsed:.1f}分钟**。冻结库 **{len(labels)}** 个唯一母体，其中真实活性 **{result["actives"]}** 个，未知背景 **{result["unknown_background"]}** 个；实际Vina成功 **{result["raw_vina_completed"]}** 个，通过最佳姿势立体质检 **{result["docked"]}** 个；姿势质检失败 **{result["pose_stereo_qc_failures"]}** 个。',
        '', '这是回顾性混库回收验证，不是新分子发现。背景没有实验阴性标签；未计算ROC-AUC/PR-AUC。',
        '',f'协议修订记录：{result.get("protocol_amendment")}。前期版本在解盲前分别发现互变异构处理丢失E/Z构型和12个三维构象未保留指定构型，原始记录保留。正式版本保留来源互变异构形式、只作pH7.4处理，并从实际三维坐标反推构型核验；不合格构象统一采用ETKDGv3/MMFF94s重新生成，再次失败则保留失败记录。另外，统一剔除最佳对接姿势发生构型翻转的母体，保留其完整评价分母，不用次优姿势替代。该技术质检规则在标签解盲前冻结，但发生于原始对接之后，因此标记pilot，不宣称完全前瞻预注册。没有按活性标签改变处理。复用严格满足Vina输入字节和全部对接参数完全相同的{result["reused_docking_jobs"]}个既有结果，其余真实计算；未按参考排名调参。总耗时包含旧版审计与作废计算。',
        '', '## 参考分子结果','', '|配置|实验臂|分子|最终排名|首次淘汰阶段|','|---|---|---|---:|---|']
    for r in recovery:lines.append(f'|{r["profile"]}|{r["arm"]}|{r["name"]}|{r["rank"] if r["rank"] is not None else "无排名"}|{r.get("first_exclusion_stage") or "保留"}|')
    lines+=['','## 完整库指标','','|配置|臂|k|命中|Recall@k|EF@k|','|---|---|---:|---:|---:|---:|']
    for profile in ('current_product','ppi_benchmark'):
        for arm in ('A0','A1','A2'):
            for k,m in result['profiles'][profile][arm]['top_k'].items():lines.append(f'|{profile}|{arm}|{k}|{m["active_hits"]}|{m["Recall"]:.3f}|{m["EF"]:.3f}|')
    lines+=['','所有Recall分母均为冻结库全部活性分子，包括准备失败、被过滤或未完成者。EF=(命中数/k)/(全部活性数/全部母体数)，不足k的空名额视为未命中。',
        '', '## 两套配置的含义','',
        'current_product调用现有风险规则和可购买/合成证据门控。未提供有效供应目录及逆合成模型，不能生成产品合格结果。原项目MW 180–500为配置中的偏好，并非当前执行路径的硬过滤；因此不虚报Nutlin-3a被MW规则删除。',
        'ppi_benchmark采用预先固定的MW 150–800、cLogP -1–7.5、TPSA≤180、RotB≤18、HBD≤6、HBA≤15、电荷-2至+2；风险和可获取性仍记录。研究范围不是通用开发标准。根目录scientific_top50为ppi_benchmark/A1；deliverable仍必须通过真实产品证据门控。',
        '', '## 结构、盲法和计算','',
        '正式受体仅来自1YCR；按UniProt实体映射分离MDM2与p53，从肽接触重原子确定搜索盒后移除多肽。Phe19、Trp23、Leu26界面检查见target_audit.json。补建末端OXT，不改变已观测重原子坐标。5Z02、5TRF、8PWC及CCD仅用于独立身份核查和解释，其共晶坐标不用于配体3D或搜索盒。',
        '所有分子统一标准化、盐片段处理、保留来源互变异构形式、pH7.4酸碱状态枚举与固定规则选一状态、ETKDG/MMFF94s三维化、Meeko PDBQT。SAR405838的7HC描述缺少环己醇相对立体化学，已与原始Figure1a及CHEMBL2381408交叉核验。',
        '筛选代码只读匿名ID、结构及计算属性；排序后封存anonymous_rankings.json，再打开evaluation/labels.json。冻结文件、源代码与软件哈希见inputs/preregistered_protocol.json。',
        'Vina 1.2.7，exhaustiveness=8，num_modes=5，seed=20261023；5个worker×4 CPU。GPU不参与标准Vina。GPU遥测可能显示100%利用率但0 MiB显存且无可见GPU进程，该异常不能归因于本次CPU计算，独立记录见evaluation/gpu_telemetry_audit.json。A0按Vina排序，A1按既有多肽残基接触召回再按Vina排序，A2对接触原子作限定刚体扰动后同样排序；三组复用同一最佳Vina姿势、同一计算预算和Top-k名额。',
        'A2保持特征数量、元素类型及内部距离，位置在盒内并排除严重重叠；参数由几何条件确定，未按测试排名调参。同时导出RDKit原子化学规则识别的芳香、疏水、供受体等药效团描述及其A2刚体扰动坐标。描述用于界面解释，排名仍为现有距离接触指纹；不能宣称完整化学药效团模型已获性能验证。',
        '',f'资源实测（整个容器配额）：`{result["resources"]}`。50分子预检使用修订前状态策略，仅用于CPU吞吐估计；耗时{result["throughput_preflight"]["elapsed_seconds"]:.1f}秒，平均单分子{result["throughput_preflight"]["mean_job_seconds"]:.1f}秒；失败率{result["throughput_preflight"]["failure_rate_including_preparation"]:.1%}。',
        '', '![回收指标](recovery_metrics.png)', '', '![性质匹配](property_matching.png)',
        '', '## 实验依据与边界','',
        '17个稳定D系列分子：原始J Med Chem 2014 Table6，DOI 10.1021/jm501541j，human MDM2(1–118)竞争FP Ki；逐条核对ChEMBL3380268与compound_record，排除不稳定C异构体。',
        'Nutlin-3a Ki=43.6 nM、SAR405838 Ki=0.88 nM：Cancer Res原始Figure1b，DOI 10.1158/0008-5472.CAN-14-0799。Brigimadlin生化IC50=2 nM：原始Supplementary TableS4及Methods p2，DOI 10.1158/1535-7163.MCT-23-0783；GST-MDM2(23–117)，不是细胞实验。',
        '17/20活性来自同一研究系列，因此属于限定化学空间的回顾性验证；不能据此推断广泛骨架泛化。性质匹配诊断见evaluation/matching_diagnostics.json；匹配并非证明背景无活性。单种子仅支持本次点估计比较，不能证明稳定、可泛化的肽引导增益。',
        '', '## 验收与问题','',f'工程验收记录：`{result["engineering_acceptance"]}`。',
        '已实际运行：结构实体解析、肽界面提取、匿名库准备、50分子吞吐预检、真实CPU Vina、共享姿势三组重排、两套规则审计与完整分母评价。是否成功回收和是否优于对照，应以上表为准。',
        '尚未验证：独立系列/多种子泛化、实验命中率、完整药效团功能等价、真实采购或合成可行性。未运行MD、MM/PBSA或FEP，未改变现有MD默认值。',
        '各阶段活性淘汰数见metrics.json，逐分子路径见molecule_stage_ledger.csv，骨架分层见scaffold_recall.csv。原始姿势、命令、失败记录及哈希均保留，可断点续跑，不重复已成功且身份一致的计算。',
        '', '## 复现','',f'Git commit: `{frozen.get("git_commit")}`。配置冻结时间：{frozen["frozen_at"]}。',
        '在同一安装环境依次运行`python yuanhe benchmark run --root RUN_ROOT`和`python yuanhe benchmark evaluate --root RUN_ROOT`；冻结哈希不符将拒绝执行。复制到其他目录需保存原始证据并创建明确的新运行清单，不能静默改写冻结路径。',
        '原始结构：https://www.rcsb.org/structure/1YCR；参考5Z02/5TRF/8PWC。原始文献和来源哈希位于evaluation/literature、assay_evidence.csv。', '']
    (root/'validation_report.md').write_text('\n'.join(lines))
    write_csv(root/'scaffold_recall.csv',result['scaffold_recall'])
    return result
