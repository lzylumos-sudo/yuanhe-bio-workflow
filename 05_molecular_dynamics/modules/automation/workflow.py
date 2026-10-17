"""Unified local input/decision/state interface for the preserved OpenMM backend."""
from __future__ import annotations
import argparse, contextlib, copy, csv, importlib.metadata, io, json, os, platform
import re, shutil, subprocess, sys, traceback, uuid
from datetime import datetime, timezone
from pathlib import Path
from . import __version__
from .state import digest, fingerprint, read, write, update, files_under, verify, job_lock, now, LockBusy
from .input_resolver import resolve_input, SelectionRequired

ROOT=Path(__file__).resolve().parents[2]

def load_mapping(path):
    import yaml
    data=yaml.safe_load(Path(path).read_text(encoding='utf-8-sig'))
    if not isinstance(data,dict): raise ValueError('配置/决定文件必须是键值映射。')
    return data

def merge(a,b):
    for k,v in b.items():
        if isinstance(v,dict) and isinstance(a.get(k),dict): merge(a[k],v)
        else: a[k]=copy.deepcopy(v)
    return a

def implementation():
    paths=list((ROOT/'modules').rglob('*.py'))+[ROOT/'scripts/run_workflow.py']
    return {p.relative_to(ROOT).as_posix():digest(p) for p in paths}

def settings(args, interactive):
    cfg=load_mapping(ROOT/'config/default.yaml')
    merge(cfg,{'MD':{'production':'100ps','cpu_threads':2},'replicates':1})
    if args.config: merge(cfg,load_mapping(Path(args.config).expanduser().resolve()))
    decisions=cfg.pop('decisions',{})
    if not isinstance(decisions,dict): raise ValueError('decisions 必须是映射。')
    if args.smoke_test:
        merge(cfg,{'MD':{'platform':'CPU','cpu_threads':2,'production':'10ps','equilibration':'2ps',
             'nvt':'1ps','npt':'1ps','report_interval':'1ps','checkpoint_interval':'2ps'},
             'replicates':1,'analysis':{'stride':1,'discard_ns':0,'sasa_stride':1},'mmgbsa':{'max_frames':4,'start_ns':0}})
    elif interactive and not args.config:
        from .forcefields import list_presets
        print(json.dumps(list_presets(),ensure_ascii=False,indent=2))
        print('当前启用整套预设：AMBER ff14SB + GAFF 2.11/AM1-BCC + TIP3P/NaCl。其他组合待验证。默认是短预运行。')
        cfg['MD']['production']=input('生产时长 [100ps]：').strip() or '100ps'
        cfg['replicates']=int(input('独立重复数 [1]：').strip() or '1')
        cfg['MD']['cpu_threads']=int(input('CPU 线程数 [2]：').strip() or '2')
        cfg['MD']['temperature']=float(input('温度 K [300]：').strip() or '300')
        cfg['protein']['ph']=float(input('蛋白默认加氢 pH（已明确的 HIS 变体优先）[7]：').strip() or '7')
        cfg['MD']['ionic_strength_molar']=float(input('NaCl 浓度 mol/L [0.15]：').strip() or '0.15')
        cfg['mmgbsa']['salt_molar']=cfg['MD']['ionic_strength_molar']
    cfg.update(workflow='md_from_complex_v2',configuration_version=2,project_name='Nexus_MD_v2')
    return cfg,decisions

def new_job(source,cfg,decisions,args):
    sha=digest(source)
    jobid=args.run_id or datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ_')+sha[:8]+'_'+uuid.uuid4().hex[:6]
    if not re.fullmatch('[A-Za-z0-9_-]+',jobid): raise ValueError('任务 ID 只使用字母、数字、下划线和连字符。')
    job=ROOT/'runs'/jobid
    job.mkdir(parents=True,exist_ok=False)
    raw=job/'raw';raw.mkdir(); copy=raw/source.name;shutil.copyfile(source,copy)
    if digest(copy)!=sha: raise ValueError('输入复制校验失败。')
    templates={};embedded_manifest={}
    from modules.utils.embedded_chemistry import read_embedded
    embedded=read_embedded(source.read_text())
    side_manifest={}
    for name in ['manifest.json','complex_manifest.json','manifest.yaml','complex_manifest.yaml']:
        if (source.parent/name).is_file(): side_manifest=load_mapping(source.parent/name);break
    if embedded:
        embedded_manifest=embedded['manifest']
        if side_manifest and side_manifest!=embedded_manifest: raise ValueError('内嵌与旁置 manifest 冲突。')
        if embedded.get('ligand_sdf'):
            (raw/'embedded-ligand.sdf').write_text(embedded['ligand_sdf'])
            templates['default']=str((raw/'embedded-ligand.sdf').relative_to(job))
    ligand_meta=(embedded_manifest or side_manifest).get('ligand') or {}
    name=Path(ligand_meta.get('template_sdf','ligand.sdf'))
    if name.is_absolute() or '..' in name.parts: raise ValueError('旁置模板路径不允许越过输入文件目录。')
    side=(source.parent/name).resolve()
    if not side.is_relative_to(source.parent): raise ValueError('旁置模板链接越过输入目录。')
    if side.is_file():
        if embedded and embedded.get('ligand_sdf') and side.read_text()!=embedded['ligand_sdf']:
            raise ValueError('内嵌与旁置 SDF 冲突。')
        dest=raw/'provided-ligand.sdf';shutil.copyfile(side,dest);templates['default']=str(dest.relative_to(job))
    meta=embedded_manifest or side_manifest
    if meta: write(raw/'provided-manifest.json',meta)
    write(job/'source.json',{'original_filename':source.name,'source_path':str(source),'sha256':sha,
          'staged':copy.relative_to(job).as_posix(),'job_id':jobid,'templates':templates,
          'provided_manifest':meta,'frontend_version':__version__,'configuration_version':2})
    write(job/'settings.json',cfg);write(job/'decisions.json',decisions)
    write(job/'raw_hashes.json',files_under(job,['raw']))
    update(job,'NEW',job_id=jobid,source_sha256=sha,scope='engineering_smoke' if args.smoke_test else 'user_configured')
    return job

def waiting(job,pending,status='WAITING_FOR_DECISION'):
    write(job/'pending.json',{'status':status,'items':pending,
         'resume_command':f'bash run_md.sh --resume runs/{job.name} --decisions /path/to/decisions.yaml',
         'note':'原始输入和已完成阶段保留；不会等待无人值守的键盘输入。'})
    update(job,status,pending_count=len(pending))
    print(json.dumps({'status':status,'job':str(job),'pending':pending},ensure_ascii=False,indent=2))
    return 4 if status=='DEPENDENCY_BLOCKED' else 3

def selected_components(inv,normalized):
    return normalized['components']

def relocate(value, job, restore=False):
    """Keep runtime artifact references portable; preserve source provenance separately."""
    if isinstance(value,dict):return {k:relocate(v,job,restore) for k,v in value.items()}
    if isinstance(value,list):return [relocate(v,job,restore) for v in value]
    if isinstance(value,str):
        if restore and value.startswith('@job/'):
            path=(job/value[5:]).resolve()
            if not path.is_relative_to(job.resolve()):raise ValueError('任务路径越界。')
            return str(path)
        if not restore and value.startswith(str(job)+'/'):
            return '@job/'+value[len(str(job))+1:]
    return value

def map_analysis_groups(groups, residues):
    chains={}
    for residue in residues:
        if residue['kind']=='protein':
            chains.setdefault(residue['original']['chain'],set()).add(residue['chain'])
    result=copy.deepcopy(groups)
    for role,definition in result.items():
        if role in ('ligand','glue'):continue
        if not isinstance(definition,dict) or not isinstance(definition.get('chains'),list):
            raise ValueError('分析组必须包含 chains 列表：'+str(role))
        mapped=[]
        for chain in definition['chains']:
            if chain not in chains or len(chains[chain])!=1:
                raise ValueError('原始链标识无法唯一映射到分析组：'+repr(chain))
            mapped.extend(chains[chain])
        definition['chains']=mapped
    return result

def preflight(job,cfg,decisions,interactive):
    from .pdb_inventory import inventory_pdb, normalize_pdb, InventoryDecisionRequired
    from .ligand_chemistry import recover_component
    source=read(job/'source.json'); raw=job/source['staged']
    inv=inventory_pdb(raw);write(job/'inventory.json',inv)
    if interactive:
        models=sorted({a['model'] for a in inv['atoms']})
        if len(models)>1 and 'model' not in decisions:
            decisions['model']=int(input('选择 MODEL '+str(models)+'：'))
        for item in inv.get('required_decisions',[]):
            if item.get('type')=='altloc' and 'altloc' not in decisions:
                value=input('存在替代构象，请输入统一 altLoc 标签（如 A）；不能混合不兼容构象：').strip()
                if value:decisions['altloc']=value
        if 'keep_components' not in decisions and (sum(c['kind']=='ligand' for c in inv['components'])>1 or any(c['kind'] not in ('protein','ligand','water') for c in inv['components'])):
            print(json.dumps([{k:c.get(k) for k in ('id','kind','atom_count','residue_ids')} for c in inv['components']],ensure_ascii=False,indent=2))
            choice=input('输入需要保留的组分 ID，用逗号分隔；回车全部保留并检查支持范围：').strip()
            if choice:decisions['keep_components']=[s.strip() for s in choice.split(',')]
    # Each decision attempt preserves its normalization artifacts for review.
    attempt=job/'normalized'/('attempt_'+uuid.uuid4().hex[:12])
    try: normalized=normalize_pdb(inv,decisions,attempt)
    except InventoryDecisionRequired as exc:
        write(job/'decisions.json',decisions)
        return None,waiting(job,exc.decisions)
    write(job/'normalization.json',normalized)
    comps=selected_components(inv,normalized)
    ligands=[c for c in comps if c['kind']=='ligand']
    protein=[c for c in comps if c['kind']=='protein']
    pending=[]
    if not protein:pending.append({'id':'protein','message':'未识别到受支持蛋白，不能构建普通蛋白复合物。'})
    others=[c for c in comps if c['kind'] not in ('protein','ligand','water')]
    if others:pending.append({'id':'unsupported_components','components':[{k:c.get(k) for k in ['id','kind','residue_ids']} for c in others],
          'message':'保留了离子/金属/修饰组分的原始记录；当前普通非共价后端没有该组分的已验证参数路径。明确选择构建体后才可排除。'})
    waters=[c for c in comps if c['kind']=='water']
    if waters and decisions.get('water_policy')!='replace_bulk':
        if interactive and input('输入含水；当前预设将移除输入水并重新溶剂化。确认不是必须保留的结合位点水？输入 yes：').strip().lower()=='yes':
            decisions['water_policy']='replace_bulk'
        else:pending.append({'id':'water_policy','key':'water_policy','options':['replace_bulk'],
                    'message':'当前后端只支持重新溶剂化；有功能性晶体水时不要确认移除。'})
    if len(ligands)>1:
        pending.append({'id':'multiple_ligands','components':[c['id'] for c in ligands],
          'message':'已识别多个配体/拷贝，未自动丢弃任何一个。当前模拟后端只验证了单配体；如研究体系确实只需一个，明确设置 keep_components，否则需要多配体后端。'})
    if not inv.get('sequence_verified',False) and decisions.get('construct_confirmed') is not True:
        if interactive and input('尚未验证完整目标序列。确认已有蛋白链就是目标构建体，内部缺口已经核对？输入 yes：').strip().lower()=='yes':
            decisions['construct_confirmed']=True
        else:pending.append({'id':'construct_confirmed','key':'construct_confirmed','options':[True],
              'message':'尚未验证完整目标序列，不能证明无缺失残基；确认实际构建体。不会自动补造起始编号之前的残基。'})
    # Chemical input anomalies are never repaired during lossless normalization.
    for item in normalized.get('pending_decisions',[]):
        kind=str(item.get('type','')).lower(); ident=item.get('id','')
        if any(x in kind for x in ['termin','hxt']):
            termini=decisions.get('termini',{})
            selected=termini.get(item.get('residue_id')) if isinstance(termini,dict) else termini
            if selected!='standard_charged':
                if interactive and input(str(item.get('message','末端化学异常'))+'\n选择常规带电末端（可能删除HXT并补OXT），输入 standard_charged；否则回车保留待决策：').strip()=='standard_charged':
                    decisions.setdefault('termini',{})[item['residue_id']]='standard_charged'
                else:pending.append({**item,'key':'termini.'+item['residue_id'],'options':['standard_charged'],
                     'message':str(item.get('message',''))+' 需选择末端处理；封端或特殊末端当前没有已验证路径。'})
        else:
            pending.append(item)
    chemistry={}; target=None
    if len(ligands)==1:
        component=copy.deepcopy(ligands[0]);target=component['id']
        if not isinstance(component.get('atoms'),list):
            ids=set(component['atom_ids']);component['atoms']=[a for a in inv['atoms'] if a['id'] in ids]
        if 'bonds' not in component:
            ids=set(component['atom_ids']);component['bonds']=[b for b in inv['bonds'] if set(b['atom_ids'])<=ids]
        # Only selected MODEL/altLoc atoms may contribute chemistry evidence.
        if normalized.get('selected_atom_ids'):
            ids=set(normalized['selected_atom_ids'])
            component['atoms']=[a for a in component['atoms'] if a['id'] in ids]
            component['atom_ids']=[a['id'] for a in component['atoms']]
            allowed=set(component['atom_ids']);component['bonds']=[b for b in component['bonds'] if set(b['atom_ids'])<=allowed]
        chem_decision=decisions.get('chemistry',{}).get(target,{})
        template=None
        if source['templates'].get('default'):
            template={'path':str(job/source['templates']['default']),'trusted':True}
            mapping=source.get('provided_manifest',{}).get('ligand',{}).get('atom_mapping',{})
            if mapping:
                original={a['id']:a for a in inv['atoms']}
                mapped={}
                for atom in component['atoms']:
                    name=original[atom['id']]['atom_name']
                    if name in mapping:mapped[atom['id']]=int(mapping[name])-1
                if len(set(mapped.values()))!=len(mapped):
                    raise ValueError('模板 atom_mapping 在重复原子名下不唯一；使用 stable ID 映射或确认候选。')
                template['atom_map']=mapped
        if chem_decision.get('smiles'):template={'smiles':chem_decision['smiles'],'trusted':True}
        if chem_decision.get('template_sdf'):
            external=Path(chem_decision['template_sdf']).expanduser().resolve()
            local=job/'chemical_definitions'/(target+'.sdf');local.parent.mkdir(exist_ok=True)
            if not local.exists() or digest(local)!=digest(external):shutil.copyfile(external,local)
            template=str(local)
        chemistry=recover_component(component,job/'chemistry'/target,decision=chem_decision,template=template,timeout=30)
        if chemistry.get('status')=='WAITING' and interactive and chemistry.get('candidate_sha256'):
            print(json.dumps(chemistry,ensure_ascii=False,indent=2))
            if input('核对候选的键级、手性、氢与电荷；明确接受这一个候选？输入 yes：').strip().lower()=='yes':
                chem_decision={**chem_decision,'accept_candidate_sha256':chemistry['candidate_sha256']}
                decisions.setdefault('chemistry',{})[target]=chem_decision
                chemistry=recover_component(component,job/'chemistry'/target,decision=chem_decision,template=template,timeout=30)
        if chemistry.get('status')!='READY':pending.extend(chemistry.get('required_decisions') or [{'id':'chemistry','component':target,**chemistry}])
        write(job/'chemistry_result.json',chemistry)
    write(job/'decisions.json',decisions)
    if pending:
        return None,waiting(job,pending,'DEPENDENCY_BLOCKED' if chemistry.get('status')=='DEPENDENCY_BLOCKED' and len(pending)==1 else 'WAITING_FOR_DECISION')
    plan={**normalized,'decisions':decisions,'sequence_verified':bool(inv.get('has_seqres')),
          'inventory_path':str(job/'inventory.json'),'normalized':normalized,
          'chemistry':chemistry,'ligand_component':target,'source_sha256':source['sha256']}
    if target:
        plan['ligand_sdf']=chemistry.get('sdf_path') or chemistry.get('sdf') or str(job/'chemistry'/target/'ligand.sdf')
    plan['sequence_verified']=bool(inv.get('sequence_verified',False))
    write(job/'plan.json',relocate(plan,job))
    return plan,None

class Tee:
    def __init__(self,*handles):self.handles=handles
    def write(self,s):
        for h in self.handles:h.write(s);h.flush()
        return len(s)
    def flush(self):
        for h in self.handles:h.flush()

def execute(job,args,interactive):
    source=read(job/'source.json');cfg=read(job/'settings.json');decisions=read(job/'decisions.json')
    verify(job,read(job/'raw_hashes.json'))
    if args.decisions:merge(decisions,load_mapping(Path(args.decisions).expanduser().resolve()))
    if args.config and args.resume:raise ValueError('不能在续跑时替换配置。为新设置建立新任务；不会复用旧拓扑。')
    current=read(job/'state.json')
    locked=(job/'execution_fingerprint.json').exists()
    if locked:
        frozen=read(job/'execution_fingerprint.json')
        if frozen['settings_sha256']!=fingerprint(cfg) or frozen['decisions_sha256']!=fingerprint(decisions):
            raise ValueError('上游结构/科学决定变化，下游结果已失去可复用资格。请建立新任务，旧结果完整保留。')
        if frozen['implementation']!=implementation():raise ValueError('实现变化，拒绝静默续跑。请用原版本恢复或建立新任务。')
        verify(job,frozen['preflight_files'])
        plan=relocate(read(job/'plan.json'),job,True);cfg=read(job/'effective_config.json')
    else:
        update(job,'RUNNING',active_stage='input_preflight')
        plan,code=preflight(job,cfg,decisions,interactive)
        if code is not None:return code
        cfg.update(complex={'file':str((job/source['staged']).relative_to(ROOT))},ligand=None)
        if plan.get('ligand_sdf'):
            sdf=Path(plan['ligand_sdf']).resolve()
            cfg['ligand']={'resname':'LIG','template_sdf':sdf.relative_to(ROOT).as_posix()}
        chains={r['chain'] for r in plan.get('residues',[]) if r.get('kind','protein')=='protein'}
        cfg['system_type']='ternary' if cfg['ligand'] and len(chains)>1 else ('protein_ligand' if cfg['ligand'] else 'protein_protein')
        # Upstream roles refer to original chains; normalize chain identities first.
        provided=source.get('provided_manifest',{}).get('analysis',{}).get('groups',{})
        if provided and not cfg['analysis'].get('groups'):
            cfg['analysis']['groups']=provided
        cfg['analysis']['groups']=map_analysis_groups(cfg['analysis'].get('groups',{}),plan.get('residues',[]))
        from .backend import validate_config
        validate_config(cfg)
        write(job/'effective_config.json',cfg)
        import yaml
        (job/'config.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False,allow_unicode=True))
    print('任务：'+str(job));print('实际设置：'+json.dumps({'MD':cfg['MD'],'replicates':cfg['replicates']},ensure_ascii=False))
    if args.inspect_only:
        if not locked:update(job,'READY',active_stage='input_preflight',note='仅输入准备检查；未参数化或运行 MD。')
        return 0
    from .backend import dependency_report,prepare,build,run,analyze,audit_run
    dependencies=dependency_report(cfg);write(job/'dependencies.json',dependencies)
    # Dependency report is expected to advertise a boolean ready flag.
    ready=dependencies.get('ready',dependencies.get('available',dependencies.get('status') in ['READY','PASS']))
    if not ready:return waiting(job,[{'id':'dependencies','report':dependencies}],'DEPENDENCY_BLOCKED')
    if locked and frozen['dependencies'].get('compatibility_fingerprint')!=dependencies.get('compatibility_fingerprint'):
        raise ValueError('运行依赖版本或参数文件变化，拒绝复用检查点。请恢复原环境或建立新任务。')
    if platform.system()!='Linux' or platform.machine() not in ('x86_64','AMD64'):
        return waiting(job,[{'id':'target_platform','message':'原生动态验证仅在 x86-64 Ubuntu/Linux 执行；Mac 输入测试不是 Linux MD 验证。'}],'DEPENDENCY_BLOCKED')
    free=shutil.disk_usage(job).free
    if free<2*1024**3:
        return waiting(job,[{'id':'free_space','free_bytes':free,'message':'可用空间低于 2 GiB，停止建系。正式长模拟需要按体系和输出频率另行估算。'}],'DEPENDENCY_BLOCKED')
    if not locked:
        write(job/'execution_fingerprint.json',{'settings_sha256':fingerprint(read(job/'settings.json')),
             'decisions_sha256':fingerprint(decisions),'implementation':implementation(),
             'dependencies':dependencies,'preflight_files':files_under(job,['raw','normalized','chemistry','plan.json','config.yaml','effective_config.json'])})
    paths={k:job/k for k in ['md','analysis','reports']}
    for p in paths.values():p.mkdir(parents=True,exist_ok=True)
    (job/'logs').mkdir(exist_ok=True)
    from .analysis_mapping import add_original_identifiers
    def mapped_build():
        build(cfg,paths)
        add_original_identifiers(job,plan,paths)
    def mapped_analysis():
        analyze(cfg,paths)
        add_original_identifiers(job,plan,paths)
    phases=[('prepare',lambda:prepare(cfg,paths,plan),['md/prepared']),
            ('build',mapped_build,['md/system','mappings/topology_atom_map.json']),
            ('md',lambda:run(cfg,paths),['md']),
            ('analysis',mapped_analysis,['analysis','mappings/analysis_original_identifiers.json']),
            ('audit',lambda:audit_run(cfg,paths,reload_checkpoint=True),['reports','acceptance.json'])]
    for name,action,folders in phases:
        state=read(job/'state.json');stages=state.get('stages',{});record=stages.get(name,{})
        if record.get('status')=='COMPLETE': verify(job,record['files'])
        else:
            stages[name]={'status':'RUNNING','started':now()};update(job,'RUNNING',stages=stages,active_stage=name)
            with (job/'logs'/(name+'.log')).open('a',encoding='utf-8') as logfile:
                with contextlib.redirect_stdout(Tee(sys.stdout,logfile)),contextlib.redirect_stderr(Tee(sys.stderr,logfile)):
                    try:result=action()
                    except BaseException:
                        traceback.print_exc();raise
            if name=='audit':
                if not result.get('passed'):raise ValueError('动态验收未通过，不能标记 COMPLETE。')
                write(job/'acceptance.json',result)
            stages[name].update(status='COMPLETE',completed=now(),files=files_under(job,folders))
            update(job,'RUNNING',stages=stages,active_stage=name)
        if args.prepare_only and name=='build':
            if current['status']!='COMPLETE':update(job,'PREPARED',active_stage='build')
            print('结构、参数化与体系构建已完成；本次未推进模拟。');return 0
    # Original identifiers remain available beside the analyzed topology.
    update(job,'COMPLETE',active_stage=None,completed=now())
    (job/'REPORT.txt').write_text('任务完成：'+job.name+'\n配置：config.yaml\n验收：acceptance.json\n原始编号映射：analysis/input_atom_map.json、input_residue_map.json\n短测只验证工程链路，不证明充分平衡、收敛、结合亲和力或协同性。\n')
    (job/'pending.json').unlink(missing_ok=True)
    print(json.dumps({'status':'COMPLETE','job':str(job),'acceptance':str(job/'acceptance.json')},ensure_ascii=False))
    return 0

def main(argv=None):
    p=argparse.ArgumentParser(description='Nexus MD：无损输入、显式科学决定、OpenMM 自动运行')
    p.add_argument('--input');p.add_argument('--input-dir',default=str(ROOT/'inputs'))
    p.add_argument('--config');p.add_argument('--decisions');p.add_argument('--run-id')
    p.add_argument('--resume');p.add_argument('--non-interactive',action='store_true')
    mode=p.add_mutually_exclusive_group();mode.add_argument('--prepare-only',action='store_true');mode.add_argument('--smoke-test',action='store_true');mode.add_argument('--inspect-only',action='store_true')
    p.add_argument('--list-presets',action='store_true')
    args=p.parse_args(argv);interactive=sys.stdin.isatty() and not args.non_interactive
    job=None
    try:
        if args.list_presets:
            from .forcefields import list_presets
            print(json.dumps(list_presets(),ensure_ascii=False,indent=2));return 0
        if args.resume:
            if args.input or args.run_id or args.smoke_test:raise ValueError('续跑不接受新输入、任务 ID 或新 smoke 设置。')
            value=Path(args.resume).expanduser()
            job=(value if value.is_absolute() else ROOT/value).resolve()
            if not job.is_relative_to(ROOT/'runs') or not (job/'source.json').is_file():
                raise ValueError('新入口续跑路径必须为本项目 runs/<任务ID>；旧 results 任务继续使用 run.sh --run-id ID --resume。')
        else:
            source=resolve_input(args.input,args.input_dir,interactive)
            cfg,decisions=settings(args,interactive)
            if args.decisions:merge(decisions,load_mapping(Path(args.decisions).expanduser().resolve()))
            job=new_job(source,cfg,decisions,args)
        with job_lock(job):return execute(job,args,interactive)
    except SelectionRequired as exc:
        job=ROOT/'runs'/('selection_'+uuid.uuid4().hex[:12]);job.mkdir(parents=True)
        return waiting(job,[{'id':'input','message':str(exc),'candidates':exc.candidates}])
    except LockBusy as exc:
        print('ERROR: '+str(exc),file=sys.stderr);return 5
    except KeyboardInterrupt:
        if job and (job/'state.json').exists():update(job,'INTERRUPTED',message='用户中断；后端从最后已提交检查点续跑。')
        return 130
    except Exception as exc:
        if job and (job/'state.json').exists():update(job,'FAILED',error_type=type(exc).__name__,message=str(exc))
        print('ERROR: '+str(exc),file=sys.stderr)
        if job: print('任务记录：'+str(job),file=sys.stderr)
        return 2
