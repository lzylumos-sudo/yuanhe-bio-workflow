"""Single-device generation leases, with explicit CPU fallback reporting."""
from .resources import detect,gpu_lease,disk_guard

def sample_reinvent(adapter,kwargs,path,cfg):
    from webui.workflow import atomic
    from yuanhe.runtime import runtime_env
    disk_guard(path,cfg)
    environment=runtime_env(kwargs['executable'],for_overlay=True)
    available=[g for g in detect(path)['gpus'] if g['allocatable']]
    if cfg['reinvent_device']!='cpu':
        for gpu in available:
            try:
                with gpu_lease(gpu):
                    return adapter(**kwargs,device='cuda:0',environment={**environment,'CUDA_VISIBLE_DEVICES':gpu['uuid']})
            except BlockingIOError:continue
        if cfg['reinvent_device']=='cuda':raise RuntimeError('No authorized free GPU for REINVENT4')
    atomic(path/'reinvent_cpu_fallback.json',{'reason':'CPU explicitly selected or no authorized unallocated GPU','device':'cpu'})
    return adapter(**kwargs,device='cpu',environment={**environment,'CUDA_VISIBLE_DEVICES':''})

def sample_peptide(adapter,kwargs,path):
    for gpu in detect(path)['gpus']:
        if not gpu['allocatable']:continue
        try:
            with gpu_lease(gpu):return adapter(**kwargs,environment={'CUDA_VISIBLE_DEVICES':gpu['uuid']})
        except BlockingIOError:continue
    raise RuntimeError('No authorized free GPU for Peptide2Mol; no uncoordinated GPU launch')

def guarded_sample(adapter,kwargs,path,cfg,kind):
    """Completed batches are reused; uncertain partial sampling never silently restarts."""
    from pathlib import Path
    from webui.workflow import atomic,read,digest
    from .config import fingerprint
    out=Path(kwargs['output_dir']);ledger=path/'generation_ledger';ledger.mkdir(exist_ok=True)
    signature={k:str(v) for k,v in kwargs.items() if k!='output_dir'}
    key=fingerprint({'kind':kind,'inputs':signature});record=ledger/(key+'.json');old=read(record)
    if old:
        if old['status']!='completed':raise RuntimeError('Interrupted batch reservation retained; raw-attempt accounting must be reviewed before additional sampling')
        if any(not (out/name).is_file() or digest(out/name)!=sha for name,sha in old['files'].items()):raise ValueError('Generation cache integrity failed')
        return read(out/'manifest.json')
    atomic(record,{'status':'reserved','budget':kwargs['budget'],'inputs':signature})
    result=sample_reinvent(adapter,kwargs,path,cfg) if kind=='reinvent4' else sample_peptide(adapter,kwargs,path)
    if result['report']['completed']>kwargs['budget']:raise ValueError('Generator exceeded reserved raw sampling budget')
    atomic(record,{'status':'completed','budget':kwargs['budget'],'inputs':signature,
                   'files':{str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file()}})
    return result
