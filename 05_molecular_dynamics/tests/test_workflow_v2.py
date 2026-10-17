"""Portable CLI/state regression. No mock outcome is native MD validation."""
import contextlib
import io
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from modules.automation import workflow as w
from modules.automation.input_resolver import candidates,resolve_input,SelectionRequired
from modules.automation.state import write,read,verify,digest,job_lock,LockBusy


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='Nexus Path (test).')
        self.root=Path(self.tmp.name).resolve()
        shutil.copytree(ROOT/'config',self.root/'config')
        self.override=patch.object(w,'ROOT',self.root);self.override.start()
    def tearDown(self):
        self.override.stop();self.tmp.cleanup()
    def call(self,*args):
        with contextlib.redirect_stdout(io.StringIO()),contextlib.redirect_stderr(io.StringIO()):
            return w.main([*args,'--non-interactive'])
    def fixture(self,name='1fkj'):
        dst=self.root/'my_complex (Input)';shutil.copytree(ROOT/'examples'/name,dst)
        source=dst/'complex.pdb';dest=dst/'multiple_points.Structure 1.PDB';source.rename(dest)
        return dest
    def decisions(self,**values):
        file=self.root/'my-decisions.json';write(file,values);return str(file)
    def test_discovery_explicit_case_unicode_and_no_history(self):
        d=self.root/'inputs';d.mkdir();a=d/'Complex (A).v2.PDB';a.write_text('END\n')
        (d/'runs').mkdir();(d/'runs/history.pdb').write_text('END\n')
        self.assertEqual(resolve_input(None,d),a)
        b=d/'second.pdb';b.write_text('END\n')
        with self.assertRaises(SelectionRequired):resolve_input(None,d)
        self.assertEqual(resolve_input(str(a),d),a)
    def test_noninteractive_multiple_outputs_selection_without_hang(self):
        d=self.root/'inputs';d.mkdir()
        for name in ['a.pdb','b.PDB']:(d/name).write_text('END\n')
        self.assertEqual(self.call('--input-dir',str(d),'--inspect-only'),3)
        report=read(next((self.root/'runs').glob('*/pending.json')))
        self.assertEqual(len(report['items'][0]['candidates']),2)
    def test_source_hash_and_unique_jobs_no_overwrite(self):
        source=self.fixture();before=digest(source)
        args=['--input',str(source),'--inspect-only','--decisions',self.decisions(construct_confirmed=True)]
        self.assertEqual(self.call(*args),0);self.assertEqual(self.call(*args),0)
        self.assertEqual(len(list((self.root/'runs').iterdir())),2)
        self.assertEqual(digest(source),before)
        for job in (self.root/'runs').iterdir():
            provenance=read(job/'source.json')
            self.assertEqual(provenance['original_filename'],source.name)
            self.assertEqual(digest(job/provenance['staged']),before)
            self.assertEqual(read(job/'state.json')['status'],'READY')
    def test_real_r2_waits_preserves_and_can_resume_input_inspection(self):
        source=ROOT/'tests/fixtures/complex_R2.pdb';before=digest(source)
        self.assertEqual(self.call('--input',str(source),'--run-id','r2','--inspect-only'),3)
        job=self.root/'runs/r2';inv=read(job/'inventory.json')
        self.assertEqual(inv['atom_count'],3094)
        norm=read(job/'normalization.json');self.assertEqual(norm['retained_atom_count'],3094)
        self.assertEqual(len(norm['variants']),6)
        pending=read(job/'pending.json')['items'];self.assertEqual(len(pending),3)
        chem=read(job/'chemistry_result.json')
        self.assertEqual(chem['formal_charge'],2)
        decisions=self.decisions(construct_confirmed=True,termini={'r0188':'standard_charged'},
            chemistry={'c002':{'accept_candidate_sha256':chem['candidate_sha256']}})
        # Simulated acceptance exercises mechanics; this is not scientific approval.
        self.assertEqual(self.call('--resume',str(job),'--inspect-only','--decisions',decisions),0)
        self.assertEqual(self.call('--resume',str(job),'--inspect-only'),0)
        self.assertEqual(len(list((job/'normalized').iterdir())),3)
        self.assertEqual(digest(source),before)
        self.assertFalse((job/'md').exists())
    def test_missing_native_runtime_is_blocked_not_complete(self):
        source=self.fixture()
        with patch('modules.automation.backend.dependency_report',return_value={'ready':False,'status':'BLOCKED'}):
            result=self.call('--input',str(source),'--run-id','blocked','--smoke-test','--decisions',self.decisions(construct_confirmed=True))
        self.assertEqual(result,4)
        self.assertEqual(read(self.root/'runs/blocked/state.json')['status'],'DEPENDENCY_BLOCKED')
    def test_lock_cannot_change_another_running_job_state(self):
        source=self.fixture();self.call('--input',str(source),'--run-id','locked','--inspect-only')
        job=self.root/'runs/locked';before=(job/'state.json').read_bytes()
        with job_lock(job):
            self.assertEqual(self.call('--resume',str(job),'--inspect-only'),5)
            self.assertEqual((job/'state.json').read_bytes(),before)
    def test_raw_tampering_is_rejected(self):
        source=self.fixture();self.call('--input',str(source),'--run-id','tamper','--inspect-only')
        job=self.root/'runs/tamper';staged=job/read(job/'source.json')['staged'];staged.write_text('END\n')
        self.assertEqual(self.call('--resume',str(job),'--inspect-only'),2)
        self.assertEqual(read(job/'state.json')['status'],'FAILED')
    def test_path_mapping_portable_and_rejects_escape(self):
        old=self.root/'runs/old';new=self.root/'runs/new'
        data={'file':str(old/'normalized/protein.pdb'),'original':'/original/file.pdb'}
        stored=w.relocate(data,old);self.assertEqual(stored['file'],'@job/normalized/protein.pdb')
        self.assertEqual(w.relocate(stored,new,True)['file'],str(new/'normalized/protein.pdb'))
        with self.assertRaises(ValueError):w.relocate('@job/../../escape',new,True)
    def test_analysis_roles_remap_original_blank_chain(self):
        residues=[{'kind':'protein','original':{'chain':''},'chain':'B'}]
        self.assertEqual(w.map_analysis_groups({'receptor':{'chains':['']}},residues)['receptor']['chains'],['B'])
        residues.append({'kind':'protein','original':{'chain':''},'chain':'C'})
        with self.assertRaises(ValueError):w.map_analysis_groups({'receptor':{'chains':['']}},residues)
    def test_incompatible_forcefield_fails_before_native(self):
        source=self.fixture();cfg=self.root/'incompatible.json'
        write(cfg,{'MD':{'water_model':'opc'}})
        code=self.call('--input',str(source),'--run-id','badff','--inspect-only','--config',str(cfg),'--decisions',self.decisions(construct_confirmed=True))
        self.assertEqual(code,2);self.assertFalse((self.root/'runs/badff/md').exists())

    def test_completed_state_resume_and_upstream_invalidation_without_native_md(self):
        source=self.fixture();(self.root/'scripts').mkdir()
        shutil.copyfile(ROOT/'scripts/run_workflow.py',self.root/'scripts/run_workflow.py')
        job=self.root/'runs/mocked';calls=[]
        def stage(name,folder):
            def action(*args,**kwargs):
                calls.append(name);(job/folder).mkdir(parents=True,exist_ok=True)
                (job/folder/(name+'.txt')).write_text(name)
                return {'passed':True,'scope':'MOCK_TEST_ONLY'}
            return action
        def mapping(*args):
            write(job/'mappings/topology_atom_map.json',{'scope':'MOCK_TEST_ONLY'})
            write(job/'mappings/analysis_original_identifiers.json',{'exports':['analyzed'] if 'analysis' in calls else []})
        with contextlib.ExitStack() as stack:
            stack.enter_context(patch('modules.automation.backend.dependency_report',return_value={'ready':True,'compatibility_fingerprint':'mock'}))
            stack.enter_context(patch.object(w.platform,'system',return_value='Linux'))
            stack.enter_context(patch.object(w.platform,'machine',return_value='x86_64'))
            stack.enter_context(patch('modules.automation.analysis_mapping.add_original_identifiers',side_effect=mapping))
            for name,folder in [('prepare','md/prepared'),('build','md/system'),('run','md'),('analyze','analysis'),('audit_run','reports')]:
                label={'run':'md','analyze':'analysis','audit_run':'audit'}.get(name,name)
                stack.enter_context(patch('modules.automation.backend.'+name,side_effect=stage(label,folder)))
            self.assertEqual(self.call('--input',str(source),'--run-id','mocked','--decisions',self.decisions(construct_confirmed=True)),0)
            self.assertEqual(calls,['prepare','build','md','analysis','audit'])
            self.assertEqual(self.call('--resume',str(job)),0)
            self.assertEqual(calls,['prepare','build','md','analysis','audit'])
            self.assertEqual(read(job/'state.json')['status'],'COMPLETE')
            self.assertEqual(self.call('--resume',str(job),'--decisions',self.decisions(construct_confirmed=False)),2)
            self.assertEqual(calls,['prepare','build','md','analysis','audit'])

if __name__=='__main__':unittest.main()
