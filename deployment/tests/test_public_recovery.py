"""Fault injection for public-run orchestration; no scientific execution claimed."""
import json
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch
from types import SimpleNamespace

PACKAGE = Path(__file__).resolve().parents[1] / 'public_minimal'
sys.path.insert(0, str(PACKAGE))
from recovery import Recovery, lease, execute, inventory, write


class RecoveryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.out = Path(self.temp.name)
        self.identity = {'version': 'test'}

    def save(self, name, text='original'):
        path = self.out / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)

    def failed(self):
        raise RuntimeError('injected failure')

    def prepared(self):
        recovery = Recovery(self.out, self.identity)
        recovery.phase('prepare', lambda: self.save('inputs/input.pdb'))
        return recovery

    def test_completed_stages_skip_after_verification(self):
        self.prepared()
        resumed = Recovery(self.out, self.identity, resume=True)
        resumed.phase('prepare', self.failed)
        self.assertEqual(resumed.state['completed'], ['prepare'])

    def test_failure_is_archived_and_partial_outputs_removed(self):
        recovery = self.prepared()
        def failing():
            self.save('screening/partial')
            self.save('screen.log', 'failed log')
            self.failed()
        with self.assertRaises(RuntimeError):
            recovery.phase('screen', failing)
        resumed = Recovery(self.out, self.identity, resume=True)
        self.assertFalse((self.out / 'screening').exists())
        self.assertTrue(list((self.out / '.recovery').glob('interrupted-*/screen.log')))
        self.assertTrue(list((self.out / '.recovery').glob('interrupted-*/FAILED.json')))
        resumed.phase('screen', lambda: self.save('screening/result'))
        resumed.finish()
        self.assertFalse((self.out / 'FAILED.json').exists())
        self.assertEqual(resumed.state['status'], 'completed')

    def test_identity_change_refused_without_rewriting_state(self):
        self.prepared()
        original = (self.out / 'RUN_STATE.json').read_bytes()
        with self.assertRaises(ValueError):
            Recovery(self.out, {'version': 'changed'}, resume=True)
        self.assertEqual(original, (self.out / 'RUN_STATE.json').read_bytes())

    def test_completed_input_tampering_refused(self):
        self.prepared()
        self.save('inputs/input.pdb', 'tampered')
        with self.assertRaises(ValueError):
            Recovery(self.out, self.identity, resume=True)

    def test_interrupted_upstream_tampering_refused(self):
        recovery = self.prepared()
        with self.assertRaises(RuntimeError):
            recovery.phase('dock', self.failed)
        self.save('inputs/input.pdb', 'tampered')
        with self.assertRaises(ValueError):
            Recovery(self.out, self.identity, resume=True)
        self.assertEqual((self.out / 'inputs/input.pdb').read_text(), 'tampered')

    def test_dynamics_checkpoint_retained_for_native_backend(self):
        recovery = self.prepared()
        checkpoint = 'dynamics/md/replica_01/checkpoint.chk'
        def failing():
            self.save(checkpoint, 'checkpoint bytes')
            self.failed()
        with self.assertRaises(RuntimeError):
            recovery.phase('dynamics', failing)
        resumed = Recovery(self.out, self.identity, resume=True)
        self.assertEqual((self.out / checkpoint).read_text(), 'checkpoint bytes')
        self.assertTrue(list((self.out / '.recovery').glob('interrupted-*/' + checkpoint)))
        resumed.phase('dynamics', lambda: self.save('md_execution.json', '{}'))

    def test_dynamics_does_not_allow_upstream_changes(self):
        recovery = self.prepared()
        with self.assertRaises(RuntimeError):
            recovery.phase('dynamics', self.failed)
        self.save('inputs/input.pdb', 'tampered')
        with self.assertRaises(ValueError):
            Recovery(self.out, self.identity, resume=True)

    def test_interrupted_analysis_preserves_committed_md(self):
        recovery = self.prepared()
        recovery.phase('dynamics', lambda: self.save('dynamics/md/replica_01/final.chk'))
        def failing():
            self.save('dynamics/analysis/partial.csv')
            self.failed()
        with self.assertRaises(RuntimeError):
            recovery.phase('analyze', failing)
        resumed = Recovery(self.out, self.identity, resume=True)
        resumed.phase('dynamics', self.failed)
        self.assertTrue((self.out / 'dynamics/md/replica_01/final.chk').is_file())
        self.assertFalse((self.out / 'dynamics/analysis').exists())

    def test_lock_excludes_second_process(self):
        code = "import sys; from pathlib import Path; sys.path.insert(0,sys.argv[1]); from recovery import lease\nwith lease(Path(sys.argv[2])): pass"
        with lease(self.out):
            result = subprocess.run([sys.executable, '-c', code, str(PACKAGE), str(self.out)], capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(b'Another process', result.stderr)

    def test_timeout_stops_phase_before_releasing_lock(self):
        with lease(self.out) as fd:
            with self.assertRaises(subprocess.TimeoutExpired):
                execute([sys.executable, '-c', 'import time; time.sleep(30)'], self.out / 'timeout.log', fd, timeout=.05)
        with lease(self.out):
            pass

    def test_sigterm_records_interruption_and_can_resume(self):
        code = """
import sys
from pathlib import Path
sys.path.insert(0,sys.argv[1])
from recovery import Recovery, lease, interruption_handlers, execute
out=Path(sys.argv[2])
with lease(out) as fd, interruption_handlers():
    r=Recovery(out,{'version':'test'})
    r.phase('prepare',lambda:execute([sys.executable,'-c','import time; print(\"ready\",flush=True); time.sleep(30)'],out/'prepare.log',fd))
"""
        child = subprocess.Popen([sys.executable, '-c', code, str(PACKAGE), str(self.out)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        try:
            deadline = time.monotonic() + 5
            while time.monotonic() < deadline:
                if (self.out / 'prepare.log').exists() and 'ready' in (self.out / 'prepare.log').read_text():
                    break
                time.sleep(.01)
            else:
                self.fail('phase did not start')
            child.send_signal(signal.SIGTERM)
            child.wait(timeout=15)
            state = json.loads((self.out / 'RUN_STATE.json').read_text())
            self.assertEqual(state['status'], 'interrupted')
            with lease(self.out):
                resumed = Recovery(self.out, self.identity, resume=True)
                resumed.phase('prepare', lambda: self.save('prepared.pdb'))
        finally:
            if child.poll() is None:
                child.kill()
                child.wait()

    def test_public_runner_resume_end_to_end_with_stubbed_science(self):
        import run_public
        class FakeSystem:
            def addParticle(self, value): pass
        class FakeContext:
            def __init__(self, *args): pass
            def setPositions(self, value): pass
            def getPlatform(self): return self
            def getName(self): return 'CUDA'
            def getPropertyValue(self, *args): return 'test GPU'
        fake_openmm = SimpleNamespace(System=FakeSystem, Context=FakeContext,
            VerletIntegrator=lambda value: None, unit=SimpleNamespace(picoseconds=1),
            Platform=SimpleNamespace(getPlatformByName=lambda name: None))
        calls = []
        fail = [True]
        def phase(command, log, fd):
            name = command[-1]
            calls.append(name)
            self.save(log.name, name)
            if name == 'dynamics':
                self.save('dynamics/md/replica_01/native.chk', 'retained')
                if fail[0]:
                    fail[0] = False
                    raise RuntimeError('injected MD failure')
            if name == 'analyze':
                write(self.out / 'md_audit.json', {'passed': True, 'formal_duration_met': False})
        with patch.dict(sys.modules, {'openmm': fake_openmm}), \
             patch.object(run_public.subprocess, 'run', return_value=SimpleNamespace(returncode=0, stdout='ok', stderr='')), \
             patch.object(run_public, 'fetch_assets', side_effect=lambda out: self.save('sources/input')), \
             patch.object(run_public, 'execute', side_effect=phase):
            with lease(self.out) as fd:
                with self.assertRaises(RuntimeError):
                    run_public.main(self.out, False, fd)
            with lease(self.out) as fd:
                run_public.main(self.out, True, fd)
            self.assertEqual(calls, ['prepare', 'dock', 'screen', 'build', 'dynamics', 'dynamics', 'analyze'])
            self.assertTrue((self.out / 'COMPLETED.json').exists())
            self.assertFalse((self.out / 'FAILED.json').exists())
            with lease(self.out) as fd:
                run_public.main(self.out, True, fd)
            self.assertEqual(len(calls), 7)

    def test_second_interruption_during_restore_is_recoverable(self):
        recovery = self.prepared()
        with self.assertRaises(RuntimeError):
            recovery.phase('screen', self.failed)
        recovery.state['active']['restoring'] = True
        recovery.save()
        (self.out / 'inputs/input.pdb').unlink()
        resumed = Recovery(self.out, self.identity, resume=True)
        self.assertEqual((self.out / 'inputs/input.pdb').read_text(), 'original')
        resumed.phase('screen', lambda: self.save('screening/result'))

    def test_unrecorded_and_legacy_outputs_refused(self):
        self.save('old.pdb')
        with self.assertRaises(ValueError):
            Recovery(self.out, self.identity)
        with self.assertRaises(FileNotFoundError):
            Recovery(self.out, self.identity, resume=True)

    def test_symlink_refused(self):
        self.prepared()
        (self.out / 'link').symlink_to(self.out / 'inputs/input.pdb')
        with self.assertRaises(ValueError):
            Recovery(self.out, self.identity, resume=True)


if __name__ == '__main__':
    unittest.main()
