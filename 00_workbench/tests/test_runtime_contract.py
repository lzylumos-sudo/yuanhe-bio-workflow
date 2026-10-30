"""Synthetic deployment contracts; stdlib only, no molecule or model computation."""
import copy
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / '00_workbench'))
from yuanhe import registry, runtime
from webui import readiness, workflow, studio


class RuntimeContractTests(unittest.TestCase):
    def test_legacy_reinvent_is_resolved_without_mutating_config(self):
        with tempfile.TemporaryDirectory() as folder:
            cfg = {'reinvent_executable': str(Path(folder) / 'bin/reinvent'), 'research': {'seed': 17}}
            before = copy.deepcopy(cfg)
            self.assertEqual(runtime.resolve_python(cfg, 'reinvent_python'), str(Path(folder) / 'bin/python'))
            self.assertEqual(cfg, before)
            cfg['reinvent_python'] = '/explicit/model/python'
            self.assertEqual(runtime.resolve_python(cfg, 'reinvent_python'), '/explicit/model/python')

    def test_missing_bare_reinvent_command_is_not_guessed_from_cwd(self):
        with patch('yuanhe.runtime.shutil.which', return_value=None):
            self.assertEqual(runtime.resolve_python({'reinvent_executable': 'missing-reinvent'}, 'reinvent_python'), '')

    def test_reinvent_explicit_environment_mismatch_is_required_diagnostic_failure(self):
        with tempfile.TemporaryDirectory() as folder:
            cfg = {'reinvent_executable': str(Path(folder) / 'model/bin/reinvent')}
            launcher = Path(cfg['reinvent_executable']); launcher.parent.mkdir(parents=True)
            launcher.write_text('#!' + (launcher.parent / 'python').as_posix() + '\n')
            self.assertTrue(runtime.reinvent_environment_check(cfg)['ok'])
            cfg['reinvent_python'] = str(Path(folder) / 'other/bin/python')
            before = dict(cfg)
            self.assertFalse(runtime.reinvent_environment_check(cfg)['ok'])
            with patch.object(readiness, 'probe', return_value={'ok': False, 'error': 'synthetic probe'}):
                result = readiness.collect(cfg, Path(folder))
            self.assertFalse(result['checks']['reinvent_environment']['ok'])
            self.assertTrue(result['checks']['reinvent_environment']['required'])
            self.assertEqual(cfg, before)

    @unittest.skipIf(os.name == 'nt', 'Linux console-script interpreter contract')
    def test_reinvent_copied_launcher_cannot_hide_old_environment(self):
        with tempfile.TemporaryDirectory() as folder:
            bindir = Path(folder) / 'new/bin'; bindir.mkdir(parents=True)
            launcher = bindir / 'reinvent'
            current, old = str(bindir / 'python'), str(Path(folder) / 'old/bin/python')
            Path(current).write_text('synthetic interpreter'); Path(current).chmod(0o755)
            cfg = {'reinvent_python': current, 'reinvent_executable': str(launcher)}
            before = dict(cfg)
            for text, expected, kind in [
                ('#!' + current + '\n', True, 'python_shebang'),
                ('#!' + old + '\n', False, 'python_shebang'),
                ('#!/bin/sh\n' + "'''exec' " + current + ' "$0" "$@"\n' + "' '''\n", True, 'distlib_shell_wrapper'),
                ('#!/bin/sh\n' + "'''exec' " + old + ' "$0" "$@"\n' + "' '''\n", False, 'distlib_shell_wrapper'),
                ('#!/usr/bin/env python\n', True, 'env_python_shebang'),
                ('#!/usr/bin/env -S python -s\n', True, 'env_python_shebang'),
                ('#!/usr/bin/env python -s\n', False, 'unrecognized_interpreter'),
                ('#!/bin/sh\nexec ' + old + ' -m reinvent "$@"\n', False, 'unrecognized_shell_wrapper'),
            ]:
                with self.subTest(text=text):
                    launcher.write_text(text); launcher.chmod(0o755)
                    checked = runtime.reinvent_environment_check(cfg)
                    self.assertEqual(checked['ok'], expected, checked)
                    self.assertEqual(checked['launcher_kind'], kind)
                    self.assertEqual(cfg, before)

    @unittest.skipIf(os.name == 'nt', 'Linux console-script interpreter contract')
    def test_reinvent_standard_wrapper_supports_spaces_in_prefix(self):
        with tempfile.TemporaryDirectory(prefix='runtime contract ') as folder:
            bindir = Path(folder) / 'bin'; bindir.mkdir()
            launcher = bindir / 'reinvent'; python = str(bindir / 'python')
            launcher.write_text('#!/bin/sh\n' + "'''exec' " + repr(python) + ' "$0" "$@"\n' + "' '''\n")
            launcher.chmod(0o755)
            self.assertTrue(runtime.reinvent_environment_check({
                'reinvent_python': python, 'reinvent_executable': str(launcher)})['ok'])

    def test_legacy_and_current_schema_work_future_schema_blocks(self):
        for cfg in ({}, {'schema_version': runtime.SCHEMA_VERSION}):
            self.assertEqual(runtime.validate_runtime_config(cfg), cfg)
        for cfg in ({'schema_version': 'yuanhe.runtime/99'}, {'schema_version': None}, []):
            with self.assertRaises(ValueError):
                runtime.validate_runtime_config(cfg)

    def test_probe_and_workers_share_environment_helper(self):
        self.assertIs(readiness.runtime_env, runtime.runtime_env)
        with tempfile.TemporaryDirectory() as folder:
            prefix = Path(folder)
            (prefix / 'bin').mkdir()
            (prefix / 'dat/leap').mkdir(parents=True)
            with patch.dict(os.environ, {'PYTHONHOME': '/wrong', 'PYTHONPATH': '/wrong', 'AMBERHOME': '/wrong'}):
                env = runtime.runtime_env(prefix / 'bin/python')
                self.assertEqual(env['PATH'].split(os.pathsep)[0], str(prefix / 'bin'))
                self.assertEqual(env['AMBERHOME'], str(prefix))
                self.assertNotIn('PYTHONHOME', env)
                self.assertNotIn('PYTHONPATH', env)
                self.assertEqual(env['PYTHONNOUSERSITE'], '1')
                merged = {**os.environ, **runtime.runtime_env(prefix / 'bin/python', for_overlay=True)}
                self.assertEqual(merged['PYTHONHOME'], '')
                self.assertEqual(merged['PYTHONPATH'], '')
                other = runtime.runtime_env(prefix/'other/bin/python', for_overlay=True)
                self.assertEqual({**os.environ, **other}['AMBERHOME'], '')

    def test_real_import_probe_ignores_contaminated_python_environment(self):
        with patch.dict(os.environ, {'PYTHONHOME': '/invalid/python', 'PYTHONPATH': '/invalid/packages'}):
            result = registry.probe_python(sys.executable, ['json'])
        self.assertTrue(result['available'], result)

    def test_tool_found_in_selected_environment_when_service_path_has_none(self):
        with tempfile.TemporaryDirectory() as folder:
            bindir = Path(folder) / 'bin'; bindir.mkdir()
            name = 'antechamber.exe' if os.name == 'nt' else 'antechamber'
            tool = bindir / name; tool.write_text('synthetic executable fixture'); tool.chmod(0o755)
            with patch.dict(os.environ, {'PATH': ''}):
                self.assertEqual(Path(runtime.resolve_executable({}, 'antechamber', bindir / 'python')), tool)
                self.assertIsNone(runtime.resolve_executable({'antechamber': str(bindir / 'missing')}, 'antechamber', bindir / 'python'))

    def test_vina_versions(self):
        for text in ('AutoDock Vina v1.2.7', 'AutoDock Vina 1.2.7\n', 'Vina v1.2.7\n'):
            self.assertEqual(runtime.vina_version(text), '1.2.7')
        self.assertEqual(runtime.vina_version('AutoDock Vina v1.2.6'), '1.2.6')
        self.assertIsNone(runtime.vina_version('unknown program 1.2.7'))

    def test_vina_registry_requires_matching_version_and_success(self):
        real = registry.load()
        catalog = {'modules': [], 'backends': [next(b for b in real['backends'] if b['backend_id'] == 'vina')]}
        for output, code, expected in [('AutoDock Vina v1.2.7', 0, True), ('AutoDock Vina 1.2.7', 0, True),
                                       ('AutoDock Vina v1.2.6', 0, False), ('AutoDock Vina v1.2.7', 1, False)]:
            with self.subTest(output=output, code=code), patch.object(registry, 'load', return_value=catalog), \
                    patch.object(registry, 'probe_python', return_value={'available': True}), \
                    patch.object(registry, 'resolve_executable', return_value='/synthetic/vina'), \
                    patch.object(registry.subprocess, 'run', return_value=subprocess.CompletedProcess([], code, output, '')):
                result = registry.snapshot({'chem_python': sys.executable})['backends'][0]
                self.assertEqual(result['runnable'], expected)
                self.assertFalse(result['verified'])

    def test_registry_uses_md_python_and_optional_pymol_does_not_block_analysis(self):
        catalog = registry.load()
        by_id = {b['backend_id']: b for b in catalog['backends']}
        self.assertEqual(by_id['openmm_system_builder']['environment'], 'md_python')
        self.assertNotIn('pymol', by_id['zlis']['executables'])
        self.assertIn('pymol', by_id['zlis']['optional_executables'])
        only = {'modules': [], 'backends': [by_id['openmm_system_builder'], by_id['zlis']]}
        def tool(cfg, name, python=''):
            return str(Path(python).parent / name) if name == 'antechamber' else None
        with patch.object(registry, 'load', return_value=only), \
                patch.object(registry, 'probe_python', return_value={'available': True}), \
                patch.object(registry, 'resolve_executable', side_effect=tool):
            result = registry.snapshot({'md_python': '/md/bin/python', 'prep_python': '/prep/bin/python',
                                        'chem_python': '/chem/bin/python'})['backends']
        self.assertEqual(result[0]['resolved_environment']['python'], '/md/bin/python')
        self.assertEqual(result[0]['executable_paths']['antechamber'], str(Path('/md/bin/antechamber')))
        self.assertTrue(result[1]['installed'])
        self.assertIsNone(result[1]['optional_executable_paths']['pymol'])

    def test_module_reports_capability_availability_without_requiring_optional_backend(self):
        catalog = registry.load()
        catalog['modules'] = [m for m in catalog['modules'] if m['module_id'] == 'MD']
        catalog['backends'] = [b for b in catalog['backends'] if b['module_id'] == 'MD']
        # A synthetic missing optional soluble backend must not hide installed membrane capabilities.
        next(b for b in catalog['backends'] if b['backend_id'] == 'openmm_soluble')['environment'] = 'missing_python'
        with patch.object(registry, 'load', return_value=catalog), \
                patch.object(registry, 'probe_python', side_effect=lambda exe, packages: {'available': bool(exe)}), \
                patch.object(registry, 'resolve_executable', side_effect=lambda cfg, name, exe='': '/synthetic/'+name if exe else None):
            result = registry.snapshot({'md_python': '/md/python'})
        md = result['modules'][0]
        self.assertTrue(md['installed'])
        self.assertTrue(md['capability_status']['membrane_system_building']['installed'])
        self.assertFalse(md['verified'])
        self.assertEqual(md['readiness'], 'Needs Validation')
        # With only GROMACS installed, production is visible but construction stays missing.
        next(b for b in catalog['backends'] if b['backend_id'] == 'openmm_system_builder')['environment'] = 'missing_python'
        with patch.object(registry, 'load', return_value=catalog), \
                patch.object(registry, 'probe_python', side_effect=lambda exe, packages: {'available': bool(exe)}), \
                patch.object(registry, 'resolve_executable', side_effect=lambda cfg, name, exe='': '/synthetic/'+name if exe else None):
            partial = registry.snapshot({'md_python': '/md/python'})['modules'][0]
        self.assertEqual(partial['readiness'], 'Partially Installed')
        self.assertFalse(partial['capability_status']['membrane_system_building']['installed'])
        self.assertTrue(partial['capability_status']['production']['installed'])

    def test_defaults_accept_explicit_interpreter_and_environment_save_preserves_history(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); (root / 'deployment').mkdir()
            site = root / 'deployment/runtime.json'
            site.write_text(json.dumps({'schema_version': runtime.SCHEMA_VERSION, 'reinvent_python': '/explicit/python'}))
            historical = root / 'frozen-config.json'; historical.write_text('{"original":true}')
            with patch.object(workflow, 'ROOT', root):
                self.assertEqual(workflow.defaults()['reinvent_python'], '/explicit/python')
                studio.environment_save({'reinvent_python': '/new/python'})
                self.assertEqual(workflow.defaults()['reinvent_python'], '/new/python')
                self.assertEqual(historical.read_text(), '{"original":true}')
                site.write_text('{"schema_version":"yuanhe.runtime/99"}')
                with self.assertRaises(ValueError):
                    workflow.defaults()
                with self.assertRaises(ValueError):
                    studio.environment_save({'reinvent_python': '/new/python'})

    def test_migration_creates_new_file_preserves_source_and_rejects_overwrite(self):
        spec = importlib.util.spec_from_file_location('migrate_runtime_test', ROOT / 'deployment/migrate_runtime.py')
        module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as folder:
            source, target = Path(folder) / 'old.json', Path(folder) / 'new.json'
            source.write_text(json.dumps({'reinvent_executable': str(Path(folder) / 'model/bin/reinvent'),
                                         'research': {'custom_policy': 'retain'}}))
            original = source.read_bytes()
            result = module.migrate(source, target)
            self.assertEqual(source.read_bytes(), original)
            self.assertEqual(result['schema_version'], runtime.SCHEMA_VERSION)
            self.assertEqual(result['research'], {'custom_policy': 'retain'})
            self.assertEqual(result['migration']['source_schema'], 'legacy-unversioned')
            self.assertEqual(result['migration']['execution_validation'], 'not_run')
            with self.assertRaises(FileExistsError):
                module.migrate(source, target)
            with self.assertRaises(ValueError):
                module.migrate(source, source)
            source.write_text('{"schema_version":"yuanhe.runtime/99"}')
            future = Path(folder) / 'future.json'
            with self.assertRaises(ValueError):
                module.migrate(source, future)
            self.assertFalse(future.exists())


if __name__ == '__main__':
    unittest.main()
