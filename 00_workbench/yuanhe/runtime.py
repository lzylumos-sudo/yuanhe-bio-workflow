"""Resolve execution environments consistently without mutating saved run settings."""
import copy
import os
from pathlib import Path
import re
import shlex
import shutil
import sys

SCHEMA_VERSION = 'yuanhe.runtime/1'


def validate_runtime_config(config):
    """Legacy files have no schema; unknown schemas must never be silently applied."""
    if not isinstance(config, dict):
        raise ValueError('Runtime configuration must be a JSON object')
    if 'schema_version' in config and config['schema_version'] != SCHEMA_VERSION:
        raise ValueError('Unsupported runtime schema: ' + str(config['schema_version']))
    return copy.deepcopy(config)


def migrated_runtime_config(config):
    """Return a new configuration; never rewrite a site file or frozen run."""
    migrated = validate_runtime_config(config)
    migrated['schema_version'] = SCHEMA_VERSION
    if not migrated.get('reinvent_python') and migrated.get('reinvent_executable'):
        migrated['reinvent_python'] = resolve_python(migrated, 'reinvent_python')
    return migrated


def resolve_python(config, key):
    """Accept explicit paths and the legacy REINVENT executable-only configuration."""
    validate_runtime_config(config)
    if key == 'current':
        return sys.executable
    if config.get(key):
        return str(config[key])
    if key in ('chem_python', 'screen_python'):
        return str(config.get('chem_python') or sys.executable)
    if key == 'reinvent_python' and config.get('reinvent_executable'):
        command = str(config['reinvent_executable'])
        executable = shutil.which(command)
        if not executable and not Path(command).is_absolute() and len(Path(command).parts) == 1:
            return ''
        executable = executable or command
        name = 'python.exe' if Path(executable).suffix.lower() == '.exe' else 'python'
        return str(Path(executable).absolute().parent / name)
    return ''


def runtime_env(executable, *, for_overlay=False):
    """Use the selected environment's tools, never another Python's user packages."""
    env = dict(os.environ)
    if executable:
        exe = shutil.which(str(executable)) or str(executable)
        bindir = Path(exe).absolute().parent
        env['PATH'] = str(bindir) + os.pathsep + env.get('PATH', '')
        # An inherited AMBERHOME can silently select another installation.
        env.pop('AMBERHOME', None)
        if (bindir.parent / 'dat/leap').is_dir():
            env['AMBERHOME'] = str(bindir.parent)
    env.pop('PYTHONHOME', None)
    env.pop('PYTHONPATH', None)
    if for_overlay:
        # Some standalone adapters overlay this mapping on os.environ. Empty values
        # prevent that merge from restoring a different Python's import settings.
        env.update(PYTHONHOME='', PYTHONPATH='')
        env.setdefault('AMBERHOME', '')
    env['PYTHONNOUSERSITE'] = '1'
    return env


def resolve_executable(config, name, python=''):
    """Honor configured tool paths; search names in the execution environment PATH."""
    return shutil.which(str(config.get(name) or name), path=runtime_env(python).get('PATH', ''))


def _launcher_python(executable, python):
    """Inspect pip/uv console scripts without executing arbitrary shell wrappers."""
    try:
        with Path(executable).open('rb') as handle:
            lines = handle.read(8192).decode('utf-8', errors='replace').splitlines()
        if not lines or not lines[0].startswith('#!'):
            return '', 'unrecognized'
        command = shlex.split(lines[0][2:].strip())
        kind = 'python_shebang'
        if command and Path(command[0]).name == 'env':
            command = command[1:]
            if command and command[0] == '-S':
                command = shlex.split(' '.join(command[1:]))
            elif len(command) != 1:
                # Without -S the kernel passes the entire tail as one env argument.
                return '', 'unrecognized_interpreter'
            if not command:
                return '', 'unrecognized_interpreter'
            resolved = shutil.which(command[0], path=runtime_env(python)['PATH'])
            command = [resolved] if resolved else []
            kind = 'env_python_shebang'
        if command and Path(command[0]).name in ('sh', 'bash'):
            # distlib uses this standard polyglot wrapper for long or spaced paths.
            # An arbitrary shell script can redirect Python elsewhere, so do not guess.
            if len(lines) < 3 or lines[2].strip() != "' '''":
                return '', 'unrecognized_shell_wrapper'
            command = shlex.split(lines[1])
            if len(command) != 4 or command[0] != 'exec' or command[2:] != ['$0', '$@']:
                return '', 'unrecognized_shell_wrapper'
            command = command[1:2]
            kind = 'distlib_shell_wrapper'
        if (not command or not Path(command[0]).is_absolute() or
                not re.fullmatch(r'python(?:\d+(?:\.\d+)*)?(?:\.exe)?', Path(command[0]).name)):
            return '', 'unrecognized_interpreter'
        return command[0], kind
    except (OSError, UnicodeError, ValueError):
        return '', 'unreadable'


def reinvent_environment_check(config):
    """Check the launcher interpreter too: copied launchers can retain an old shebang."""
    python = resolve_python(config, 'reinvent_python')
    python = (shutil.which(python) or python) if python else ''
    executable = resolve_executable(config, 'reinvent_executable', python) or config.get('reinvent_executable', '')
    launcher_python, launcher_kind = _launcher_python(executable, python) if executable else ('', 'missing')
    ok = bool(python and executable and launcher_python and
              Path(python).absolute().parent == Path(executable).absolute().parent ==
              Path(launcher_python).absolute().parent)
    return {'ok': ok, 'python': python, 'executable': executable,
            'launcher_python': launcher_python, 'launcher_kind': launcher_kind,
            'detail': 'Use the REINVENT launcher and Python from the same environment bin directory; '
                      'the launcher shebang must also select that environment. Standard pip/uv long-path '
                      'wrappers and env Python shebangs use the selected runtime. Reinstall the launcher '
                      'if it still points to an old environment. Arbitrary shell wrappers are not automatically certified. '
                      'No configuration was rewritten.'}


def vina_version(output):
    """Official Vina prints a v-prefixed version; retain compatibility without v."""
    match = re.search(r'\bVina\s+v?(\d+\.\d+\.\d+)\b', output, re.I)
    return match[1] if match else None
