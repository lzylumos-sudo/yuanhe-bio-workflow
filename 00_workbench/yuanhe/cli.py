"""Repository-native CLI; no installation or computation on validation."""
import argparse
import json
from pathlib import Path
import sys
from .version import identity
from .registry import load, snapshot, validate
from .configuration import protocol_for, fingerprint, merge_layers
from .errors import envelope

def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == 'benchmark':
        from scientific.mdm2.cli import main as benchmark_main
        return benchmark_main(argv[1:])
    p = argparse.ArgumentParser(prog='yuanhe')
    sub = p.add_subparsers(dest='command', required=True)
    for command in ('version', 'doctor', 'modules', 'protocols'): sub.add_parser(command)
    sub.add_parser('benchmark', help='Frozen soluble-PPI benchmark and development throughput')
    for command in ('validate', 'run'):
        q = sub.add_parser(command); q.add_argument('--config', type=Path)
        q.add_argument('--target-config', type=Path)
        q.add_argument('--override', type=Path)
    q = sub.add_parser('status'); q.add_argument('run_id')
    args = p.parse_args(argv)
    try:
        if args.command == 'version': result = identity()
        elif args.command == 'protocols': result = load('protocol_registry')
        elif args.command in ('modules', 'doctor'):
            from webui.workflow import defaults
            result = snapshot(defaults())
            if args.command == 'doctor': result['registry_validation'] = validate()
        elif args.command == 'status':
            from webui.workflow import status
            result = status(args.run_id)
        else:
            from webui.workflow import validate_config, create
            layers = [(name, json.loads(path.read_text())) for name, path in
                      [('target', args.target_config), ('task', args.config), ('user_override', args.override)] if path]
            resolution = merge_layers(layers)
            cfg, _, _ = validate_config(resolution['config'])
            result = {'valid': True, 'config': cfg, 'config_hash': fingerprint(cfg),
                      'protocol': protocol_for(cfg), 'resolution': resolution}
            if args.command == 'run':
                result = {**create(resolution['config']), 'status': 'created',
                          'note': 'No compute started. Existing start/review APIs control execution.'}
        print(json.dumps(result, indent=2, allow_nan=False)); return 0
    except (ValueError, OSError, RuntimeError) as exc:
        print(json.dumps(envelope(str(exc))), file=sys.stderr); return 2

if __name__ == '__main__': raise SystemExit(main())
