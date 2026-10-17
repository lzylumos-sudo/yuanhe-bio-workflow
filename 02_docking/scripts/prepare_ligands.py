import sys
from sbvs.cli import main
sys.argv.insert(1,'prepare-ligands')
raise SystemExit(main())
