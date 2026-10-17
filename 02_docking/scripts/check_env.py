import sys
from sbvs.cli import main
sys.argv.insert(1,'doctor')
raise SystemExit(main())
