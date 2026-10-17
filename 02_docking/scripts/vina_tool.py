import sys
from sbvs.cli import main
sys.argv.insert(1,'screen')
raise SystemExit(main())
