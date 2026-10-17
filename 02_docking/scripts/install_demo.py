import sys,subprocess,json
from pathlib import Path
from sbvs.paths import ROOT,write_json
from sbvs.preparation.receptor import prepare
from sbvs.preparation.box import make_box
sys.path.insert(0,str(ROOT/'scripts'))
from download_data import main as download
download();prepare(ROOT/'data/raw/1DRF.pdb',ROOT/'data/prepared/demo_1drf',classroom=True);make_box(ROOT/'data/reference/FOL_crystal.sdf',5,ROOT/'data/prepared/demo_1drf/box.json')
subprocess.run([sys.executable,str(ROOT/'scripts/cco_smoke.py')],check=True)
write_json(ROOT/'reports/class_ready.json',{'ready':True,'evidence':'reports/cco_cli_smoke.json','browser_acceptance':'separate test required'})
