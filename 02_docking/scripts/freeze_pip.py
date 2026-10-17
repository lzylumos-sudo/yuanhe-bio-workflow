import sys,importlib.metadata
from pathlib import Path
rows=[]
for d in importlib.metadata.distributions():
 name=d.metadata['Name']
 if name.lower().startswith('nvidia-') or name.lower() in {'playwright','pyee','greenlet','molscrub'}:rows.append(f'{name}=={d.version}')
(Path(sys.argv[1])/'env/requirements-pip.lock').write_text('# pip supplements only; conda packages are pinned separately\n'+'\n'.join(sorted(rows,key=str.lower))+'\n')
