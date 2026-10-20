from pathlib import Path
import json
r=Path(__file__).resolve().parents[1]
lines=(r/'public/_headers').read_text().splitlines()
headers=[]
for l in lines[1:6]:
 k,v=l.strip().split(': ',1);headers.append({'key':k,'value':v})
(r/'vercel.json').write_text(json.dumps({'buildCommand':'pnpm run build','outputDirectory':'dist','headers':[{'source':'/(.*)','headers':headers}]},indent=2)+'\n')
