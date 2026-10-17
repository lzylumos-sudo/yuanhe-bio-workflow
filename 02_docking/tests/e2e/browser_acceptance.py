import sys,os,time,json,re,csv,math
from pathlib import Path
from playwright.sync_api import sync_playwright,expect
from sbvs.paths import ROOT,write_json,read_json
from sbvs.export import read_rows
evidence=ROOT/'reports/browser';evidence.mkdir(parents=True,exist_ok=True);audits=[]
with sync_playwright() as p:
 browser=p.chromium.launch(headless=True,args=['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader'])
 context=browser.new_context(viewport={'width':1440,'height':1100})
 # Browser internet is blocked; local assets and computations must remain usable.
 context.route('**/*',lambda route:route.continue_() if route.request.url.startswith('http://127.0.0.1:5902/') else route.abort())
 page=context.new_page();errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
 for i in range(3):
  page.goto('http://127.0.0.1:5902/');expect(page.locator('#results')).to_be_hidden();expect(page.locator('#ranking tr')).to_have_count(0)
  if i==0:page.screenshot(path=str(evidence/'initial.png'),full_page=True)
  page.get_by_label('Ligand SMILES',exact=True).fill('CCO');start=time.monotonic();page.get_by_role('button',name='Start screening',exact=True).click()
  expect(page.locator('#ranking tr')).to_have_count(1,timeout=180000);elapsed=time.monotonic()-start
  rid=page.locator('#jobid').inner_text().split(': ',1)[1];rundir=ROOT/'runs'/rid;rows=read_rows(rundir/'ranking.csv');score=float(rows[0]['vina_affinity_kcal_mol']);assert math.isfinite(score)
  assert page.locator('#ranking tr td').nth(2).inner_text()==rows[0]['vina_affinity_kcal_mol']
  log=next((rundir/'entries').glob('*/stdout.log'));assert 'AutoDock Vina' in log.read_text()
  m=read_json(rundir/'run_manifest.json');assert m['config']['engine']=='vina' and m['config']['device']=='cpu'
  page.wait_for_function('window.sbvsViewer !== undefined');page.screenshot(path=str(evidence/f'cco_{i+1}.png'),full_page=True)
  audits.append({'run_id':rid,'elapsed_seconds':elapsed,'affinity_kcal_mol':score,'screenshot':str((evidence/f'cco_{i+1}.png').relative_to(ROOT)),'status':'passed'})
 assert len({a['run_id'] for a in audits})==3
 write_json(ROOT/'reports/browser_acceptance.json',{'tests':audits,'browser_external_network':'blocked','page_errors':errors})
 if errors:raise AssertionError(errors)
 browser.close()
print(json.dumps(audits,indent=2))
