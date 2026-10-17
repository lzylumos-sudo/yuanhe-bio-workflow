from playwright.sync_api import sync_playwright,expect
from sbvs.paths import ROOT,read_json,write_json
import time
results=[]
with sync_playwright() as p:
 b=p.chromium.launch(headless=True,args=['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']);page=b.new_page(viewport={'width':1440,'height':1100});errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
 for route in ('classroom','history','environment','research'):
  page.goto('http://127.0.0.1:5902/'+route);page.wait_for_load_state('networkidle');assert page.locator('h1').count()==1
  if route=='classroom':
   page.locator('[data-evidence="box"] summary').click();expect(page.locator('[data-evidence="box"] pre')).to_contain_text('30.22840625')
  if route=='research':
   rid=read_json(ROOT/'reports/research_acceptance.json')['custom_pdb_review']['receptor_id'];page.locator('#receptorSelect').select_option(rid);expect(page.locator('#pocketPreview')).to_be_visible();page.wait_for_timeout(500)
  page.screenshot(path=str(ROOT/f'reports/browser/{route}.png'),full_page=True);results.append({'page':route,'status':'passed'})
 assert not errors,errors
 b.close()
write_json(ROOT/'reports/pages_acceptance.json',{'pages':results,'page_errors':errors});print('All four secondary pages rendered without JS errors')
