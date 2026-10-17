import json
from playwright.sync_api import sync_playwright,expect
from sbvs.paths import ROOT,read_json,write_json
from sbvs.export import report
# Import the completed CLI class run into the same queue catalogue for review only.
out=ROOT/read_json(ROOT/'reports/class20_result.json')['run'];m=read_json(out/'run_manifest.json');rid=m['run_id']
with sync_playwright() as p:
 b=p.chromium.launch(headless=True,args=['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']);page=b.new_page(viewport={'width':1440,'height':1100});page.goto('http://127.0.0.1:5902/?run='+rid);expect(page.locator('#ranking tr')).to_have_count(20);page.wait_for_function('window.sbvsViewer !== undefined')
 for i in (0,1):
  page.locator('#ranking tr').nth(i).click();page.wait_for_timeout(600);page.get_by_role('button',name='保存当前姿态 PNG',exact=True).click();expect(page.locator('#poseLabel')).to_contain_text('PNG 已保存',timeout=10000)
 page.screenshot(path=str(ROOT/'reports/browser/class20.png'),full_page=True);b.close()
images=list((out/'images').glob('*.png'));assert len(images)>=2 and images[0].read_bytes()!=images[1].read_bytes();report(out)
write_json(ROOT/'reports/class_report_images.json',{'run':str(out.relative_to(ROOT)),'images':[str(x.relative_to(ROOT)) for x in images],'status':'passed'});print(str(out/'report.html'))
