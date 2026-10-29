const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path');
const root=path.resolve(__dirname,'../..');
(async()=>{
  const browser=await chromium.launch({headless:true,channel:process.env.CHROME_CHANNEL||'chrome'});
  const context=await browser.newContext({acceptDownloads:true});
  const page=await context.newPage(),errors=[];
  page.on('pageerror',error=>errors.push(String(error)));
  const base=process.env.UI_FIXTURE_URL;
  async function select(kind){
    await page.locator('#moduleSelect').selectOption(kind);
    await page.locator('#scopeAck').check();
  }
  async function save(){
    const response=page.waitForResponse(r=>r.url().endsWith('/api/modules/jobs')&&r.request().method()==='POST');
    await page.locator('#saveModule').click();const r=await response;
    assert.equal(r.status(),201,await r.text());const body=await r.json();
    await page.waitForFunction(id=>document.getElementById('moduleJobIdentity').textContent.includes(id),body.id);
    await page.waitForFunction(()=>!document.getElementById('moduleStart').disabled);
    return body.id;
  }
  async function execute(expected='completed'){
    await page.locator('#moduleStart').click();
    await page.waitForFunction(()=>/run: (completed|blocked|failed)|build: (completed|failed)/.test(document.getElementById('moduleJobState').textContent),null,{timeout:180000});
    assert.match(await page.locator('#moduleJobState').textContent(),new RegExp(expected),await page.locator('#moduleOperations').textContent());
    assert.equal(await page.locator('#moduleError').textContent(),'');
  }
  try{
    await page.goto(base+'modules.html');
    await page.waitForFunction(()=>document.getElementById('moduleSelect').options.length===11);
    // Every independent form must be accessible with zero research runs.
    for(const kind of ['prepare','peptide2mol','reinvent4','risk','availability','docking','docking_control','screening','interaction','md','md_tpr']){
      await select(kind);assert.ok((await page.locator('#moduleDescription').textContent()).length>20);
      assert.equal(await page.locator('#moduleFields [data-key="executable"]').count(),0);
    }
    await select('risk');
    await page.locator('[data-role="library"]').setInputFiles({name:'ethanol.smi',mimeType:'text/plain',buffer:Buffer.from('id\tsmiles\nPUBLIC_ETHANOL\tCCO\n')});
    await save();await execute();
    assert.match(await page.locator('#moduleOperations').textContent(),/"risk"/);
    const riskLink=page.locator('#moduleArtifacts a').filter({hasText:'risk/summary.json'}).first();
    const download=page.waitForEvent('download');await riskLink.click();const artifact=await download;
    assert.equal(await artifact.failure(),null);
    const saved=await artifact.path();assert.ok(JSON.parse(fs.readFileSync(saved,'utf8')).input===1);
    await select('prepare');
    await page.locator('[data-role="complex"]').setInputFiles(path.join(root,'data/structures/9nu4.pdb'));
    await save();await execute();
    await select('interaction');
    await page.locator('[data-role="complex"]').setInputFiles(path.join(root,'04_interaction_analysis/tests/data/3PTB.pdb'));
    await page.locator('[data-key="ligand"]').fill('BEN');
    await save();await execute();
    if(process.env.MODULE_UI_GMX==='1'){
      await select('md_tpr');
      await page.locator('[data-role="production"]').setInputFiles(path.join(process.env.UI_FIXTURE_ROOT,'md_fixture/production.tpr'));
      await page.locator('[data-key="device"]').selectOption('cpu');await page.locator('[data-key="threads"]').fill('1');
      await page.locator('[data-key="review_reason"]').fill('Public two-argon engineering fixture; no biological claim.');
      await save();await execute();
      const link=page.locator('#moduleArtifacts a').filter({hasText:'prepared_md/completion.json'}).first();
      const waiting=page.waitForEvent('download');await link.click();const d=await waiting;
      const proof=JSON.parse(fs.readFileSync(await d.path(),'utf8'));
      assert.equal(proof.plan.final_step,100);assert.equal(proof.candidate_qualified,false);
    }
    await page.screenshot({path:process.env.UI_SCREENSHOT||'module-ui.png',fullPage:true});
    await page.setViewportSize({width:390,height:844});
    assert.equal(await page.evaluate(()=>document.documentElement.scrollWidth>window.innerWidth+1),false);
    assert.deepEqual(errors,[]);
    console.log('PASS: 11 independent forms; real risk, preparation, interaction'+(process.env.MODULE_UI_GMX==='1'?', GROMACS completion':'')+'; artifact download; mobile layout; no browser errors');
  }finally{await context.close();await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
