// Browser regression for UI state only. Synthetic API fixtures do not certify engines.
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
const assert=require('node:assert/strict'),http=require('node:http'),fs=require('node:fs'),path=require('node:path');
const root=path.resolve(__dirname,'../..'),staticRoot=path.join(root,'00_workbench/webui/static');
const version=fs.readFileSync(path.join(root,'VERSION'),'utf8').trim();
const makeJob=(id,module,report)=>({id,module,parameters:{},inputs:{},scope:'Synthetic UI regression only',systems:[],artifacts:[],operations:report?[{action:'run',status:'completed',result:{result:{report},candidate_qualified:false},log:''}]:[]});
const jobs={A:makeJob('A','risk'),B:makeJob('B','peptide2mol',{completed:1,valid:0,unique:0}),C:makeJob('C','reinvent4',{completed:1,valid:1,unique:1})};
jobs.E={...makeJob('E','md'),systems:[{id:'MD_fixture',status:'building'}],operations:[{id:'build-E',action:'build',status:'running',log:''}]};
let serviceVersion=version,available=true,startMode='busy',startCount=0,saveCount=0,holdRead=null,holdStart=null;
const catalog=Object.fromEntries(['risk','peptide2mol','reinvent4','md'].map(k=>[k,{label:k,description:'Synthetic public UI fixture; no chemical computation.',files:{},fields:[]}]))
const json=(res,data,status=200)=>{res.writeHead(status,{'Content-Type':'application/json'});res.end(JSON.stringify(data));};
const server=http.createServer(async(req,res)=>{
  const url=new URL(req.url,'http://127.0.0.1');
  if(url.pathname==='/api/v1/versions')return json(res,{software_version:serviceVersion,source:{commit:'0123456789abcdef0123456789abcdef01234567'},git:{dirty:false}});
  if(url.pathname==='/api/project')return json(res,{version:{software_version:serviceVersion,source:{commit:'0123456789abcdef0123456789abcdef01234567'}},notice:'UI regression fixture',registered_modules:[],runs:[],defaults:{structure:'public.pdb',receptor_range:[3,500],receptor_chain:'A',peptide_chain:'B',budget:1,seed:17,predocking:{}}});
  if(url.pathname==='/api/modules/jobs'){
    if(req.method==='POST'){saveCount++;jobs.D=makeJob('D','risk');return json(res,jobs.D,201);}
    return json(res,{modules:catalog,jobs:Object.values(jobs).map(({id,module})=>({id,module})),resource_status:{available,policy:'exclusive',queued:false}});
  }
  const match=url.pathname.match(/^\/api\/modules\/jobs\/([^/]+)(\/(?:start|cancel))?$/);
  if(match?.[2]==='/cancel'){
    const job=jobs[match[1]];job.operations[0].status='cancelled';json(res,job);
    setTimeout(()=>{job.systems[0].status='cancelled';},200);return;
  }
  if(match){const id=match[1];if(match[2]){startCount++;if(startMode==='late'){holdStart=res;return;}available=false;return json(res,{error:'Resource lock occupied',error_code:'RESOURCE_BUSY',recoverable:true,details:{queued:false}},409);}if(holdRead===id){holdRead=res;return;}return json(res,jobs[id]);}
  let file=path.resolve(staticRoot,'.'+(url.pathname==='/'?'/index.html':url.pathname));
  if(!file.startsWith(staticRoot+path.sep)||!fs.existsSync(file)){res.writeHead(404);return res.end();}
  res.setHeader('Content-Type',file.endsWith('.js')?'text/javascript':file.endsWith('.css')?'text/css':'text/html');res.end(fs.readFileSync(file));
});
const until=async(test)=>{for(let n=0;n<100;n++){if(test())return;await new Promise(r=>setTimeout(r,20));}throw Error('Fixture did not receive expected request');};
(async()=>{
 await new Promise(r=>server.listen(0,'127.0.0.1',r));const base='http://127.0.0.1:'+server.address().port;
 const browser=await chromium.launch({headless:true,...(process.env.CHROME_EXECUTABLE?{executablePath:process.env.CHROME_EXECUTABLE}:{channel:process.env.CHROME_CHANNEL||'chrome'})});
 const page=await browser.newPage(),errors=[];page.on('pageerror',e=>errors.push(e.message));
 async function selected(id){await page.locator('#moduleJobSelect').selectOption(id);await page.waitForFunction(id=>document.getElementById('moduleJobIdentity').textContent.endsWith(' · '+id),id);}
 try{
  await page.goto(base+'/modules.html');await page.waitForFunction(()=>!document.getElementById('saveModule').disabled);
  assert.match(await page.locator('#modulePlatformVersion').innerText(),/0123456789ab/);
  await selected('A');await page.locator('#moduleStart').click();await page.waitForFunction(()=>document.getElementById('moduleError').textContent.includes('not been queued'));
  assert.equal(await page.locator('#moduleError').getAttribute('role'),'status');assert.equal(await page.locator('#moduleStart').isDisabled(),true);
  await page.waitForTimeout(3200);assert.equal(startCount,1,'Busy jobs must not automatically resubmit');
  available=true;await page.locator('#refreshModules').click();await page.waitForFunction(()=>!document.getElementById('moduleStart').disabled);assert.equal(await page.locator('#moduleError').innerText(),'');
  await selected('B');assert.match(await page.locator('#moduleOperations').innerText(),/Valid molecules: 0.*Unique molecules: 0/);assert.match(await page.locator('#moduleOperations').innerText(),/No usable molecules/);
  await selected('C');assert.match(await page.locator('#moduleOperations').innerText(),/does not establish candidate qualification/);
  // A delayed response for a previous selection may not replace the new job.
  holdRead='A';await page.locator('#moduleJobSelect').selectOption('A');await until(()=>holdRead&&typeof holdRead!=='string');await selected('B');const lateRead=holdRead;holdRead=null;json(lateRead,jobs.A);await page.waitForTimeout(100);assert.match(await page.locator('#moduleJobIdentity').innerText(),/ · B$/);
  // The same protection applies to errors from an earlier start request.
  await selected('A');startMode='late';await page.locator('#moduleStart').click();await until(()=>holdStart);await selected('B');json(holdStart,{error:'Old job start failed'},500);holdStart=null;await page.waitForTimeout(150);assert.equal(await page.locator('#moduleError').innerText(),'');assert.match(await page.locator('#moduleJobIdentity').innerText(),/ · B$/);
  // Successful input creation selects the new immutable snapshot.
  await page.locator('#scopeAck').check();await page.locator('#saveModule').click();await page.waitForFunction(()=>document.getElementById('moduleJobIdentity').textContent.endsWith(' · D'));assert.equal(saveCount,1);
  // A cancelled operation can arrive before the worker's final system metadata.
  await selected('E');await page.locator('#moduleCancel').click();
  await page.waitForFunction(()=>document.getElementById('moduleJobState').textContent==='build: cancelled');
  await page.waitForFunction(()=>document.getElementById('moduleSystem').textContent.includes('cancelled'),null,{timeout:10000});
  assert.equal(await page.locator('#moduleRunMD').isDisabled(),true);
  await selected('D');
  serviceVersion='0.0.0-mismatched';const previousStarts=startCount;await page.locator('#moduleStart').click();await page.waitForFunction(()=>document.getElementById('modulePlatformVersion').textContent.includes('differs'));assert.equal(startCount,previousStarts,'A service upgrade while the page is open must block mutation before POST');
  assert.equal(await page.locator('#saveModule').isDisabled(),true);assert.equal(await page.locator('#moduleStart').isDisabled(),true);await selected('B');assert.match(await page.locator('#moduleOperations').innerText(),/No usable molecules/);
  // Main workbench also refuses new work with mismatched frontend/backend versions.
  await page.goto(base+'/');await page.waitForFunction(()=>document.getElementById('connection').textContent==='Version mismatch');await page.locator('#newRun').click();assert.equal(await page.locator('#runForm button[type="submit"]').isDisabled(),true);await page.locator('#closeDialog').click();
  assert.deepEqual(errors,[]);console.log('PASS: recoverable busy state; no auto queue; stale response/error isolation; cancellation metadata refresh; zero outputs; version identity and mismatch guards; historical viewing');
 }finally{await browser.close();await new Promise(r=>server.close(r));}
})().catch(e=>{console.error(e);server.close();process.exitCode=1;});
