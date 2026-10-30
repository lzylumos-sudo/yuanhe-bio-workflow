// Delayed synthetic API responses must never replace the user's current view.
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
const assert=require('node:assert/strict'),http=require('node:http'),fs=require('node:fs'),path=require('node:path');
const root=path.resolve(__dirname,'../..'),staticRoot=path.join(root,'00_workbench/webui/static');
const version=fs.readFileSync(path.join(root,'VERSION'),'utf8').trim();
const config={structure:'synthetic.pdb',receptor_chain:'A',peptide_chain:'B',receptor_range:[1,2],budget:1,seed:17,predocking:{}};
const stageNames={prepare:'Complex preparation',peptide2mol:'Peptide2Mol',reinvent4:'REINVENT4',risk:'Substructure risk screening',availability:'Pre-Docking Screening',docking:'AutoDock Vina',screening:'Virtual screening',gromacs:'GROMACS'};
const makeRun=id=>({id,config,reviews:{},stages:Object.entries(stageNames).map(([key,label])=>({key,label,status:key==='risk'?'pending':'completed',counts:{},description:'Synthetic UI fixture',artifacts:{}}))});
const holds=new Set(),pending=new Map();
const json=(res,data,status=200)=>{res.writeHead(status,{'Content-Type':'application/json'});res.end(JSON.stringify(data));};
function respond(key,res,data){if(holds.delete(key)){pending.set(key,{res,data});return;}json(res,data);}
function release(key,error=false){const held=pending.get(key);assert(held,'Missing delayed response '+key);pending.delete(key);json(held.res,error?{error:'Stale request failure'}:held.data,error?500:200);}
const until=async(test)=>{for(let n=0;n<100;n++){if(test())return;await new Promise(r=>setTimeout(r,20));}throw Error('Expected delayed request was not received');};
const server=http.createServer((req,res)=>{
 const url=new URL(req.url,'http://127.0.0.1'),p=url.pathname;
 if(p==='/api/project')return json(res,{version:{software_version:version,git:{commit:'synthetic'}},notice:'Synthetic UI state test; no execution',registered_modules:[],runs:['A','B'],defaults:config});
 if(p==='/api/v1/versions')return json(res,{software_version:version,git:{commit:'synthetic'}});
 const match=p.match(/^\/api\/project\/run\/([AB])(?:\/(rows|log)\/([^/]+)|\/(structure))?$/);
 if(match){const [,id,kind,key,structure]=match;if(structure)return json(res,{atoms:[]});if(!kind)return respond('detail:'+id,res,makeRun(id));
  if(kind==='rows'){const offset=url.searchParams.get('offset')||'0',search=url.searchParams.get('search')||'';return respond('rows:'+id+':'+key+':'+offset+':'+search,res,{rows:[{id:[id,key,offset,search||'all'].join('-')}],total:250});}
  return respond('log:'+id+':'+key,res,{text:id+' '+key+' log'});
 }
 if(/^\/api\/studio\/run\/[AB]\/operations$/.test(p))return json(res,{operations:[]});
 let file=path.resolve(staticRoot,'.'+(p==='/'?'/index.html':p));
 if(!file.startsWith(staticRoot+path.sep)||!fs.existsSync(file)){res.writeHead(404);return res.end();}
 if(p==='/js/studio.js'&&process.env.STUDIO_SCRIPT)file=process.env.STUDIO_SCRIPT;
 res.setHeader('Content-Type',file.endsWith('.js')?'text/javascript':file.endsWith('.css')?'text/css':file.endsWith('.svg')?'image/svg+xml':'text/html');res.end(fs.readFileSync(file));
});
(async()=>{
 await new Promise(r=>server.listen(0,'127.0.0.1',r));
 const browser=await chromium.launch({headless:true,...(process.env.CHROME_EXECUTABLE?{executablePath:process.env.CHROME_EXECUTABLE}:{channel:process.env.CHROME_CHANNEL||'chrome'})});
 const page=await browser.newPage(),errors=[],posts=[];page.on('pageerror',e=>errors.push(e.message));page.on('request',r=>{if(r.method()==='POST')posts.push(r.url());});
 const table=()=>page.locator('#tableBody').innerText();
 async function select(id){await page.locator('#runs').selectOption(id);await page.waitForFunction(id=>document.querySelector('#runs').value===id&&document.querySelector('#tableBody').textContent.includes(id+'-'),id);}
 async function stage(key){await page.getByRole('button',{name:new RegExp(stageNames[key])}).first().click();}
 try{
  await page.goto('http://127.0.0.1:'+server.address().port+'/');await page.getByText('Service connected',{exact:true}).waitFor();await select('A');
  // A new selection stops exposing A as an actionable run while B is loading.
  holds.add('detail:B');await page.locator('#runs').selectOption('B');await until(()=>pending.has('detail:B'));
  assert.equal(await page.locator('#execute').isDisabled(),true);
  const oldStage=page.getByRole('button',{name:new RegExp(stageNames.risk)}).first();assert.equal(await oldStage.isDisabled(),true);
  await oldStage.dispatchEvent('click');await page.locator('#execute').dispatchEvent('click');
  assert.equal(await page.locator('#execute').isDisabled(),true);assert.deepEqual(posts,[]);
  await page.getByRole('button',{name:'① Prepare and review receptor / Docking box',exact:true}).click();
  assert.match(await page.locator('#error').innerText(),/Create a run first/);assert.deepEqual(posts,[]);
  release('detail:B');await page.waitForFunction(()=>document.querySelector('#tableBody').textContent.includes('B-risk'));
  await stage('prepare');await select('A');
  // A refresh of A may finish after the user has selected B.
  holds.add('detail:A');await page.locator('#refresh').click();await until(()=>pending.has('detail:A'));await select('B');release('detail:A');await page.waitForTimeout(100);
  assert.equal(await page.locator('#runs').inputValue(),'B');assert.match(await table(),/B-prepare/);
  // Late request failures from an old run must not appear against the new run.
  await select('A');holds.add('detail:A');await page.locator('#refresh').click();await until(()=>pending.has('detail:A'));await select('B');release('detail:A',true);await page.waitForTimeout(100);
  assert.equal(await page.locator('#error').innerText(),'');assert.equal(await page.locator('#runs').inputValue(),'B');
  // Changing stages protects both the record table and execution log.
  holds.add('rows:B:prepare:0:');await page.locator('#refresh').click();await until(()=>pending.has('rows:B:prepare:0:'));await stage('risk');await page.waitForFunction(()=>document.querySelector('#tableBody').textContent.includes('B-risk'));release('rows:B:prepare:0:');await page.waitForTimeout(100);assert.match(await table(),/B-risk/);
  holds.add('log:B:prepare');await stage('prepare');await until(()=>pending.has('log:B:prepare'));await stage('risk');await page.waitForFunction(()=>document.querySelector('#log').textContent==='B risk log');release('log:B:prepare');await page.waitForTimeout(100);assert.equal(await page.locator('#log').textContent(),'B risk log');
  // An old page must not replace a newer search or reset its pagination.
  holds.add('rows:B:risk:100:');await page.locator('#next').click();await until(()=>pending.has('rows:B:risk:100:'));await page.locator('#search').fill('filtered');await page.waitForFunction(()=>document.querySelector('#tableBody').textContent.includes('B-risk-0-filtered'));release('rows:B:risk:100:');await page.waitForTimeout(100);
  assert.equal(await page.locator('#search').inputValue(),'filtered');assert.match(await table(),/B-risk-0-filtered/);assert.match(await page.locator('#pageInfo').innerText(),/^1–/);assert.equal(await page.locator('#prev').isDisabled(),true);
  assert.deepEqual(errors,[]);if(process.env.UI_SCREENSHOT)await page.screenshot({path:process.env.UI_SCREENSHOT});
  console.log('PASS: main workbench run selection, stale errors, stage records/logs, pagination and search remain isolated from delayed responses');
 }finally{await browser.close();server.closeAllConnections();await new Promise(r=>server.close(r));}
})().catch(e=>{console.error(e);server.closeAllConnections();server.close();process.exitCode=1;});
