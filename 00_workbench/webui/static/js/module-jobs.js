import {PLATFORM_VERSION} from './platform-version.js';
const $=id=>document.getElementById(id);
const element=(tag,text)=>{const e=document.createElement(tag);if(text!==undefined)e.textContent=text;return e;};
let catalog={},current=null,busy=false,resourceStatus=null,versionCompatible=false;
let selectionEpoch=0,detailRequest=0,listRequest=0,errorVersion=0;
const resourceNote=element('p'),versionNote=element('p','Checking service version…');
resourceNote.id='moduleResourceStatus';resourceNote.setAttribute('role','status');
versionNote.id='modulePlatformVersion';versionNote.setAttribute('role','status');
$('moduleError').after(resourceNote);$('moduleDescription').before(versionNote);
async function api(path,body){
  if(body!==undefined&&!path.endsWith('/cancel')){await checkVersion();if(!versionCompatible)throw Error('Service version is not verified. Refresh the page before starting or saving a job.');}
  const r=await fetch(path,body===undefined?{}:{method:'POST',headers:{'Content-Type':'application/json','X-Yuanhe-UI-Version':PLATFORM_VERSION},body:JSON.stringify(body)}),data=await r.json();
  if(!r.ok){const error=Error(data.error||data.message||('Request failed ('+r.status+')'));Object.assign(error,{code:data.error_code,recoverable:data.recoverable,details:data.details,suggestedAction:data.suggested_action});throw error;}return data;
}
function clearError(expected=errorVersion){if(expected!==errorVersion)return;$('moduleError').textContent='';$('moduleError').removeAttribute('data-code');errorVersion++;}
function fail(e,epoch=selectionEpoch,expected=errorVersion){
  if(epoch!==selectionEpoch||expected!==errorVersion)return;
  if(e.code==='UI_VERSION_MISMATCH'){versionCompatible=false;versionNote.textContent='Service version changed. Reload this page or deploy matching frontend and backend files before submitting jobs.';}
  const waiting=e.code==='RESOURCE_BUSY';
  $('moduleError').textContent=waiting?'Server busy. This request has not been queued. Wait for the active task to finish, then start again.':e.message||String(e);
  $('moduleError').dataset.code=e.code||'';$('moduleError').setAttribute('role',waiting?'status':'alert');errorVersion++;
  if(waiting){resourceStatus={available:false};render();}
}
function renderAvailability(){
  const occupied=resourceStatus?.available===false;
  resourceNote.textContent=occupied?'Compute resource occupied. No new job will start automatically. Refresh or wait for availability, then start the job.':'';
  $('saveModule').disabled=busy||!versionCompatible;
}
async function checkVersion(){
  try{const version=await api('/api/v1/versions'),commit=version.source?.commit||version.git?.commit||'unknown';versionCompatible=version.software_version===PLATFORM_VERSION;
    versionNote.textContent='Platform '+version.software_version+' · Source '+commit.slice(0,12)+(version.git?.dirty?' · Local changes':'')+(versionCompatible?'':' · UI version '+PLATFORM_VERSION+' differs. Reload the page or deploy matching frontend and backend files before submitting jobs.');
    versionNote.title='Server source: '+commit+'; UI version: '+PLATFORM_VERSION;
  }catch{versionCompatible=false;versionNote.textContent='Service version unavailable. Refresh before submitting or starting jobs.';}render();
}
function form(){
  const spec=catalog[$('moduleSelect').value];if(!spec)return;
  $('moduleDescription').textContent=spec.description;
  $('moduleFiles').replaceChildren();$('moduleFields').replaceChildren();$('scopeAck').checked=false;
  for(const [role,types] of Object.entries(spec.files)){
    const label=element('label',role+' ('+types.replace('?','')+(types.endsWith('?')?', optional':'')+')');
    const file=element('input');file.type='file';file.accept=types.replace('?','');file.required=!types.endsWith('?');file.dataset.role=role;file.setAttribute('aria-label',role+' input');
    label.append(file);$('moduleFiles').append(label);
  }
  for(const f of spec.fields){
    const label=element('label',f.label),input=element(f.type==='choice'?'select':'input');input.dataset.key=f.key;input.required=true;
    if(f.type==='choice')for(const choice of f.choices){const option=element('option',choice);option.value=choice;input.append(option);}
    else{input.type=['integer','number'].includes(f.type)?'number':'text';if(input.type==='number'){input.min=f.min;input.max=f.max;input.step=f.type==='integer'?'1':'any';}}
    input.value=f.type==='vector'?f.default.join(', '):f.default;input.setAttribute('aria-label',f.label);label.append(input);$('moduleFields').append(label);
  }
}
async function refresh({automatic=false}={}){
  const request=++listRequest,epoch=selectionEpoch,previousError=errorVersion;
  let info;try{info=await api('/api/modules/jobs');}catch(e){if(request===listRequest)fail(e,epoch,previousError);return;}
  if(request!==listRequest||epoch!==selectionEpoch)return;catalog=info.modules;resourceStatus=info.resource_status||null;
  if(!$('moduleSelect').options.length){for(const [key,spec] of Object.entries(catalog)){const option=element('option',spec.label);option.value=key;$('moduleSelect').append(option);}form();}
  const selected=current?.id||$('moduleJobSelect').value;$('moduleJobSelect').replaceChildren(element('option','Select a saved job'));$('moduleJobSelect').options[0].value='';
  for(const job of info.jobs){const option=element('option',(catalog[job.module]?.label||job.module)+' · '+job.id);option.value=job.id;$('moduleJobSelect').append(option);}
  $('moduleJobSelect').value=selected;
  if(selected)await load(selected,{clearFailure:!automatic||$('moduleError').dataset.code==='RESOURCE_BUSY'&&resourceStatus?.available===true});
  else if(!automatic)clearError(previousError);render();
}
async function load(id,{clearFailure=true}={}){
  const epoch=selectionEpoch,request=++detailRequest,previousError=errorVersion;
  try{const job=await api('/api/modules/jobs/'+id);if(epoch!==selectionEpoch||request!==detailRequest||$('moduleJobSelect').value!==id)return;current=job;render();if(clearFailure)clearError(previousError);}
  catch(e){if(request===detailRequest)fail(e,epoch,previousError);}
}
function render(){
  renderAvailability();
  const job=current;if(!job)return;$('moduleJobDetail').hidden=false;
  $('moduleJobIdentity').textContent=(catalog[job.module]?.label||job.module)+' · '+job.id;
  $('moduleJobState').textContent=job.operations[0]?job.operations[0].action+': '+job.operations[0].status:'Input saved; computation has not started';
  $('moduleInputs').textContent=JSON.stringify({parameters:job.parameters,inputs:job.inputs,scope:job.scope},null,2);
  const running=job.operations.some(x=>x.status==='running'),md=job.module==='md',last=job.operations[0];
  $('moduleStart').textContent=md?'Build membrane system':last?.status==='completed'?'Completed':job.module==='md_tpr'&&last?'Resume prepared MD':'Start module';
  $('moduleStart').disabled=!versionCompatible||resourceStatus?.available===false||running||busy||(job.operations.length>0&&(md||job.module!=='md_tpr'||last.status==='completed'));
  $('moduleCancel').disabled=!running||busy;
  $('moduleSystemReview').hidden=!md||!job.systems.length;
  const previous=$('moduleSystem').value;$('moduleSystem').replaceChildren();
  for(const system of job.systems){const option=element('option',system.id+' · '+system.status+' · '+(system.atom_count??'?')+' atoms');option.value=system.id;$('moduleSystem').append(option);}
  if([...$('moduleSystem').options].some(x=>x.value===previous))$('moduleSystem').value=previous;
  const system=job.systems.find(x=>x.id===$('moduleSystem').value);
  const reviewed=$('moduleReviewAck').checked&&$('moduleReviewReason').value.trim();
  $('moduleRunMD').hidden=!md;$('moduleAnalyzeMD').hidden=!md;
  $('moduleRunMD').disabled=!versionCompatible||resourceStatus?.available===false||running||busy||!reviewed||!system||!['built','completed'].includes(system.status)||job.operations.some(x=>x.action==='run'&&x.status==='completed');
  $('moduleAnalyzeMD').disabled=!versionCompatible||resourceStatus?.available===false||running||busy||!reviewed||system?.status!=='completed'||job.operations.some(x=>x.action==='analyze'&&x.status==='completed');
  const operationRows=new Map([...$('moduleOperations').children].map(row=>[row.dataset.operation,row])),activeOperations=new Set();
  for(const [index,op] of job.operations.entries()){
    const key=job.id+':'+(op.id??index),previous=operationRows.get(key),signature=JSON.stringify([op.action,op.status,op.error,op.result]);activeOperations.add(key);
    const logText=op.log||'Engine logs are available in the artifact list.';
    if(previous?._signature===signature){const log=previous.querySelector('details pre');if(log.textContent!==logText)log.textContent=logText;if($('moduleOperations').children[index]!==previous)$('moduleOperations').insertBefore(previous,$('moduleOperations').children[index]||null);continue;}
    const opened=previous?[...previous.querySelectorAll('details')].map(details=>details.open):null;
    const row=element('div');row.className='operation';row.dataset.operation=key;row._signature=signature;row.append(element('strong',op.action+': '+op.status));if(op.error)row.append(element('p',op.error));
    const result=op.result?.result;
    if(result?.plan)row.append(element('p','Final step: '+result.plan.final_step+' · Final time: '+result.plan.final_time_ps+' ps · Checkpoint and trajectory verified'));
    else if(result?.runs)row.append(element('p',result.runs.length+' reference runs · Top-1 pose recovery: '+Math.round(result.top1_recovery_fraction*100)+'%'));
    else if(result?.ligands)row.append(element('p',result.ligands.length+' ligand analyses · '+result.ligands.reduce((n,r)=>n+r.n_events,0)+' interaction events'));
    else if(result?.status==='blocked_by_strict_preflight')row.append(element('p','Strict pre-docking screening did not pass. No docking was started.'));
    else if(result?.report){const report=result.report;row.append(element('p','Actual sampling attempts: '+report.completed+' · Valid molecules: '+(report.valid??'unknown')+' · Unique molecules: '+(report.unique??'unknown')));if(report.valid===0)row.append(element('p','No usable molecules were produced. Execution finished; output acceptance has not passed.'));else row.append(element('p','Generated molecules require screening and review; execution completion does not establish candidate qualification.'));}
    else if(result?.risk)row.append(element('p','Risk screening finished for '+result.risk.input+' input records. Inspect detailed decisions below.'));
    const details=element('details');details.open=opened?opened[0]:['running','failed','blocked'].includes(op.status);details.append(element('summary','Execution log'),element('pre',logText));row.append(details);
    const proof=element('details');proof.open=opened?opened[1]:false;proof.append(element('summary','Detailed result and evidence'),element('pre',op.result?JSON.stringify(op.result,null,2):'No completed result'));row.append(proof);if(previous)previous.replaceWith(row);else $('moduleOperations').append(row);if($('moduleOperations').children[index]!==row)$('moduleOperations').insertBefore(row,$('moduleOperations').children[index]||null);}
  for(const [key,row] of operationRows)if(!activeOperations.has(key))row.remove();
  const artifactLinks=new Map([...$('moduleArtifacts').children].map(link=>[link.getAttribute('href'),link])),activeArtifacts=new Set();
  for(const file of job.artifacts){const href='/module-artifacts/'+job.id+'/'+file.path.split('/').map(encodeURIComponent).join('/'),label=file.path+' ('+file.size+' bytes)';activeArtifacts.add(href);let link=artifactLinks.get(href);if(!link){link=element('a');link.href=href;link.setAttribute('download','');$('moduleArtifacts').append(link);}if(link.textContent!==label)link.textContent=label;}
  for(const [href,link] of artifactLinks)if(!activeArtifacts.has(href))link.remove();
}
function encoded(file){return new Promise((resolve,reject)=>{const reader=new FileReader();reader.onerror=()=>reject(Error('Could not read input file'));reader.onload=()=>resolve(reader.result.split(',')[1]);reader.readAsDataURL(file);});}
$('moduleForm').onsubmit=async event=>{
  event.preventDefault();if(busy||!versionCompatible)return;busy=true;clearError();const epoch=selectionEpoch,previousError=errorVersion;renderAvailability();
  try{
    const kind=$('moduleSelect').value,spec=catalog[kind],parameters={},files={};
    for(const f of spec.fields){const value=$('moduleFields').querySelector('[data-key="'+f.key+'"]').value;
      parameters[f.key]=f.type==='vector'?value.split(',').map(x=>Number(x.trim())):['number','integer'].includes(f.type)?Number(value):value;}
    for(const input of $('moduleFiles').querySelectorAll('input')){const file=input.files[0];if(!file)continue;if(file.size>64*1024*1024)throw Error('Each file must be <= 64 MiB');files[input.dataset.role]={name:file.name,base64:await encoded(file)};}
    const saved=await api('/api/modules/jobs',{module:kind,parameters,files,acknowledged:$('scopeAck').checked});if(epoch!==selectionEpoch)return;selectionEpoch++;current=saved;await refresh();
  }catch(e){fail(e,epoch,previousError);}finally{busy=false;render();}
};
async function start(action){
  if(!current||busy||!versionCompatible||resourceStatus?.available===false)return;busy=true;clearError();const id=current.id,epoch=selectionEpoch,previousError=errorVersion;render();
  try{const job=await api('/api/modules/jobs/'+id+'/start',{action,system_id:$('moduleSystem').value,system_reviewed:$('moduleReviewAck').checked,review_reason:$('moduleReviewReason').value});if(epoch===selectionEpoch&&$('moduleJobSelect').value===id)current=job;}
  catch(e){fail(e,epoch,previousError);}finally{busy=false;if(epoch===selectionEpoch&&$('moduleJobSelect').value===id)await load(id,{clearFailure:false});render();}
}
$('moduleStart').onclick=()=>start(current?.module==='md'?'build':'run');
$('moduleRunMD').onclick=()=>start('run');$('moduleAnalyzeMD').onclick=()=>start('analyze');
$('moduleCancel').onclick=async()=>{if(!current||busy)return;const id=current.id,epoch=selectionEpoch;clearError();const previousError=errorVersion;busy=true;render();try{const job=await api('/api/modules/jobs/'+id+'/cancel',{});if(epoch===selectionEpoch&&$('moduleJobSelect').value===id)current=job;}catch(e){fail(e,epoch,previousError);}finally{busy=false;render();}};
$('moduleSelect').onchange=()=>{clearError();form();};$('refreshModules').onclick=async()=>{await checkVersion();await refresh();};
$('moduleJobSelect').onchange=()=>{selectionEpoch++;detailRequest++;clearError();const id=$('moduleJobSelect').value;current=null;$('moduleJobDetail').hidden=true;if(id){$('moduleReviewAck').checked=false;$('moduleReviewReason').value='';load(id);}renderAvailability();};
$('moduleReviewAck').onchange=render;$('moduleReviewReason').oninput=render;$('moduleSystem').onchange=()=>{$('moduleReviewAck').checked=false;render();};
renderAvailability();checkVersion();refresh();
setInterval(()=>{if(busy)return;if(resourceStatus?.available===false)refresh({automatic:true});else if(current&&current.operations.some(x=>x.status==='running'))load(current.id);},3000);
