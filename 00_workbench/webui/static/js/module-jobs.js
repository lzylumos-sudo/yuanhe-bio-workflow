const $=id=>document.getElementById(id);
const element=(tag,text)=>{const e=document.createElement(tag);if(text!==undefined)e.textContent=text;return e;};
let catalog={},current=null,busy=false;
async function api(path,body){const r=await fetch(path,body===undefined?{}:{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});const data=await r.json();if(!r.ok)throw Error(data.error||data.message||('Request failed ('+r.status+')'));return data;}
function fail(e){$('moduleError').textContent=e.message||String(e);}
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
async function refresh(){
  const info=await api('/api/modules/jobs');catalog=info.modules;
  if(!$('moduleSelect').options.length){for(const [key,spec] of Object.entries(catalog)){const option=element('option',spec.label);option.value=key;$('moduleSelect').append(option);}form();}
  const selected=current?.id||$('moduleJobSelect').value;$('moduleJobSelect').replaceChildren(element('option','Select a saved job'));$('moduleJobSelect').options[0].value='';
  for(const job of info.jobs){const option=element('option',(catalog[job.module]?.label||job.module)+' · '+job.id);option.value=job.id;$('moduleJobSelect').append(option);}
  $('moduleJobSelect').value=selected;if(selected)await load(selected);
}
async function load(id){current=await api('/api/modules/jobs/'+id);render();}
function render(){
  const job=current;if(!job)return;$('moduleJobDetail').hidden=false;
  $('moduleJobIdentity').textContent=(catalog[job.module]?.label||job.module)+' · '+job.id;
  $('moduleJobState').textContent=job.operations[0]?job.operations[0].action+': '+job.operations[0].status:'Input saved; computation has not started';
  $('moduleInputs').textContent=JSON.stringify({parameters:job.parameters,inputs:job.inputs,scope:job.scope},null,2);
  const running=job.operations.some(x=>x.status==='running'),md=job.module==='md',last=job.operations[0];
  $('moduleStart').textContent=md?'Build membrane system':job.module==='md_tpr'&&last?'Resume prepared MD':'Start module';
  $('moduleStart').disabled=running||busy||(job.operations.length>0&&(md||job.module!=='md_tpr'||last.status==='completed'));
  $('moduleCancel').disabled=!running||busy;
  $('moduleSystemReview').hidden=!md||!job.systems.length;
  const previous=$('moduleSystem').value;$('moduleSystem').replaceChildren();
  for(const system of job.systems){const option=element('option',system.id+' · '+system.status+' · '+(system.atom_count??'?')+' atoms');option.value=system.id;$('moduleSystem').append(option);}
  if([...$('moduleSystem').options].some(x=>x.value===previous))$('moduleSystem').value=previous;
  const system=job.systems.find(x=>x.id===$('moduleSystem').value);
  const reviewed=$('moduleReviewAck').checked&&$('moduleReviewReason').value.trim();
  $('moduleRunMD').hidden=!md;$('moduleAnalyzeMD').hidden=!md;
  $('moduleRunMD').disabled=running||busy||!reviewed||!system||!['built','completed'].includes(system.status)||job.operations.some(x=>x.action==='run'&&x.status==='completed');
  $('moduleAnalyzeMD').disabled=running||busy||!reviewed||system?.status!=='completed'||job.operations.some(x=>x.action==='analyze'&&x.status==='completed');
  $('moduleOperations').replaceChildren();
  for(const op of job.operations){const row=element('div');row.className='operation';row.append(element('strong',op.action+': '+op.status));if(op.error)row.append(element('p',op.error));
    const details=element('details');details.open=op===last;details.append(element('summary','Logs and result'),element('pre',op.log||'No log output yet'),element('pre',op.result?JSON.stringify(op.result,null,2):'No completed result'));row.append(details);$('moduleOperations').append(row);}
  $('moduleArtifacts').replaceChildren();
  for(const file of job.artifacts){const link=element('a',file.path+' ('+file.size+' bytes)');link.href='/module-artifacts/'+job.id+'/'+file.path.split('/').map(encodeURIComponent).join('/');link.setAttribute('download','');$('moduleArtifacts').append(link);}
}
function encoded(file){return new Promise((resolve,reject)=>{const reader=new FileReader();reader.onerror=()=>reject(Error('Could not read input file'));reader.onload=()=>resolve(reader.result.split(',')[1]);reader.readAsDataURL(file);});}
$('moduleForm').onsubmit=async event=>{
  event.preventDefault();if(busy)return;busy=true;$('saveModule').disabled=true;$('moduleError').textContent='';
  try{
    const kind=$('moduleSelect').value,spec=catalog[kind],parameters={},files={};
    for(const f of spec.fields){const value=$('moduleFields').querySelector('[data-key="'+f.key+'"]').value;
      parameters[f.key]=f.type==='vector'?value.split(',').map(x=>Number(x.trim())):['number','integer'].includes(f.type)?Number(value):value;}
    for(const input of $('moduleFiles').querySelectorAll('input')){const file=input.files[0];if(!file)continue;if(file.size>64*1024*1024)throw Error('Each file must be <= 64 MiB');files[input.dataset.role]={name:file.name,base64:await encoded(file)};}
    current=await api('/api/modules/jobs',{module:kind,parameters,files,acknowledged:$('scopeAck').checked});await refresh();
  }catch(e){fail(e);}finally{busy=false;$('saveModule').disabled=false;render();}
};
async function start(action){
  if(!current||busy)return;busy=true;render();$('moduleError').textContent='';
  try{current=await api('/api/modules/jobs/'+current.id+'/start',{action,system_id:$('moduleSystem').value,system_reviewed:$('moduleReviewAck').checked,review_reason:$('moduleReviewReason').value});}
  catch(e){fail(e);}finally{busy=false;if(current)await load(current.id).catch(fail);}
}
$('moduleStart').onclick=()=>start(current?.module==='md'?'build':'run');
$('moduleRunMD').onclick=()=>start('run');$('moduleAnalyzeMD').onclick=()=>start('analyze');
$('moduleCancel').onclick=async()=>{try{current=await api('/api/modules/jobs/'+current.id+'/cancel',{});render();}catch(e){fail(e);}};
$('moduleSelect').onchange=form;$('refreshModules').onclick=()=>refresh().catch(fail);
$('moduleJobSelect').onchange=()=>{const id=$('moduleJobSelect').value;if(id){$('moduleReviewAck').checked=false;$('moduleReviewReason').value='';load(id).catch(fail);}else{current=null;$('moduleJobDetail').hidden=true;}};
$('moduleReviewAck').onchange=render;$('moduleReviewReason').oninput=render;$('moduleSystem').onchange=()=>{$('moduleReviewAck').checked=false;render();};
refresh().catch(fail);
setInterval(()=>{if(current&&!busy&&current.operations.some(x=>x.status==='running'))load(current.id).catch(fail);},3000);
