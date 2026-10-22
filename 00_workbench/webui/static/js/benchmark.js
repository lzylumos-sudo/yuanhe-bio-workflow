const node=(tag,text,cls)=>{const e=document.createElement(tag);if(text!==undefined)e.textContent=text;if(cls)e.className=cls;return e};
const panel=document.querySelector('#benchmarkPanel');
let latest,choice='ppi_benchmark',arm='A1';
function table(headers,rows){const wrap=node('div',undefined,'table-wrap'),t=node('table'),head=node('thead'),tr=node('tr'),body=node('tbody');headers.forEach(x=>tr.append(node('th',x)));head.append(tr);rows.forEach(row=>{const r=node('tr');row.forEach(x=>r.append(node('td',x??'—')));body.append(r)});t.append(head,body);wrap.append(t);return wrap}
function link(name,label){const a=node('a',label||name,'subtle');a.href='/api/benchmark/file?name='+encodeURIComponent(name);a.target='_blank';a.rel='noopener';return a}
function render(d){
 panel.hidden=!d.configured;if(!d.configured)return;
 panel.replaceChildren();const title=node('div',undefined,'panel-title');title.append(node('h2','MDM2–p53 能力验证'),node('span',d.results_ready?'结果已保存':'真实CPU筛选进行中','pill'));panel.append(title,node('p',d.scope));
 panel.append(node('p',`完整库 ${d.population??'—'} 个母体 · 构象准备失败 ${d.preparation_failures} 个 · 已有合格结果复用 ${d.reused_jobs} 个`));
 const p=d.progress,bar=node('progress');bar.max=p.planned||1;bar.value=p.completed||0;bar.style.width='100%';panel.append(bar,node('p',`已处理 ${p.completed||0} / ${p.planned||0} 个可计算分子；成功 ${p.success||0} 个。总耗时 ${d.elapsed_minutes.toFixed(1)} / ${d.hard_budget_minutes} 分钟。`));
 panel.append(node('p',`科学版本 ${d.scientific_git_commit?.slice(0,12)||'未冻结'} · ${d.run_id}`,'muted'));
 if(!d.results_ready){panel.append(node('p','计算与评价尚未完成，界面不提前显示活性标签和参考分子排名。每5秒更新。'));return}
 const controls=node('div',undefined,'actions'),profile=node('select'),select=node('select');profile.setAttribute('aria-label','验证过滤配置');select.setAttribute('aria-label','验证实验臂');
 for(const [key,label] of [['ppi_benchmark','PPI研究配置'],['current_product','现有产品规则']])profile.append(new Option(label,key));profile.value=choice;profile.onchange=()=>{choice=profile.value;render(latest)};
 for(const [key,label] of [['A0','A0 无肽基线'],['A1','A1 真实p53接触特征'],['A2','A2 扰动特征对照']])select.append(new Option(label,key));select.value=arm;select.onchange=()=>{arm=select.value;render(latest)};controls.append(profile,select);panel.append(controls);
 const m=d.metrics,stats=m.profiles[choice][arm];panel.append(node('h3',`${choice} / ${arm}`));
 panel.append(table(['Top-k','活性命中 / 全部活性','Recall','EF'],[10,20,50].map(k=>{const x=stats.top_k[k];return[k,`${x.active_hits} / ${m.actives}`,(100*x.Recall).toFixed(1)+'%',x.EF.toFixed(2)]})));
 panel.append(table(['参考分子','排名','首次淘汰阶段'],d.references.filter(r=>r.profile===choice&&r.arm===arm).map(r=>[r.name,r.rank||'无排名',r.first_exclusion_stage||'保留'])));
 panel.append(node('p',`真实对接完成 ${m.docked} / ${m.total_parents}。所有指标仍以完整库和全部 ${m.actives} 个活性分子为分母。产品规则与研究结果分别报告。`));
 if(!m.complete_docking)panel.append(node('p','存在未完成计算的分子，不能宣称完整库全部通过验证。','error'));
 const image=node('img');image.src='/api/benchmark/file?name=recovery_metrics.png';image.alt='三组回收率和富集因子比较';image.style.cssText='max-width:100%;height:auto';panel.append(image);
 const downloads=node('div',undefined,'actions');for(const [name,label]of [['validation_report.md','中文验证报告'],['scientific_top50.csv','研究Top-50（PPI/A1）'],['deliverable_top50.csv','可交付Top-50'],['all_rankings.csv','完整六组排名'],['molecule_stage_ledger.csv','逐阶段记录'],['reference_recovery.csv','参考分子回收']])downloads.append(link(name,label));panel.append(downloads);
}
async function update(){try{const r=await fetch('/api/benchmark/status');if(!r.ok)throw Error('任务状态读取失败');latest=await r.json();render(latest)}catch(e){panel.hidden=false;panel.replaceChildren(node('h2','MDM2验证状态暂不可用'),node('p',e.message));}}
update();setInterval(update,5000);
