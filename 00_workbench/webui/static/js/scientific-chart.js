/* UI-only plotting: values and units come from the existing analysis response. */
export function scientificChart(canvas,rows,xKey,series,xLabel,yLabel){
 canvas._ybChartResize?.disconnect();
 const ctx=canvas.getContext('2d'),w=Math.max(280,Math.round(canvas.getBoundingClientRect().width||900)),h=280,left=64,right=20,top=w<500?90:60,bottom=50,dpr=window.devicePixelRatio||1;
 canvas.width=w*dpr;canvas.height=h*dpr;canvas.style.height=h+'px';ctx.setTransform(dpr,0,0,dpr,0,0);canvas.setAttribute('role','img');canvas.setAttribute('aria-label',`${yLabel} versus ${xLabel}. ${series.map(s=>s.label).join(', ')}`);canvas.tabIndex=0;
 canvas._ybChartResize=new ResizeObserver(()=>{const next=Math.max(280,Math.round(canvas.getBoundingClientRect().width));if(Math.abs(next-w)>1)scientificChart(canvas,rows,xKey,series,xLabel,yLabel)});canvas._ybChartResize.observe(canvas);
 let tip=canvas.nextElementSibling;if(!tip?.classList.contains('chart-tip')){tip=document.createElement('p');tip.className='chart-tip';canvas.after(tip)}
 const points=rows.filter(r=>Number.isFinite(r[xKey])&&series.some(s=>Number.isFinite(r[s.key])));
 const reset=()=>tip.textContent=points.length?'Hover, tap, or use ← → to inspect recorded values.':'No finite analysis values available.';
 reset();ctx.clearRect(0,0,w,h);ctx.font='12px Segoe UI';
 if(!points.length){ctx.fillStyle='#98a6b8';ctx.fillText('No analysis values available',left,120);return}
 let xmin=Math.min(...points.map(r=>r[xKey])),xmax=Math.max(...points.map(r=>r[xKey]));
 const ys=points.flatMap(r=>series.map(s=>r[s.key])).filter(Number.isFinite);let ymin=Math.min(...ys),ymax=Math.max(...ys);
 if(xmin===xmax){xmin-=.5;xmax+=.5}if(ymin===ymax){const pad=Math.abs(ymin)*.05||1;ymin-=pad;ymax+=pad}else{const pad=(ymax-ymin)*.08;ymin-=pad;ymax+=pad}
 const X=x=>left+(x-xmin)/(xmax-xmin)*(w-left-right),Y=y=>h-bottom-(y-ymin)/(ymax-ymin)*(h-top-bottom),fmt=v=>Math.abs(v)>=10000?v.toExponential(2):Number(v.toFixed(3)).toString();
 ctx.fillStyle='#98a6b8';ctx.fillText(yLabel,left,18);ctx.textAlign='center';
 for(let i=0;i<=4;i++){const y=ymin+(ymax-ymin)*i/4,px=left+(w-left-right)*i/4,py=Y(y);ctx.strokeStyle='#293545';ctx.beginPath();ctx.moveTo(left,py);ctx.lineTo(w-right,py);ctx.stroke();ctx.fillStyle='#98a6b8';ctx.textAlign='right';ctx.fillText(fmt(y),left-10,py+4);ctx.textAlign='center';ctx.fillText(fmt(xmin+(xmax-xmin)*i/4),px,h-bottom+20)}
 ctx.fillText(xLabel,(left+w-right)/2,h-9);ctx.textAlign='left';
 series.forEach((s,i)=>{ctx.strokeStyle=s.color;ctx.lineWidth=2;ctx.beginPath();let continuous=false;for(const r of points){if(!Number.isFinite(r[s.key])){continuous=false;continue}const x=X(r[xKey]),y=Y(r[s.key]);if(continuous)ctx.lineTo(x,y);else ctx.moveTo(x,y);continuous=true}ctx.stroke();for(const r of points){if(!Number.isFinite(r[s.key]))continue;ctx.fillStyle=s.color;ctx.beginPath();ctx.arc(X(r[xKey]),Y(r[s.key]),2,0,Math.PI*2);ctx.fill()}ctx.fillStyle=s.color;const lx=w<500?left:left+i*230,ly=w<500?34+i*18:34;ctx.fillRect(lx,ly,16,2);ctx.fillText(s.label,lx+23,ly+5)});
 let index=0;const show=i=>{index=Math.max(0,Math.min(points.length-1,i));const r=points[index];tip.textContent=`${xLabel}: ${fmt(r[xKey])} · `+series.map(s=>`${s.label}: ${Number.isFinite(r[s.key])?fmt(r[s.key]):'missing'} ${s.unit||''}`).join(' · ')};
 canvas.onpointermove=e=>{const b=canvas.getBoundingClientRect(),v=xmin+((e.clientX-b.left)/b.width*w-left)/(w-left-right)*(xmax-xmin);show(points.reduce((best,r,i)=>Math.abs(r[xKey]-v)<Math.abs(points[best][xKey]-v)?i:best,0))};
 canvas.onpointerleave=reset;canvas.onfocus=()=>show(index);canvas.onkeydown=e=>{if(['ArrowLeft','ArrowRight'].includes(e.key)){e.preventDefault();show(index+(e.key==='ArrowLeft'?-1:1))}};
}
