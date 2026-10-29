const {spawn}=require('node:child_process'),fs=require('node:fs'),os=require('node:os'),path=require('node:path');
const root=fs.mkdtempSync(path.join(process.env.UI_TEST_TEMP||os.tmpdir(),'yuanhe-module-ui-'));
const args=[path.join(__dirname,'module_ui_fixture.py'),'--root',root,'--port','0'];
if(process.env.MODULE_TEST_GMX)args.push('--gmx',process.env.MODULE_TEST_GMX);
const server=spawn(process.env.PYTHON_EXECUTABLE||'python',args,{windowsHide:true,stdio:['ignore','pipe','pipe']});
let log='',started=false,browser;
const timer=setTimeout(()=>{console.error(log);process.exitCode=1;server.kill();},60000);
server.stderr.on('data',x=>{log+=x.toString();});
server.stdout.on('data',x=>{
 log+=x.toString();const match=log.match(/FIXTURE_URL=(http:\/\/127\.0\.0\.1:\d+\/)/);if(!match||started)return;
 started=true;clearTimeout(timer);
 browser=spawn(process.execPath,[path.join(__dirname,'test_module_ui.cjs')],{windowsHide:true,stdio:'inherit',
  env:{...process.env,UI_FIXTURE_ROOT:root,UI_FIXTURE_URL:match[1],MODULE_UI_GMX:process.env.MODULE_TEST_GMX?'1':'0'}});
 browser.on('exit',code=>{process.exitCode=code??1;if(code)console.error(log);server.kill();});
 browser.on('error',e=>{console.error(e);process.exitCode=1;server.kill();});
});
server.on('error',e=>{clearTimeout(timer);console.error(e);process.exitCode=1;});
server.on('exit',()=>{clearTimeout(timer);if(!started){console.error(log);process.exitCode=1;}});
for(const signal of ['SIGTERM','SIGINT'])process.on(signal,()=>{browser?.kill();server.kill();process.exitCode=1;});
