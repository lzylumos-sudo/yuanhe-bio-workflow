// Own both fixture and browser lifetimes; never leave a local server running.
const {spawn} = require('node:child_process');
const fs = require('node:fs'), os = require('node:os'), path = require('node:path');
const directory = fs.mkdtempSync(path.join(process.env.UI_TEST_TEMP || os.tmpdir(), 'yuanhe-ui-evidence-'));
const root = path.join(directory, 'runs');
const server = spawn(process.env.PYTHON_EXECUTABLE || 'python',
  [path.join(__dirname, 'evidence_ui_fixture.py'), '--root', root, '--port', '0'],
  {windowsHide:true, stdio:['ignore','pipe','pipe']});
let log = '', started = false, browser;
const timer = setTimeout(()=>{console.error('Fixture startup timed out:',log);process.exitCode=1;server.kill();},30000);
server.stderr.on('data', data=>{log+=data.toString();});
server.stdout.on('data', data=>{
  log+=data.toString();const match=log.match(/FIXTURE_URL=(http:\/\/127\.0\.0\.1:\d+\/)/);
  if(!match || started)return;started=true;clearTimeout(timer);
  browser=spawn(process.execPath,[path.join(__dirname,'test_evidence_revision_ui.cjs')],
    {windowsHide:true,stdio:'inherit',env:{...process.env,UI_FIXTURE_ROOT:root,UI_FIXTURE_URL:match[1]}});
  browser.on('exit',code=>{process.exitCode=code??1;server.kill();});
  browser.on('error',error=>{console.error(error);process.exitCode=1;server.kill();});
});
server.on('error',error=>{clearTimeout(timer);console.error(error);process.exitCode=1;});
server.on('exit',()=>{clearTimeout(timer);if(!started){console.error(log);process.exitCode=1;}});
process.on('SIGINT',()=>{browser?.kill();server.kill();process.exitCode=130;});
process.on('SIGTERM',()=>{browser?.kill();server.kill();process.exitCode=143;});
