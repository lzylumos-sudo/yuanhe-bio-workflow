// Run against evidence_ui_fixture.py only. This is UI acceptance, not model acceptance.
const {chromium} = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const root = process.env.UI_FIXTURE_ROOT;
const url = process.env.UI_FIXTURE_URL || 'http://127.0.0.1:8879/';
assert(root && new URL(url).hostname === '127.0.0.1', 'Use isolated local fixtures only');
(async () => {
  const browser = await chromium.launch({channel: process.env.CHROME_CHANNEL || 'chrome', headless: true});
  let page;
  try {
    page = await browser.newPage({viewport: {width: 1400, height: 1000}});
    const errors = [], posts = [];
    page.on('pageerror', error => errors.push(error.message));
    page.on('request', request => {if(request.method() === 'POST') posts.push(request.url());});
    await page.goto(url); await page.getByText('Service connected', {exact:true}).waitFor();
    // Selection starts an asynchronous run load. Wait until its final request
    // completes so the click cannot hit the empty-state "create run" button.
    const loaded = page.waitForResponse(r => r.url().endsWith('/api/project/run/UI_RISK_FIXTURE/structure'));
    await page.locator('#runs').selectOption('UI_RISK_FIXTURE'); await loaded;
    await page.getByRole('button', {name:/Substructure risk screening/}).click();
    await page.getByText('· 1 records', {exact:true}).waitFor();
    assert.equal(await page.locator('#error').innerText(), '');
    const records = path.join(root, 'UI_RISK_FIXTURE/risk/records.json');
    const before = fs.readFileSync(records, 'utf8');
    const open = () => page.getByRole('button', {name:'Revise pre-docking evidence (preserve records)', exact:true}).click();
    await open(); await page.locator('#evidenceRevisionDialog').waitFor();
    assert.equal(fs.readFileSync(records, 'utf8'), before); assert.equal(posts.length, 0);
    const firstDialog = page.locator('#evidenceRevisionDialog');
    await firstDialog.getByText('Add an individual evidence review', {exact:true}).click();
    await firstDialog.locator('[name=reviewer]').fill('UI_TEST_ONLY');
    await firstDialog.locator('[name=candidate]').selectOption('SYNTHETIC_ONLY');
    await firstDialog.locator('[name=review_target]').selectOption('BRENK:Aliphatic_long_chain');
    await firstDialog.locator('[name=reviewed_at]').fill('2026-10-28');
    await firstDialog.locator('[name=review_notes]').fill('Synthetic UI interaction; will cancel, not a scientific approval');
    await firstDialog.getByRole('button', {name:'Add review record', exact:true}).click();
    await firstDialog.getByText('1 review record(s) ready to submit.', {exact:true}).waitFor();
    await page.getByRole('button', {name:'Cancel', exact:true}).click();
    assert.equal(fs.readFileSync(records, 'utf8'), before); assert.equal(posts.length, 0);
    await open();
    const dialog = page.locator('#evidenceRevisionDialog');
    await dialog.locator('[name=reviewer]').fill('UI_TEST_ONLY');
    await dialog.locator('[name=reason]').fill('No scientific approval; regression of evidence editor');
    await dialog.getByText('Advanced review records (JSON)', {exact:true}).click();
    await dialog.locator('[name=risk_reviews]').fill('{invalid');
    await dialog.getByRole('button', {name:'Save evidence revision', exact:true}).click();
    await page.locator('#evidenceRevisionError').filter({hasText:/JSON/}).waitFor();
    assert.equal(posts.length, 0); assert.equal(fs.readFileSync(records, 'utf8'), before);
    await dialog.locator('[name=risk_reviews]').fill('[]');
    await dialog.getByRole('button', {name:'Save evidence revision', exact:true}).click();
    await dialog.waitFor({state:'detached'});
    await page.getByRole('button', {name:/Substructure risk screening.*Not run/s}).waitFor();
    assert.equal(posts.length, 1); assert(!fs.existsSync(records));
    const history = path.join(root, 'UI_RISK_FIXTURE/evidence_revisions');
    const revision = fs.readdirSync(history)[0];
    assert.equal(fs.readFileSync(path.join(history, revision, 'risk/records.json'),'utf8'), before);
    const original = JSON.parse(fs.readFileSync(path.join(root,'UI_RISK_FIXTURE/config.json')));
    const reviewed = JSON.parse(fs.readFileSync(path.join(root,'UI_RISK_FIXTURE/reviews/risk.json')));
    assert.deepEqual(original.predocking, reviewed.predocking);
    const switched = page.waitForResponse(r => r.url().endsWith('/api/project/run/UI_ROUTE_FIXTURE/structure'));
    await page.locator('#runs').selectOption('UI_ROUTE_FIXTURE'); await switched;
    await page.getByRole('button', {name:/Pre-Docking Screening/}).click();
    await page.getByText('· 1 records', {exact:true}).waitFor();
    await open();
    const routeDialog = page.locator('#evidenceRevisionDialog');
    await routeDialog.locator('[name=reviewer]').fill('UI_TEST_ONLY');
    await routeDialog.locator('[name=reason]').fill('Incomplete route must stay pending');
    await routeDialog.getByText('Add an individual evidence review', {exact:true}).click();
    await routeDialog.locator('[name=candidate]').selectOption('SYNTHETIC_ONLY');
    assert.equal(await routeDialog.locator('[name=review_target] option').nth(1).evaluate(option=>option.disabled), true);
    await routeDialog.getByText('Advanced review records (JSON)', {exact:true}).click();
    const routeRecords = path.join(root, 'UI_ROUTE_FIXTURE/availability/records.json');
    const routeBefore = fs.readFileSync(routeRecords, 'utf8'), row = JSON.parse(routeBefore)[0];
    await routeDialog.locator('[name=route_reviews]').fill(JSON.stringify([{
      compound_id: row.compound_id, route_sha256: row.availability.routes[0].route_sha256,
      reviewer:'UI_TEST_ONLY', reason:'Synthetic rejection test', reviewed_at:'2026-10-28', identity_stereo_selectivity_confirmed:true,
    }]));
    await routeDialog.getByRole('button', {name:'Save evidence revision', exact:true}).click();
    await page.locator('#evidenceRevisionError').filter({hasText:'cannot replace'}).waitFor();
    assert.equal(fs.readFileSync(routeRecords, 'utf8'), routeBefore);
    assert.deepEqual(errors, []);
    if (process.env.UI_SCREENSHOT) await page.screenshot({path:process.env.UI_SCREENSHOT});
    await routeDialog.getByRole('button', {name:'Cancel', exact:true}).click();
    console.log('PASS: large-field table, read-only open/cancel, invalid JSON, evidence-only save, preserved archives/policies, incomplete route rejected; no scientific jobs started.');
  } catch (error) {
    if (page && process.env.UI_SCREENSHOT) {
      await page.screenshot({path:process.env.UI_SCREENSHOT}).catch(()=>{});
      console.error('UI failure details:', await page.locator('#error').innerText().catch(()=>''));
    }
    throw error;
  } finally {await browser.close();}
})().catch(error=>{console.error(error);process.exitCode=1;});
