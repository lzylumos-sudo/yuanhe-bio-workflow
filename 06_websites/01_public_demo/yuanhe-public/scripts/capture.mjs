import {chromium} from '@playwright/test';
import fs from 'node:fs';
fs.mkdirSync('../../work/qa',{recursive:true});
const browser=await chromium.launch({channel:'chrome',headless:true,args:['--enable-unsafe-swiftshader']});
const page=await browser.newPage({viewport:{width:1440,height:1000}});
await page.goto('http://127.0.0.1:4173');await page.screenshot({path:'../../work/qa/home-desktop.png',fullPage:true});
await page.goto('http://127.0.0.1:4173/#demo');await page.getByRole('button',{name:/加载公开结构/}).click();await page.locator('.ngl canvas').waitFor();await page.locator('.viewer-message').waitFor({state:'detached'});await page.waitForTimeout(1800);await page.screenshot({path:'../../work/qa/demo-desktop.png',fullPage:true});
await page.setViewportSize({width:390,height:844});await page.goto('http://127.0.0.1:4173');await page.screenshot({path:'../../work/qa/home-mobile.png',fullPage:true});
await browser.close();
