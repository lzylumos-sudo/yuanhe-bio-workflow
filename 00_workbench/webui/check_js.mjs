// Syntax-check every front-end ES module without spawning child processes
// (the sandbox blocks piped stdio, so we parse in-process with vm.SourceTextModule).
// Usage: node --experimental-vm-modules webui/check_js.mjs
import { readdirSync, readFileSync } from 'node:fs';
import { join, relative } from 'node:path';
import { fileURLToPath } from 'node:url';
import vm from 'node:vm';

const root = fileURLToPath(new URL('./static/js/', import.meta.url));
function directories(dir) {
  return [dir, ...readdirSync(dir, { withFileTypes: true })
    .filter((entry) => entry.isDirectory())
    .flatMap((entry) => directories(join(dir, entry.name)))];
}
const dirs = directories(root);
let failures = 0;
let checked = 0;

for (const dir of dirs) {
  let files;
  try { files = readdirSync(dir).filter((f) => f.endsWith('.js')); } catch { continue; }
  for (const file of files) {
    const path = join(dir, file);
    const src = readFileSync(path, 'utf8');
    checked += 1;
    try {
      // eslint-disable-next-line no-new
      new vm.SourceTextModule(src, { identifier: path });
      console.log(`ok   ${relative(root, path)}`);
    } catch (err) {
      failures += 1;
      console.log(`FAIL ${relative(root, path)}\n     ${String(err.message).split('\n')[0]}`);
    }
  }
}
console.log(`\n${checked - failures}/${checked} modules parsed`);
process.exit(failures ? 1 : 0);
