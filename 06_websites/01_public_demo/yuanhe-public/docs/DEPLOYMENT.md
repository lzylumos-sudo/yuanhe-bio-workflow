# Cloudflare Pages deployment

Production URL: https://yuanhe-bio-demo.pages.dev/
Project: `yuanhe-bio-demo`; production branch: `main`.

Run from `06_websites/01_public_demo/yuanhe-public` with Node.js 22.12+ and pinned pnpm:

```sh
pnpm install --frozen-lockfile
pnpm build
pnpm test
pnpm exec playwright test
pnpm dlx wrangler pages deploy dist --project-name yuanhe-bio-demo --branch main
```

Use the existing Cloudflare account. Upload only `dist`, never the repository root or private workbench. The Pages dashboard also accepts a ZIP containing the files of `dist` at its root. Updating GitHub alone does not publish this Direct Upload project.

The release contains public 9NU4 coordinates and derived contact geometry. It does not execute scientific workloads. Preserve the sample-data and evidence-boundary labels.

After publication, check the dark homepage, Chinese/English switching, all routes, 3D structure loading, contact filtering and CSV/JSON export. Verify HTTPS and security headers. Use the previous production deployment for rollback if needed.

The 2026-10-20 release was prepared on Windows; its local audit hashes were independently verified before publication from macOS. Browser acceptance recorded on Windows is separate from subsequent live-site verification.
