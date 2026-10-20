# 元核生物 · 公开商业演示网站

独立 React 19 + TypeScript + Vite 静态网站，NGL 按需加载。无科研后端、登录、GPU、移动盘或私有 API 依赖。六个页面、双语界面、八阶段演示流程。

## 本地运行

需要 Node.js 22.12+（本次 Node 24 验证）、pnpm 11.19.0。

```sh
pnpm install --frozen-lockfile
pnpm dev
pnpm build
pnpm preview --port 4173
pnpm test
pnpm exec playwright test
```

浏览器测试默认使用已安装的 Chrome；如使用 Playwright Chromium，安装其浏览器并在 playwright.config.ts 移除 channel。测试配置会自动在 14173 端口启动构建预览，测试前先运行 `pnpm build`。

## 交付范围与真实状态

- 首页、产品、交互工作台、MurJ 案例、技术研究、合作页面，中文 / English 切换。
- 公开 RCSB 9NU4 原始 PDB、深色首页与结构渲染图；NGL 支持旋转、缩放、链显示、界面查看、残基定位。
- 由公开结构独立预计算的接触表：筛选、搜索、排序、CSV 下载；4 Å 下 19 个 A 链残基。
- 八阶段均可浏览；生成、优化、对接、筛选、MD 没有获准公开样例，明确标识占位，不展示伪造分子、分数、轨迹或计算进度。
- 研究简报支持浏览器打印 / 保存 PDF、JSON 证据下载。不是候选药物验证报告。
- 合作页可下载需求草稿。项目方尚未提供正式联系邮箱，因此没有伪造联系方式或发送接口；这项上线内容仍待补齐。

生产交付为 `dist/`，部署方式见 [docs/DEPLOYMENT.md](docs/DEPLOYMENT.md)。演示步骤见 [docs/DEMO.md](docs/DEMO.md)，数据边界见 [docs/DATA.md](docs/DATA.md)，仓库检查见 [docs/AUDIT.md](docs/AUDIT.md)。
