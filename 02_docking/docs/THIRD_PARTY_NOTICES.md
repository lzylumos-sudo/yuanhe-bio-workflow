# 第三方来源与许可证

本项目不把整套发行包统一标为 Apache-2.0。GNINA 包含其仓库说明的 GPL/Apache 授权部分及第三方依赖；准确原文保存在 licenses/，分发前应分别审查这些条款。NVIDIA 用户态库有独立许可，不等于显卡驱动随项目迁移。

源码许可证原文与 URL/SHA256：docs/licenses/sources.json。Conda 每个实际包的版本、build、来源和 license：env/package_provenance.json。pip 来源：env/pip_sources.json。引擎与前端资产：env/vendor.lock.json。

Ubuntu Base 来源 https://cdimage.ubuntu.com/ubuntu-base/releases/24.04/release/ ，原始发行包保存在 runtime/downloads/。RCSB PDB/CCD 来源与哈希在 data/sources.json。3Dmol.js 2.5.5 通过固定版本 jsDelivr npm 镜像本地化，不在浏览器运行时访问 CDN。
