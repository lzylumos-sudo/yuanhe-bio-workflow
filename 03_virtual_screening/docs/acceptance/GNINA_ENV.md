# GNINA 真实执行环境：已接通到哪一步（本轮实测）

## 结论

**GNINA 仍不能真实执行**，阻塞点是 CUDA 运行时库链，不是本项目代码、不是路径、也不是权限。
因此本轮**没有任何真实 GNINA 复对接证据**；`redock` 阶段在本机输出 `blocked_external_tool`，
正式主候选数为 0（缺 seed 支持）。适配器本身已按规格实现完毕（多 seed、固定受体坐标系聚类、
planned/attempted/succeeded/supporting 四类计数分离），但没有真实运行过，状态记为**未验证**。

## 已核实的事实（可复核）

| 项 | 值 |
| --- | --- |
| Linux 入口 | WSL2 `Ubuntu-22.04`（`Ubuntu 22.04.5 LTS`），`wsl.exe -d Ubuntu-22.04 -u root` |
| 上游项目 Linux 环境 | `F:\MemPilot` 记录为 `/opt/sbvs/env`（Linux glibc 2.39），本机该目录**不存在**，无法复用 |
| GPU | `NVIDIA GeForce RTX 5070`，`/usr/lib/wsl/lib` 中存在 `libcuda.so.1`、`libnvidia-*`、`nvidia-smi` |
| GNINA 二进制 | 按上游 `env/vendor.lock.json` 锁定版本下载：1.3.3 `gnina.cuda12.8.static`，2,056,131,000 字节，sha256 `3340c1f4…f72172`（与上游记录一致） |
| 位置 | `06_Virtual_Screening_Studio/vendor/gnina`，通过 `/mnt/c/...` 在 WSL 内执行 |

## 运行时依赖补齐过程（逐项实测）

| 顺序 | 报错 | 处理 | 结果 |
| --- | --- | --- | --- |
| 1 | `libcudnn.so.9: cannot open shared object file` | 从 NVIDIA redist CDN 取 cuDNN 9.8.0.87（`cudnn-linux-x86_64-9.8.0.87_cuda12-archive.tar.xz`，sha256 记录在 `vendor/cudnn_runtime.json`），解包到 `vendor/cudnn/lib` | 通过 |
| 2 | `libcudart.so.12` 缺失 | 从 PyPI `nvidia-cuda-runtime-cu12` wheel 抽取 `libcudart.so.12` 到 `vendor/nvidia-cuda-runtime-cu12/lib`（sha256 记录在 `vendor/cuda_runtime.json`） | 通过 |
| 3 | `libcublas.so.12` 缺失 | 同上抽取 `nvidia-cublas-cu12` 的 `libcublas(-Lt).so.12` | 通过 |
| 4 | **`libcusparse.so.12` 缺失** | 未处理（需要继续补齐 cuSPARSE、cuRAND、cuFFT、cuSOLVER 等整套 CUDA 12 运行时） | **阻塞** |

`gnina.cuda12.8.static` 并非完全静态链接，它仍动态依赖整套 CUDA 12 运行时；
本机 WSL 发行版只带驱动库，不带 CUDA 运行时。

## 若要解除阻塞（具体动作，不含“放宽规则”）

1. 在 `Ubuntu-22.04` 内安装 CUDA 12 运行时（`cuda-cudart-12-*`、`libcusparse-12-*`、`libcurand-12-*`、
   `libcufft-12-*`、`libcusolver-12-*`），或把上游机器上 `/opt/sbvs` 整目录复制到本机 WSL；
2. 用 `wsl -d Ubuntu-22.04 -u root -e env LD_LIBRARY_PATH=… vendor/gnina --version` 复验，应输出 `gnina v1.3.3`；
3. 再执行 `zvss run --budget configs/budgets/redock_trial.yaml`，产出 5 化合物 × 3 seed 的真实复对接证据。

在这三步完成前，`redock` 的证据字段一律为 `not_evaluated` / `needs_input`，正式主候选保持 0。
