# 构建与真实验收记录

日期：2026-10-13。项目 F:\02_docking；专用 WSL GNINA-SBVS-1e0ea3f4；UUID 1e0ea3f4-3fef-4108-ae52-4af02d120ed8。

## 阶段 1：独立环境与启动器 — completed
Ubuntu Base 24.04.5 官方 rootfs 与 SHA256 匹配。WSL2 实际导入/启动成功，VHDX 位于项目 runtime/wsl/ext4.vhdx，F 盘 exFAT。现有默认 Ubuntu-22.04 未修改，最终列表显示其仍 Stopped。
证据：logs/install/initial_setup.log、env/rootfs.lock.json、runtime/local_machine.json、reports/register_acceptance.log。

## 阶段 2：Vina、化学准备与 CCO — passed
独立 Vina 1.2.7，Python 3.11.17、RDKit 2026.9.1、Meeko 0.8.0。1DRF 的唯一 FOL 实例 A:187，共 32 重原子，CCD 恢复键级/电荷并保留晶体坐标。初次受体准备因 GLN12 缺原子失败，随后明确使用 PDBFixer 补缺失原子，不重建 loop，不删除坏残基。受体实际 Meeko PDBQT 及准备 JSON 均保留。
命令：run.sh prepare-receptor / make-box / screen；证据：logs/install/demo_install.log、reports/cco_cli_smoke.json、data/prepared/demo_1drf/preparation.json。

## 阶段 3：默认页面 — passed
可见 label 配体 SMILES 与按钮 开始筛选；单独 CCO 可直接运行，无额外确认。最新三次全新 run 耗时 1.973 / 1.891 / 1.902 秒，真实 Vina affinity -2.546，页面/CSV/原始输出一致。首次页面无预填结果。
命令：tests/e2e/browser_acceptance.py；证据：reports/browser_acceptance.json 与 reports/browser/。

## 阶段 4：完整课堂、批量、图与报告 — passed
20/20 演示分子成功，180 poses；100 分子流程库已提供。9 步页面读实际文件，四个二级页面真实浏览器渲染无 JS 错误。两张不同实际 pose 图已生成并检查。
运行：runs/20261013_104206_d78a5f9f32；报告为该目录 report.html；证据：reports/class_report_images.json、reports/pages_acceptance.json。

## 阶段 5：GNINA 与科研两阶段 — passed
GNINA 1.3.3 官方 SHA256 匹配。初次缺 CUDA 用户态库与无效模型别名的失败日志保留；最终固定 crossdock_default2018_1_3，分别完成 GPU 和 --no_gpu CPU CNN 对接。GPU trace 有 768 次成功 /dev/dxg ioctl。独立 Vina 与 GNINA 未混用。
命令：scripts/gnina_smoke.py gpu / cpu；run.sh rescore --operation score_only / redock / refinement。证据：reports/gnina_*、reports/gpu_device_trace.json、reports/research_acceptance.json。多 seed 11/23/47 已实际运行。

## 阶段 6：科研审查、QC、状态与导出 — passed with documented scope
公开自定义受体的计划/显式确认/准备/口袋预览/CCO 筛选已实测。支持 mmCIF 计划和外部 PDBQT 来源 unknown。Molscrub 0.2.2 pH 状态枚举与双 seed 状态/母体汇总实测通过。Vina 导出 180 poses 化学身份一致；FOL 同坐标系对称性 RMSD Top1 1.8458 Å。人工 MD 导出门禁、SDF 坐标和文件可解析性通过自动测试样本验证，未执行 MD。
证据：reports/artifact_acceptance.json、reports/enumeration_acceptance.json、reports/FOL_redocking.json。

## 阶段 7：恢复、离线、路径、文档与备份 — passed within tested scope
真实取消、坏分子混合批次、服务硬中断后恢复、已完成条目不重算通过。损坏输出校验失败后，故障条目不误计成功，真实重试恢复；CLI 占用时停机保护拒绝继续。
Linux 无网络命名空间 CCO、浏览器禁止外网、本地 3D、Windows 参数数组转发、中文空格目录迁移路径通过。同机现有 VHDX 登记路径/UUID 验证、doctor/CCO、安全停止及重启通过。
独立 runtime 快照 runtime/downloads/runtime_snapshot_20261013_191355.tar，13,558,968,320 字节，SHA256 清单已生成。未测试快照恢复导入、复制 VHDX 后新注册或另一台真实电脑。
证据：reports/recovery_acceptance.json、corruption_acceptance.json、stop_guard.json、offline_namespace.log、path_acceptance.json、windows_path_test.log、runtime_export.json。

## 真实限制
老师地址在宿主限时请求返回 HTTP 401（UNAUTHENTICATED），未看到案例内容；初次沙箱请求 HTTP 000 不是最终连通性结论。未做远程老师平台验收。
未开放共价对接、一般金属/复杂辅因子参数化、柔性受体、缺失长 loop 重建及可靠 MOL2 输入。未做标签富集验证。CPU --no_gpu 在当前 GPU 主机实测，另一台物理无 GPU 机器未实测。
完整矩阵和路径见 docs/实现与验收报告.md 与 reports/acceptance_matrix.json。not_tested 和 blocked 未计入 passed。

## 当前可继续状态
Web 服务由 sbvs 普通用户运行于 127.0.0.1:5902，Windows localhost health 已核验项目 UUID、ready=true。环境备份完成后已重新启动。进一步开发应保留已完成 run 和原始分数，先读实际验收证据再修改；不得自动注销发行版。
