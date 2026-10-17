# GNINA / Vina 本地 SBVS 工作站

双击项目根目录的 `START_WINDOWS.cmd`，打开 `http://localhost:5902/`。默认输入框的可见标签为“配体 SMILES”，输入 `CCO` 后点击“开始筛选”，即运行独立 Vina CPU。初次页面没有历史分数；所有分数均来自新建 run 的真实引擎输出。

本项目采用 Windows 外层目录 + 项目专用 WSL2 ext4.vhdx。当前项目身份见 `project.json`，实际 WSL 登记见 `runtime/local_machine.json`。只用专用 `/opt/sbvs/` 环境，不借用用户其它 Ubuntu。

## 入口

- 课堂快捷：单个或多行 SMILES；也可上传 SMI/TXT/CSV/SDF。
- 9 步课堂：真实环境、受体来源、准备记录、盒子、库文件，随后使用同一后端对接、排序、看姿态、保存图和报告。
- 科研：受体 PDB/mmCIF/PDBQT 上传、准备计划与明确审查、盒子预览确认、Vina/GNINA 协议、自定义搜索参数、状态枚举。
- 历史任务：独立的 run 目录；运行结果及日志不覆盖。
- 环境与导出：真实版本与设备；实际验收结果详见 `docs/实现与验收报告.md`。

关闭请双击 `STOP_WINDOWS.cmd`。有任务时先等待完成，或明确取消；脚本不会杀死其它项目进程。迁移前执行 `PREPARE_MIGRATION.ps1`，等待确认本项目已停止，再自行检查磁盘上其它程序的占用。新电脑/新盘符上执行 `REGISTER_ON_THIS_PC.ps1`，随后正常启动。不要在 WSL 运行时复制 VHDX，不要使用 `wsl --unregister` 解绑。

完整说明：`docs/使用指南.md`、`docs/科研筛选指南.md`、`docs/迁移与恢复.md`。本项目不自动运行 MD，交付接口只导出人工选择的文件。
