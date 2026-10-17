# 第三方来源与归属

本仓库包含项目编排、适配器与分析代码；不随仓库分发第三方模型源码、权重或计算环境。

| 组件 | 原始来源 | 本项目使用方式 |
|---|---|---|
| Peptide2Mol | https://github.com/BLUE-Flowing/Peptide2Mol | 独立模型子进程；本地来源带 MIT 许可证 |
| REINVENT4 / LibInvent | https://github.com/MolecularAI/REINVENT4 | 独立采样；本地来源带 Apache-2.0 许可证 |
| GNINA | https://github.com/gnina/gnina | 外部对接程序 |
| RDKit | https://www.rdkit.org | 化学身份、构象与结构处理 |
| OpenMM | https://openmm.org | 独立动力学后端 |
| GROMACS | https://www.gromacs.org | 已准备体系执行 |
| ProLIF / PyMOL / PoseBusters | 见各模块环境及来源文档 | 相互作用分析、可视化、质量检查 |

具体提交与权重 SHA-256 见 `deployment/model_dependencies.json`。引用链接用于来源识别；本次未重新联网审计全部许可条款。恢复第三方目录时必须保留其许可证与作者声明。

原论文全文与抽取文本不再分发。结构示例保留原始 PDB 字节，来源与校验见 `data/structures/manifest.json`。它们不是团队测定的结构，也不是实验药效证据。
