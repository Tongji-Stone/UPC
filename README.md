# UPC 2023 A：中文修订稿

全部权威代码、论文源码、输入副本和新结果均在本目录。根目录参考文献、UPC 和备份目录保持原样。

## 查看与复现

- `main.pdf`：最终中文论文；`main.tex`：可编辑主稿。
- `model.py`：七层大气、变马赫数阻力、显式面积、平移/姿态及有限充气求解器。
- `config.json`：质量、伞和典型设计条件；`run_analysis.py`：全部标定、检验和敏感性计算。
- `results.json` 和 `results/*.csv`：数值、逐时轨迹、误差、边界；`generated/*.tex`：从结果自动生成的正文数值及表格。
- `figures`：12 幅逻辑图，每幅含 PDF、SVG、320 DPI PNG；灰度预览在 `build/grayscale`。
- `data/source_registry.json`、`data/wyoming_soundings_provenance.json`：文献与探空来源；图线数字化说明见 `data/stratos_october_digitization.json`。

从本目录运行唯一计算入口：

```powershell
.\.venv\Scripts\python.exe reproduce.py
```

完整重算、重新制图、编译和论文一致性检查：

```powershell
.\build.ps1
```

只重新编译已更新的正文：` .\build.ps1 -SkipCalculations`。环境已含项目虚拟环境；新环境可用 Python 3.13 和 requirements.txt 安装依赖。计算无需联网，探空采用已保存的原始输入；`data/build_wyoming_soundings.py` 仅用于显式更新远端探空，不在日常复现中重新下载历史输入。

复现脚本调用已安装的 math-modeling skill 图表审计和清单工具，默认位置为用户目录 `.codex/skills/math-modeling`；编译需要本机 XeLaTeX、latexmk、pdfimages、pdftoppm。首次迁移计算环境时请同时保留这些只读工具或安装同一 skill。绘图辅助函数副本位于 utils。

## 关键物理口径

联合面积约 1.35697 m²，3/7 月标定固定 Cd=0.6；10 月只检验面积迁移。正式工况总质量 190 kg，Cd 采用 2017 式 26，M>1.25 延拓为 0.692。大气仅覆盖几何高度 0–86 km。

约 58.22 km 的答案依赖恢复温度 400 K、自由与开伞气动过载 5g、径向载荷 2g、着陆速度 6m/s 和中心姿态系数，属于典型装备设计情景。温度指标描述近壁气流；旋转系数和耐受筛选值是明确假设。10 月峰速低估 9.66%，完整曲线 RMSE 22.68m/s，保留在正文中。

## 构建与核对

`build_paper.py` 将 main.tex、generated 和正式 PDF 图件逐文件复制至 latex_project，再调用 skill 的真实 XeLaTeX 构建和验证。这样编译副本不包含 Python 虚拟环境。`main.build.json` 是编译器工具生成的源码/PDF 绑定；`build/revision/paper_source_binding.json` 将权威源码、图件与编译副本对应。`verify_translation.py` 已改为八章结构、禁用措辞、结果、边界、收敛和 PDF 源码一致性检查。

`results/复现清单.json` 保存输入 SHA-256、版本、参数、种子及唯一运行命令。旧稿和旧图只在 build/revision 归档，退出正式计算链。审计记录位于 build/revision；字体审计对 Type0 顶层 FontDescriptor 的“可能未嵌入”提示由 pdffonts 和递归字体检查复核，实际嵌入为 yes。



## 本次完成检查（2026-10-10）

最终PDF为19页：摘要第1页、目录第2页，正文从第3页开始；共11幅正式图、9张表、29个编号公式、8项参考文献。12幅候选图全部重新生成，另有灰度预览。最终PDF SHA-256：`0a41c7f6331213e462c9690fb2d5229cfa832a1f4e910696cd070cd4abf1c3fe`。

实际使用并读取了 math-modeling、建模手、编程手、论文手，以及 LaTeX、PDF、科研可视化、双引擎论文搜索入口。按用户要求保留既有中文LaTeX版式，仅交付PDF与源码。独立检查 M1、P1、P2、W1、W2 均通过；最终W2修正了2017论文卷期和一处图注说明，无遗留阻塞。

关键命令均已执行且退出码为0：

- `python test_model.py`：7项物理与数值测试。
- `python reproduce.py`：独立完整复算、12图和7张数据表；另两张正文表为参数与符号表。
- `build.ps1`：最终完整复现与编译；随后只对两处文字勘误执行 `build.ps1 -SkipCalculations`。
- `latex_paper.py doctor/build/validate`：XeLaTeX环境、真实构建、源码/PDF绑定和排版检查。
- `check_figure.py`、`figure_audit.py --questions q1 --strict`：导出及覆盖检查；Type0字体提示经所有图的pdffonts结果确认属于检查器误报。
- `python verify_translation.py`：八章结构、禁用措辞、数值、收敛及文件一致性。
- `python render_paper.py`：Poppler逐页渲染全部19页，并完成彩色/灰度图及页面视觉检查。

最终编译警告、未解析引用、未嵌入字体、空白页均为零。所有数值结论与本次结果一致；10月数据偏差、热条件与姿态假设已在正文中说明。
