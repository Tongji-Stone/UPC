# Space Diving（UPC 2023 Problem A）论文草稿

本目录以随附的同济 Team 355 论文为版式参考，提供一份重新建模、重新写作的英文 LaTeX 稿。请先阅读论文中“条件性高度筛选”的定义：**57.7 km 不是已证实的人体安全高度**。题目给定的 190 kg 不能单独决定最高可跳高度。

## 文件

- main.tex：可编辑论文；main.pdf：9 页 A4 编译稿。
- model.py：标准大气、自由落体、开伞情景、标定、数值扫描与制图。
- data/stratos_summary.csv：Stratos 科学报告公布的三次跳跃**汇总观测值**，并非原始遥测。
- results.json：运行模型生成的计算结果。
- figures 文件夹：论文图。

## 复现

在本目录运行：

    python model.py
    pdflatex -interaction=nonstopmode -halt-on-error main.tex
    pdflatex -interaction=nonstopmode -halt-on-error main.tex

Python 需要 numpy、scipy 和 matplotlib。MiKTeX 的 latexmk 在当前电脑缺少 Perl，因此使用两遍 pdflatex；已经实际编译成功。

## 数据和验证边界

1. 三次 Stratos 跳跃的高度、峰速、低于 0.1 g 的时长与开伞冲击来自 [Stratos Scientific Summit Report](https://lru.praxis.dk/Lru/microsites/hvadermatematik/hem3download/kap3a_QR20_ekstra_Report_Final.pdf)，主要是报告印刷页 7、10–11、46。3 月与 7 月峰速用于拟合；10 月事件留出检验。报告并未公开本稿所需的逐秒 CSV。
2. [Guerster 与 Walter 的 PLOS ONE 论文](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0187798)给出 Felix 连装备质量 121.2 kg 和姿态/阻力分析，但与任务报告同源，**不能算独立验证集**。
3. [FAI 的 Alan Eustace 记录](https://www.fai.org/news/yet-be-beaten-alan-eustaces-high-altitude-parachute-jump-records-still-stand-10-years)来自另一场飞行，质量接近 190 kg；该跳从头使用稳定伞，所以论文另拟合它的有效阻力，仅检查时间量，不把它当成裸跳模型的盲测。
4. 大气用 [1976 U.S. Standard Atmosphere](https://ntrs.nasa.gov/citations/19770009539)，这是参考剖面，不是某天的探空气球观测。热指标为绝热壁恢复温度近似值，不是衣料、皮肤或人体实测温度。
5. 个人主伞开启高度采用 2566.8 m MSL（约当地地面以上 5000 ft）。Stratos 报告叙述中的 5000 ft MSL 与任务记录存在基准面混用的迹象；本稿采用可与记录高度差相符的 MSL 值。载人舱回收伞不属于跳伞者个人主伞。

## 提交前

替换首页和页眉的队号占位符，并按当届 [UPC 规则](https://uphysicsc.com/contestrules.html)检查队伍身份、摘要字数、引用及 AI 使用披露。若需要把条件性筛选升格为装备安全结论，应补充宇航服热试验、开伞载荷测试、稳定控制参数和逐秒飞行测量。当前计算没有这些材料。
