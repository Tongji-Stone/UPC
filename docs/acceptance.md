# Acceptance record — coursework edition (2026-10-15)

The English paper retains its eight sections, model equations, numerical results, and single standalone-code appendix.

| Requested change | Result |
|---|---|
| Slightly larger body text and tighter spacing | PASS — 12 pt Times New Roman, 1.2 line spacing. |
| Contents matches body size | PASS — all contents entries and page numbers are 12 pt; the contents fit on one page. |
| New title | PASS — “An Aerodynamic and Thermodynamic Analysis of the Maximum Release Altitude for Space Diving.” |
| Coursework authors | PASS — 骆欣畅、魏贤东、刘奕韬 appear on the first page; romanized names appear in PDF metadata. No team-number placeholder remains. |
| Header layout and size | PASS — “Problem A: Space Diving” is on the left; “Page X of X” is on the right. |
| First-page summary | PASS — the problem, title, authors, and Summary appear. The PDF Summary contains 210 lexical words or 223 whitespace tokens, including numbers and units. |
| Formal references without URLs | PASS — 7 verified formal publications, each cited at least once in the main text; IEEE numbered format, complete author lists and report/standard identifiers. Webpage, web-archive, and product-sheet entries removed. Their data and assumptions remain in the text; the project source records retain their provenance without representing them as formal-publication evidence. |
| Larger and uniform figure text | PASS — legends, ticks, axis labels, annotations, and figure captions use an 11 pt base size. Mathematical subscripts and superscripts retain their normal smaller size. Plot PDF width matches its placement width, preventing unintended font shrinking. |
| References and equations | PASS — all numbered equations and all figures/tables have references; every bibliography entry is cited. All 39 displayed equations and all numerical results remain unchanged. |
| Independent code | PASS — complete solver remains unchanged and runs with site-packages disabled. Outputs remain 105.57896, 56.35661, and 66.28165 km. |
| Compilation and visual inspection | PASS — 32 pages, three XeLaTeX passes, no unresolved references, missing characters, or overfull boxes. All pages rendered; first-page Summary and page-two header inspected. |

`docs/checks/typography_revision.json` retains the earlier typography measurements and PDF hash; the latest reference checks are in `docs/checks/reference_audit_20261015.json`. `docs/checks/acceptance.json` and `main.build.json` record numerical, font, source, and PDF checks.

The prior model test suite (14 physical tests and 7 spin tests) is unchanged. Source hashes in `docs/numerical_baseline.json` confirm that this revision does not change the model, configuration, inputs, or independent core solver. The original repository metadata remains untouched.

Reference formatting was checked against the IEEE Reference Guide and IEEE author reference guidance. The current bibliography contains no online-only entries. Verification URLs are retained in the separate audit record, not the paper.

The revised first page and all 32 pages were rendered. The bibliography occupies its own page. Model code, inputs, and results are unchanged. This edition is for coursework rather than formal UPC submission.
