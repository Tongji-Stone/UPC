# Space Diving — UPC 2023, Problem A

This is the complete English edition of the current manuscript. The model, parameters, input measurements, and numerical methods are unchanged by translation. All paper figures are regenerated in English. Figure 12 spans geometric altitudes from 0 to 20 km; actual free-fall trajectories end at parachute deployment.

This coursework edition shows the three student names on the title page and records their romanized names in the PDF metadata. Its first-page Summary is 223 whitespace-delimited tokens, including numbers and units. The model's 56.36 km result remains conditional on the declared equipment-screen assumptions.

## Read the paper

- `main.pdf`: compiled English paper.
- `main.tex`: authoritative editable source.
- `figures/`: 14 data figures, each in PDF, SVG, and PNG; three force/geometry diagrams are embedded TikZ in `main.tex`.
- `docs/terminology.md`: translation glossary.
- `docs/acceptance.md`: acceptance record.

## Standalone solver

```powershell
python -I -S appendices/solve.py
```

Python 3.10 or newer is recommended. This file uses only the Python standard library (`math` and `bisect`); it requires no other project files, packages, or external programs. It contains the adaptive Runge–Kutta solver, interpolation and ground-event handling, the three-flight area calibration, atmospheric and gravity models, free-fall equations, recovery temperature, analytical local flat-spin response, altitude-boundary searches, parachute inflation, and fully inflated descent to landing.

The only console output is three heights in km, in the order free-fall load, thermal condition, and flat spin:

```text
105.57896
56.35661
66.28165
```

The first value uses the explicitly isolated atmospheric extrapolation to 110 km. All other calculations retain the 86 km atmospheric domain. Parachute and landing results are computed internally and are available in the dictionary returned by `main()`; they do not add console output. The appendix includes this exact file, without transcription or omitted solver routines.

## Full research reproduction

The modular research implementation uses SciPy and data files to regenerate the complete results, including the observational comparisons and sensitivity results already present in the source project. Its printed altitude boundaries are 105.58, 56.36, and 66.28 km. The standalone implementation differs by less than 2 m because it uses its own integrator and calibration search; the paper's rounded values are identical.

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe reproduce.py --compile
```

Install local XeLaTeX, Poppler (`pdftoppm`), Times New Roman (regular, italic, bold, bold italic), SimSun, and Latin Modern Mono. The build runs XeLaTeX three times and does not require Perl or latexmk. The paper uses `newtxmath` for Times-style body mathematics. Data figures and TikZ labels use actual Times New Roman. The appendix retains Latin Modern Mono at the same 12 pt nominal size as the body.

For a quick rebuild from the supplied results:

```powershell
python make_figures.py
python generate_paper_data.py
python build_paper.py
python verify_project.py
python render_paper.py
```

On Windows, `build.ps1 -SkipCalculations` compiles and verifies the supplied current inputs. Without that switch it runs the full reproduction pipeline.

All paths resolve relative to this project. Offline reproduction uses the included flight summary, digitized speed curve and source image, and selected sounding records. Source URLs and provenance are in `data/source_registry.json` and `data/parameter_provenance.json`. No other project directory is needed.

## Layout and numerical checks

`test_model.py` and `test_spin_model.py` contain 21 physical/numerical checks. `verify_appendix.py` checks standalone-solver equivalence, parachute continuity, and landing against the modular solver. `verify_project.py` checks terminology, references, all figure fonts, Figure 12's data range, numerical invariance, isolated standard-library execution, appendix font size, and PDF/source bindings. `render_paper.py` renders every page and creates contact sheets in `build/revision/`.

Generated build products stay under `build/`; the published PDF is `main.pdf`. Run `package_project.py` to create `delivery/UPC_English.zip`. The archive excludes Git metadata, environments, caches, intermediate LaTeX products, and rendered page images. The original repository history is retained locally and is not altered by the replacement.
