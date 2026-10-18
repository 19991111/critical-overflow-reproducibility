# Build the current paper

Install requirements-paper.txt in the same Python virtual environment as the core. System prerequisites are XeLaTeX, BibTeX, standard SIAM class dependencies including TeX Live algorithms / algorithm.sty, and Poppler pdfinfo, pdffonts, pdftoppm and pdftotext. Typical Debian packages are texlive-xetex, texlive-latex-extra, texlive-science and poppler-utils.

```sh
.venv/bin/python -m pip install -r requirements-paper.txt
.venv/bin/python -B tools/build_paper.py --out runs/paper
```

Use a new output directory. This verifies the original ZIP and current editable TeX, regenerates the six figures and their raw-data inputs, verifies all four tables, and starts with no retained PDF/auxiliary output. The two documents compile in alternating rounds. BibTeX runs from the output directory, respecting openout_any=p. Three rounds stabilized the tested build. The SIAM class, TeX and bibliography text are unchanged.

The original BUILD.txt remains in the authoritative ZIP. It used Tectonic 0.17.0. This environment tested native XeLaTeX from TeX Live 2025/dev, the same engine family. The original compatibility code fails under the tested pdfLaTeX despite claiming no effect there; the failure log is retained. No security policy or document layout was changed to pass the build.

The test machine lacked algorithm.sty. An exact, checksum-verified unmodified v0.1 package was supplied externally through normal TEXINPUTS; provenance is in evidence/TeX_RUNTIME_DEPENDENCY.json. A standard system installation of TeX Live algorithms supplies this dependency without that test-machine path. The repository requires no original research-workspace absolute path and selects no new license for this runtime dependency or the paper.

PDFs appear under runs/paper/paper/output/pdf/, regenerated PDF/SVG/PNG figures under runs/paper/paper/manuscript/figures/, and every page render under runs/paper/render/. Detailed numeric, input and document checks are in the new output root. Checked current PDFs are also included in paper/output/pdf/ for direct reading.

Changing authoritative TeX or figure source is a new version, not a silent replay of the frozen input. Preserve the original ZIP, document edits and update source mapping/checks explicitly before claiming a newly validated paper version.
