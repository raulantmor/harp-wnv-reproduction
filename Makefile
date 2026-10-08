PYTHON ?= .venv/bin/python
.NOTPARALLEL:
export OPENBLAS_NUM_THREADS := 1
export OMP_NUM_THREADS := 1
export MPLCONFIGDIR := /tmp/harp-matplotlib

.PHONY: setup data prepare fit validate figures report all test
setup:
	bash scripts/setup.sh
data:
	$(PYTHON) src/download.py base
	$(PYTHON) src/download.py prism
	$(PYTHON) scripts/pdsi_parallel.py
prepare:
	$(PYTHON) src/prepare.py geography
	$(PYTHON) src/prepare.py climate
fit:
	$(PYTHON) src/forecast.py
	$(PYTHON) src/postcheck.py
validate:
	$(PYTHON) src/validate.py
test:
	$(PYTHON) -m pytest -q tests
figures:
	$(PYTHON) src/plot.py
report:
	bash scripts/check_tex.sh
	$(PYTHON) src/report.py
	cd report && latexmk -pdf -interaction=nonstopmode -halt-on-error report.tex
all: data prepare fit validate test figures report
