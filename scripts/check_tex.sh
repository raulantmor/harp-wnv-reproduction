#!/usr/bin/env bash
set -euo pipefail
for package in article.cls hyperref.sty booktabs.sty natbib.sty graphicx.sty amsmath.sty xcolor.sty caption.sty subcaption.sty geometry.sty longtable.sty float.sty; do
    kpsewhich "$package"
done
