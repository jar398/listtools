#!/bin/bash
set -e
set -v

function doit {
make work/msw3.csv
make work/mdd1.10.csv
make work/col22_pteriomorphia.csv
make work/col23_pteriomorphia.csv
time src/exemplar.py --A work/msw3.csv --B work/mdd1.10.csv \
              > msw3-mdd1.10-exemplars-v5.csv
time src/plugin.py --A work/msw3.csv --B work/mdd1.10.csv \
              --exemplars msw3-mdd1.10-exemplars-v5.csv \
              > msw3-mdd1.10-plugin-v5.csv
time src/exemplar.py --A ~/g/listtools/work/col22_pteriomorphia.csv \
              --B ~/g/listtools/work/col23_pteriomorphia.csv \
              > col22_col23_pteriomorphia-exemplars-v5.csv
time src/plugin.py --A work/col22_pteriomorphia.csv \
              --B work/col23_pteriomorphia.csv \
              --exemplars col22_col23_pteriomorphia-exemplars-v5.csv \
              > col22_col23_pteriomorphia-plugin-v5.csv
}

doit 2>v5.txt
