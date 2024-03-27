#!/bin/bash
set -e

# checklistbank demo

function doit {
  set -v
  echo --------- work/msw3.csv
  rm -f work/msw3-clean.csv work/msw3.csv 
  make work/msw3.csv
  echo --------- work/mdd1.10.csv
  rm -f work/mdd1.10-clean.csv work/mdd1.10.csv
  make work/mdd1.10.csv
  echo --------- work/col22_pteriomorphia.csv
  rm -f work/col22_pteriomorphia-clean.csv work/col22_pteriomorphia.csv
  make work/col22_pteriomorphia.csv
  echo --------- work/col23_pteriomorphia.csv
  rm -f work/col23_pteriomorphia-clean.csv work/col23_pteriomorphia.csv
  make work/col23_pteriomorphia.csv
  echo --------- msw3-mdd1.10-exemplars-v5.csv
  time src/exemplar.py --A work/msw3.csv --B work/mdd1.10.csv \
                > work/msw3-mdd1.10-exemplars-v5.csv
  echo --------- msw3-mdd1.10-plugin-v5.csv
  time src/plugin.py --A work/msw3.csv --B work/mdd1.10.csv \
                --exemplars work/msw3-mdd1.10-exemplars-v5.csv \
                > work/msw3-mdd1.10-plugin-v5.csv
  echo --------- col22_col23_pteriomorphia-exemplars-v5.csv
  time src/exemplar.py --A ~/g/listtools/work/col22_pteriomorphia.csv \
                --B ~/g/listtools/work/col23_pteriomorphia.csv \
                > work/col22_col23_pteriomorphia-exemplars-v5.csv
  echo --------- col22_col23_pteriomorphia-plugin-v5.csv
  time src/plugin.py --A work/col22_pteriomorphia.csv \
                --B work/col23_pteriomorphia.csv \
                --exemplars work/col22_col23_pteriomorphia-exemplars-v5.csv \
                > work/col22_col23_pteriomorphia-plugin-v5.csv
}

doit >work/v5.txt 2>&1
