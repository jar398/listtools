#!/usr/bin/env python3
# ≥ ≤ ≳ ≲ ≂

# Actually this is RCC-29...

# intended for use with: from rcc5 import *

# RCC-5 relationship WITHIN hierarchy
# Use logical | for disjunction, & for conjunction

rcc5_symbols = {}
rcc5_relationships = {}
def def_rcc5_symbol(n, name):
  rcc5_symbols[n] = name
  rcc5_relationships[name] = n
  return n

# want rel and ~ / rel but-not ~ , eg. ≳ ≲ ≳≲ 
# Maybe I just need a completely separate set of operators for
# the 'tipes.'

EQ = def_rcc5_symbol(1 << 0, '=')
LT = def_rcc5_symbol(1 << 1, '<')
GT = def_rcc5_symbol(1 << 2, '>')
DISJOINT = def_rcc5_symbol(1 << 3, '!')
CONFLICT = def_rcc5_symbol(1 << 4, '><')

LE = def_rcc5_symbol(LT|EQ, '{< =}')       # ≤, synonym
GE = def_rcc5_symbol(GT|EQ, '{> =}')       # ≥, accepted
NOINFO = def_rcc5_symbol(LT|GT|EQ|CONFLICT|DISJOINT, '{< = > >< !}')
COMPARABLE = def_rcc5_symbol(LT|GT|EQ, '{< = >}')
OVERLAP = def_rcc5_symbol(LT|GT|EQ|CONFLICT, '{< = > ><}')  # ∩, equivalent, similar

HAS_PARENT = LT
SYNONYM = LE
MYNONYS = GE
NEQ = def_rcc5_symbol(LT|GT|CONFLICT|DISJOINT, '{< > >< !}')

PERI = def_rcc5_symbol(LT|DISJOINT, '{< !}')
IREP = def_rcc5_symbol(GT|DISJOINT, '{> !}')


# If we come across an int relation coding, there had better be such a relation

def rcc5_symbol(x):
  return rcc5_symbols[x]

def rcc5_relationship(x):
  return rcc5_relationships[x]

def reverse_relationship(rel):
  l = (rel & LT)
  g = (rel & GT)
  strip = (rel - (l|g))
  if l > 0: strip |= GT
  if g > 0: strip |= LT
  return strip
