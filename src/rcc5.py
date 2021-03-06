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

LE = def_rcc5_symbol(LT|EQ, '<=')       # ≤, synonym
GE = def_rcc5_symbol(GT|EQ, '>=')       # ≥, accepted
NOINFO = def_rcc5_symbol(LT|GT|EQ|CONFLICT|DISJOINT, '?')
COMPARABLE = def_rcc5_symbol(LT|GT|EQ, '<=>')
OVERLAP = def_rcc5_symbol(LT|GT|EQ|CONFLICT, '~')  # equivalent, similar

ACCEPTED = LT
SYNONYM = LE

"""

#     ≴ ≵  ≶  ≱  ≰  ≱  ℮
#  ≲ = Less-Than or Equivalent To
#  ⪞ = Similar or Greater Than
#  ≵ = Neither Greater-Than nor Equivalent To

ESTIMATES_EQ = def_rcc5_symbol(1 << (5+0), '℮=')       # = or similar i.e. =
ESTIMATES_LT = def_rcc5_symbol(1 << (5+1), '℮<')       # < or similar
ESTIMATES_GT = def_rcc5_symbol(1 << (5+2), '℮>')       # > or similar
ESTIMATES_DISJOINT = def_rcc5_symbol(1 << (5+3), '℮!')    # same as !
ESTIMATES_CONFLICT = def_rcc5_symbol(1 << (5+4), '℮><')   # same as ≂ in practice

ESTIMATES_LE = def_rcc5_symbol(ESTIMATES_LT|ESTIMATES_EQ, '℮≲')   # via mrca
ESTIMATES_GE = def_rcc5_symbol(ESTIMATES_GT|ESTIMATES_EQ, '℮≥')   # via mrca

ESTIMATES_NLE = def_rcc5_symbol(ESTIMATES_GT|ESTIMATES_DISJOINT|ESTIMATES_CONFLICT, '℮≰')   # via mrca
ESTIMATES_NGE = def_rcc5_symbol(ESTIMATES_LT|ESTIMATES_DISJOINT|ESTIMATES_CONFLICT, '℮≱')   # via mrca
"""



# If we come across an int relation coding, there had better be such a relation

def rcc5_symbol(x):
  return rcc5_symbols[x]

def rcc5_relationship(x):
  return rcc5_relationships[x]

def reverse_relationship(rel):
  l = (rel & LT)
  g = (rel & GT)
  strip = (rel - (l|g))
  if l > 0: strip |= g
  if g > 0: strip |= l
  return strip
