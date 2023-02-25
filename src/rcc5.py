#!/usr/bin/env python3
# ≥ ≤ ≳ ≲ ≂

# Actually this is RCC-29...

# intended for use with: from rcc5 import *

# RCC-5 relationship WITHIN hierarchy
# Use logical | for disjunction, & for conjunction

rcc5_symbols = {}
rcc5_eulerxs = {}
rcc5_relationships = {}
def def_rcc5_symbol(n, name, exname=None):
  if exname==None: exname = name
  rcc5_eulerxs[n] = exname
  rcc5_symbols[n] = name
  rcc5_relationships[name] = n
  rcc5_relationships[exname] = n
  return n

EQ = def_rcc5_symbol(1 << 0, '=')
LT = def_rcc5_symbol(1 << 1, '<')
GT = def_rcc5_symbol(1 << 2, '>')
DISJOINT = def_rcc5_symbol(1 << 3, '!')
CONFLICT = def_rcc5_symbol(1 << 4, '><')

LE = def_rcc5_symbol(LT|EQ, '<=', '{< =}')       # ≤, synonym
GE = def_rcc5_symbol(GT|EQ, '>=', '{> =}')       # ≥, accepted
NOINFO = def_rcc5_symbol(LT|GT|EQ|CONFLICT|DISJOINT, '?', '{< = > >< !}')
COMPARABLE = def_rcc5_symbol(LT|GT|EQ, 'comparable', '{< = >}')
OVERLAP = def_rcc5_symbol(LT|GT|EQ|CONFLICT, 'not!', '{< = > ><}')  # ∩, equivalent, similar

HAS_PARENT = LT
SYNONYM = LE
MYNONYS = GE
NEQ = def_rcc5_symbol(LT|GT|CONFLICT|DISJOINT, 'not=', '{< > >< !}')

PERI = def_rcc5_symbol(LT|DISJOINT, '<!', '{< !}')
IREP = def_rcc5_symbol(GT|DISJOINT, '>!', '{> !}')


# If we come across an int relation coding, there had better be such a relation

def rcc5_symbol(ship):
  return rcc5_symbols[ship]

def rcc5_eulerx(ship):
  return rcc5_eulerxs[ship]

def rcc5_relationship(name):
  return rcc5_relationships[name]

def reverse_relationship(ship):
  l = (ship & LT)
  g = (ship & GT)
  strip = (ship - (l|g))
  if l > 0: strip |= GT
  if g > 0: strip |= LT
  return strip
