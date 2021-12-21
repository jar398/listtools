#!/usr/bin/env python3

import sys
from typing import NamedTuple, Any
from collections import namedtuple
from util import log

# Records, properties, and plans (field access for csv file rows
# rows)

MISSING = ''   # What is stored for missing values
_SHADOW = ["<SHADOW>"]   # What is stored when value should be MISSING
_NODEFAULT = ["<NODEFAULT>"]   # Used to detect that default is unsupplied
_CYCLING = ["<CYCLING>"]

# ----- Plans

class Plan(NamedTuple):
  propid_to_pos: Any
  props: Any

def make_plan_from_header(header):
  #log("%s" % (header,))
  props = [get_property(label) for label in header]
  #log("%s" % (tuple(((prop.id, prop.label) for prop in props)),))
  return make_plan(props)

def make_plan(props):
  if len(props) > 0:
    assert isinstance(props[0], Property), props
  propid_to_pos = [None] * _global_property_counter
  for (prop, pos) in zip(props, range(len(props))):
    propid_to_pos[prop.id] = pos
  return Plan(propid_to_pos, props)

# ----- Properties

# How about:
#   - properties that are just functions - oh yeah, they're just
#     functions
#   - properties that cache fixed functions - oh hmm

class Property(NamedTuple):
  id: Any
  label: str
  filler: Any                   # cached
  getter: Any
  setter: Any

_global_property_counter = 0
_label_to_property = {}

def get_property(label, filler=None, getter=None, setter=None, fresh=False):
  global _global_property_counter
  prop = _label_to_property.get(label)
  if fresh:
    if prop:
      log("# Creating another property with label %s" % label)
    prop = None
  if not prop:
    prop = Property(_global_property_counter, label, filler, getter, setter)
    _global_property_counter += 1
    _label_to_property[label] = prop
  return prop

# You can explicitly pass context=None to get ambient context.

def get_set(prop, context=None):
  return (getter(prop, context=context),
          setter(prop, context=context))

def getter(prop, context=None):
  if context:
    return contextual_getter(prop, context)
  elif prop.getter:
    return prop.getter
  else:
    return ambient_getter(prop)

def ambient_getter(prop):
  setit = setter(prop)
  def getit(x, default=_NODEFAULT):
    if (not x.cloned_from and
        prop.id < len(x.plan.propid_to_pos)):
      pos = x.plan.propid_to_pos[prop.id]
    else:
      pos = None
    if pos != None:
      stored = x.positional[pos]
    else:
      stored = x.lookedup.get(prop.id, MISSING)
    if stored == MISSING:
      if default != _NODEFAULT:
        val = default
      elif prop.filler:
        setit(x, _CYCLING)
        val = prop.filler(x)
        setit(x, val)
      elif x.cloned_from:
        return getit(x.cloned_from)
      else:
        assert False, \
          ("missing value for property '%s' position %s plan %s keys %s" %
           (prop.label, pos, x.plan, tuple(x.lookedup.keys())))
    elif stored == _SHADOW: val = MISSING
    elif stored == _CYCLING:
      assert False, "cycled while computing %s" % prop.label
    else: val = stored
    return val
  return getit

def setter(prop, context=None):
  if context:
    return contextual_setter(prop, context)
  else:
    return ambient_setter(prop)

def ambient_setter(prop):
  if prop.setter: return prop.setter
  def setit(x, val):
    stored = _SHADOW if val == MISSING else val
    if (not x.cloned_from and
        prop.id < len(x.plan.propid_to_pos)):
      pos = x.plan.propid_to_pos[prop.id]
    else:
      pos = None
    if pos != None:
      x.positional[pos] = stored
    else:
      x.lookedup[prop.id] = stored
  return setit

# e.g. checklist.get_parent = contextual_getter(get_property("parent"), checklist)

def contextual_getter(prop, context):
  column_mep = _get_column_mep(prop, context)
  def getit(inst, default=_NODEFAULT):
    stored = mep_get(column_mep, inst, MISSING)
    if stored == MISSING:
      if default != _NODEFAULT:
        val = default
      elif prop.filler:
        assert not "TBD"
      else:
        val = getter(prop)(inst)      # Inherit from ambient
    else:
      # what if there's a filler?  should use it?  how does it know
      # the context?
      val = stored
    return val
  return getit

def contextual_setter(prop, context):
  column_mep = _get_column_mep(prop, context)
  def setit(inst, val):
    column_mep[inst.id] = val
  return setit

# get_column_mep maps from property to column_mep (column_mep is a 'mep')
#   within a given context (something with a .columns property)

def _get_column_mep(prop, context):
  return get_column(prop, context).record_to_value

# There is one Column per property per context...

def get_column(prop, context):
  assert isinstance(context, Context), context
  column = mep_get(context.columns, prop, None)
  if column == None:
    column = Column(mep(), {})
    mep_set(context.columns, prop, column)
  return column

class Column(NamedTuple):
  record_to_value : Any
  value_to_record : Any
  # maybe more later

# ----- Completely generic csv loader for any kind of table.

def get_records(column):
  assert isinstance(column, Column)
  return column.value_to_record.values()

def get_record(column, value, default=None):
  return column.value_to_record[value]

def table_to_context(row_iterable, primary_key_prop):
  get_pk = getter(primary_key_prop) # ambient
  Q = make_context()
  pk_col = get_column(primary_key_prop, Q)
  assert isinstance(pk_col, Column)
  inverse_pk_dict = pk_col.value_to_record
  row_iterator = iter(row_iterable)
  header = next(row_iterator)
  plan = make_plan_from_header(header)
  for row in row_iterator:
    inst = construct(plan, row)
    inverse_pk_dict[get_pk(inst)] = inst  # ambient -> contextual
  return Q

def get_registrar(pk_prop, Q):
  get_pk = getter(pk_prop)      # ambient
  pk_col = get_column(pk_prop, Q)
  inverse_pk_dict = pk_col.value_to_record
  def register(inst):
    inverse_pk_dict[get_pk(inst)] = inst  # ambient -> contextual
  return register

class Context(NamedTuple):
  columns : Any  # mep()    # property -> (usage to value, value to usage)

def make_context():
  return Context(mep())

# Records

class Record(NamedTuple):
  id: int
  plan: Any
  positional: Any
  lookedup: bool
  cloned_from: Any

def constructor(*props, more=()):
  all = props + more
  plan = make_plan(all)
  def constructit(*values):
    assert len(values) == len(props), \
      (len(values), [prop.label for prop in props])
    return construct(plan, values + ((MISSING,) * len(more)))
  return constructit

def construct(plan, row):
  global _global_record_counter
  if len(row) != len(plan.props):
    print("** WNA: have %s args, expect %s" %
          (len(row),
           [prop.label for prop in plan.props]),       
          file=sys.stderr)
    assert False
  record = Record(_global_record_counter,
                      plan,
                      row,
                      {},
                      None)
  _global_record_counter += 1
  return record

def clone(inst):
  global _global_record_counter
  cloned = Record(_global_record_counter,
                    empty_plan, [], {}, inst)
  _global_record_counter += 1
  return cloned
  
empty_plan = make_plan(())

_global_record_counter = 0


# Maps keyed by record

_nodefault = []
def mep(): return {}
def mep_get(mep, inst, default=_nodefault):
  if default is _nodefault:
    return mep[inst.id]
  else:
    return mep.get(inst.id, default)
def mep_set(mep, inst, j):
  mep[inst.id] = j


# Test records

if __name__ == '__main__':
  a = get_property('a')
  b = get_property('b')
  c = get_property('c')
  m = get_property('m')
  n = get_property('n', filler=lambda x:("filled", get_a(x)))
  plan = make_plan_from_header(['a', 'b', 'c'])
  x = construct(plan, ['a value', 'b value', 'c value'])
  (get_a, set_a) = get_set(a)
  get_b = getter(b)
  print(get_a(x))
  print(get_b(x))
  (get_c, set_c) = get_set(c)
  print(get_c(x, 'no c'))
  set_a(x, 'new a value')
  print(get_a(x))
  set_c(x, 'c value')
  print(get_c(x))
  q = constructor(a, b, more=(m,))
  y = q('a it', 'b it')
  print(get_b(x, 'b is'))
  print(get_c(x, 'no c'))
  (get_m, set_m) = get_set(m)
  print(get_m(x, 'no m'))
  set_m(x, 'm value')
  print(get_m(x))
  (get_n, set_n) = get_set(n)
  print(get_n(x))
  set_a(x, "changed a")
  print(get_a(x))

  q = make_context()
  (q_get_a, q_set_a) = get_set(a, context=q)
  q_set_a(x, 'a in q')
  print(q_get_a(x))
  print(get_a(x))
