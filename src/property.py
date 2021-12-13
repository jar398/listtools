import sys
from typing import NamedTuple, Any
from collections import namedtuple

# Instances, properties, and plans (field access for csv file rows
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
  return make_plan([get_property(label) for label in header])

def make_plan(props):
  if len(props) > 0:
    assert isinstance(props[0], Property), props
  propid_to_pos = [None] * global_property_counter
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

def get_property(label, filler=None, getter=None, setter=None):
  global _global_property_counter
  prop = _label_to_property.get(label)
  if not prop:
    prop = Property(_global_property_counter, label, filler, getter, setter)
    _global_property_counter += 1
    _label_to_property[label] = prop
  return prop

def get_set(prop):
  return (getter(prop), setter(prop))

defaulting = ['defaulting']

def getter(prop):
  if prop.getter: return prop.getter
  setit = setter(prop)
  def getit(x, default=_NODEFAULT):
    pos = x.plan.propid_to_pos[prop.id]
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
      else:
        assert False, "missing value for property %s" % prop.label
    elif stored == _SHADOW: val = MISSING
    elif stored == _CYCLING:
      assert False, "cycled while computing %s" % prop.label
    else: val = stored
    return val
  return getit

def setter(prop):
  if prop.setter: return prop.setter
  def setit(x, val):
    pos = x.plan.propid_to_pos[prop.id]
    stored = _SHADOW if val == MISSING else val
    if pos != None:
      x.positional[pos] = stored
    else:
      x.lookedup[prop.id] = stored
  return setit

# Instances

class Instance(NamedTuple):
  id: int
  plan: Any
  positional: Any
  lookedup: bool

global_instance_counter = 0

def constructor(*props, more=[]):
  all = props + more
  plan = make_plan(all)
  def constructit(*values):
    assert len(values) == len(props), \
      (len(values), [prop.label for prop in props])
    return construct(plan, values + ((MISSING,) * len(more)))
  return constructit

def construct(plan, row):
  global global_instance_counter
  if len(row) != len(plan.props):
    print("** WNA: have %s args, expect %s" %
          (len(row),
           [prop.label for prop in plan.props]),       
          file=sys.stderr)
    assert False
  instance = Instance(global_instance_counter,
                      plan,
                      row,
                      {})
  global_instance_counter += 1
  return instance

def generate_rows(instance_generator, props):
  yield [prop.label for prop in props]
  getters = map(getter, props)
  for x in instance_generator:
    yield [g(x, MISSING) for g in getters]

# Test row-based instances

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
  print(get_n(x))
