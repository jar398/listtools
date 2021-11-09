import sys

def Property(key):
  return [key]

def set_property(x, prop, value):
  x[prop[0]] = value

nodefault = []

def getter(prop):
  key = prop[0]
  def get_property(x, default=nodefault):
    if default is nodefault:
      return x[key]
    else:
      return x.get(key, default)
  return get_property

def setter(prop):
  key = prop[0]
  def set_property(x, value):
    x[key] = value
  return set_property

IDENTITY_KEY = 0
global_counter = 0

def constructor(*props):
  keys = [prop[0] for prop in props]
  nkeys = len(keys)
  def cons(*values):
    global global_counter
    global_counter += 1
    if len(values) != nkeys:
      print("** Wrong number of arguments to constructor\n Have %s, want %s" %
            (len(values), keys),
            file=sys.stderr)
      assert False
    entity = {IDENTITY_KEY: global_counter}
    i = 0
    for key in keys:
      entity[key] = values[i]
      i += 1
    return entity
  return cons

_Identity = Property(IDENTITY_KEY)

get_identity = getter(_Identity)

def inspect(y):
  if isinstance(y, dict):
    for (key, val) in y.items():
      if key == IDENTITY_KEY:
        print("# Identity = %s" % (val,), file=sys.stderr)
      else:
        print("# %s = %s" % (key, inspect_scalar(val)), file=sys.stderr)
  else:
    print("# %s" % inspect_scalar(val), file=sys.stderr)

def inspect_scalar(val):
  if val == False:
    return "False"
  elif isinstance(val, str):
    return "'%s'" % val
  elif isinstance(val, int):
    return "%s" % val
  elif isinstance(val, dict) and IDENTITY_KEY in val:
    return "<%s>" % get_identity(val)
  else:
    return "..."

