
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
    assert len(values) == nkeys
    entity = {IDENTITY_KEY: global_counter}
    i = 0
    for key in keys:
      entity[key] = values[i]
      i += 1
    return entity
  return cons

_Identity = Property(IDENTITY_KEY)

get_identity = getter(_Identity)
