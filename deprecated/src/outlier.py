# -----------------------------------------------------------------------------
# Decide about relations between two trees
# This code is not used - it is an alternative to the set-based "blocks" further down

# Returns a pair (overlap, outlier)
#   overlap is in y and x {< = > ><} y (but not !)
#   outlier is in y and x {< >< !} y (but not {> =}) x ⋧ y
#   either can be None
# x is in the A checklist, y is in the B checklist.

def analyze(AB, x, y):
  assert False   # Not in use yet!
  (over, out) = (None, None)
  for d in get_inferiors(y):    # in B
    (over2, out2) = analyze(x, d)
    over = over or over2; out = out or out2
    if over and out: return (over, out)
  # TBD: if over and out are both present, and one of the two is a synonym ...
  if over or out:      # Exclude peripherals
    return (over, out)
  m = AB.record_match(y)
  j = AB.join(m, y)
  return (j, None) if m else (None, j)

def outlier(AB, x, y):
  (_, out) = analyze(x, y)
  return out

# RCC5 decision procedure, where x and y are in different sources.

def cross_ge(AB, x, y):
  return not outlier(AB, x, y)

def cross_eq(AB, x, y):
  # see if they're peers ??
  return cross_ge(x, y) and cross_ge(y, x)

def gt(AB, x, y):
  return cross_ge(x, y) and not cross_ge(y, x)

