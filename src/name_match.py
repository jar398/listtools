# Calculate match score / probability
# See wikipedia record linkage article

import util
import parse

def index_by_property(A, fn, value_set):
  index = {}
  for x in postorder_records(A):
    val = fn(x)
    if val:
      value_set.add(val)
      have = index.get(val)
      key = get_primary_key(x)
      if have:
        have.add(key)           # could cache protonymic for speed
      else:
        index[val] = {key}
  return index

def find_matches(AB):
  epithets = {}
  (A_index, B_index) = \
    = map(lambda Cop: \
          index_by_property(Cop.A,
                            lambda x: get_protonymic(x).epithet,
                            epithets),
          (AB, AB.swap()))   
  yield "A id", "B id"
  for ep in epithets:
    A_keys = A_index.get(ep, ())
    B_keys = B_index.get(ep, ())
    assert A_keys or B_keys
    n = max(len(A_keys), len(B_keys))
    cache = {}
    for Cop in (AB, AB.swap()):
      for x_key in A_keys:         # set of Record... ids?? for this epithet
        x = look_up_record(Cop.A, x_key)
        for y_key in B_keys:
          key = (x_key, y_key)            # ?
          probe = cache.get(key)
          if probe != None:
            prob = probe
          else:
            prob = prob_score(get_protonymic(...))
            cache[key] = prob
            cache[(y_key, x_key)] = prob
      for x_key in A_keys:
        x = look_up_record(Cop.A, x_key)
        candidates = [look_up_record(Cop.B, y_key)
                      for y_key in B_keys.get(ep, ())
                      if cache[(x_key, y_key)] >= threshold]
        # Now talk about candidates that don't have exactly one match!
        if len(candidates) == 1:
          set_match(x, rel(EQ, candidates[0], 'unique good match'))
        if len(candidates) == 0:
          set_match(x, rel(NOINFO, None, 'no match'))
        else:
          set_match(x_key, rel(NOINFO, None '|'.join(y_key, '')))  #.. ?????

# If n is the number of taxa with a given epithet, then:
# consider a possibility of x (external) that the true match is to
# something without an epithet.  So we need to expand from n to xn,
# say x = 1.1.  Then: there are (xn)^2 pairs,
# which by symmetry would reduce to
#   h = P(H) = xn*(xn+1)/2

# m = P(R ~ S given match) = P(E|H)
# u = P(R ~ S given no match) = P(E|¬H)

def score_protos(p, q):
  if proto_lt(q, p): (p, q) = (q, p)

  # (check to see if memoized)

  mu = (1, 1)               # TP or FN, FP or TN
  def contribute(m1, u1):
    mu[0] *= m1
    mu[1] *= u1
  if p.epithet and q.epithet:
    if p.epithet == q.epithet:
      contribute(1 - 1e-10, 1e-11) # TP, FN
    else:
      contribute(1e-9, 1 - 1e-8) # FP, TN
  if p.genus and q.genus:
    if p.genus == q.genus:
      contribute(1 - 1e-8, 1e-9)
    else:
      pass
  if p.token and q.token:
    if p.token == q.token:
      contribute(1 - 1e-6, 1e-7)
    else:
      pass
  if p.year and q.year:
    if p.year == q.year:
      contribute(1 - 1e-4, 1e-5)
    else:
      pass
  if (p.currentGenus and q.currentGenus and
      p.genus == None and q.genus == None):
    if p.currentGenus == q.currentGenus:
      contribute(1 - 1e-2, 1e-3)
  # more...

  return (m[0], u[0])

# Need a total ordering ideally... review this

def proto_lt(q, p): return q < p # Hmmmm

def prob_from_score(score, h):

  # Calculated evidence prior (?):
  # e = P(E) = P(E|H)P(H) + P(E|¬H)P(¬H)
  e = m*h + u*(1-h)

  # Score = P(H|E) = P(H) * P(E|H) / P(E)
  return h*m/e


def get_protonymic(x):
  verbatim = get_scientific(x) or get_canonical(x)
  nymic = parse.parse_protonymic(all,
                                 verbatim,
                                 gn_full=get_gn_full(x, None),
                                 gn_authorship=get_gn_authorship(x, None),
                                 gn_stemmed=get_stemmed(x, None))
  
