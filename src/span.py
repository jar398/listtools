#!/usr/bin/env python3

# Create a spanning tree for A+B based on B priority

# For each equivalent pair, only the B record will have a superior and inferiors.

# Would it be worthwhile to have a way to write a spanning tree as a
# checklist, so that it can be re-read without loss? ...
# Actually, would it be worthwhile to be able to do that with a workspace ...

import argparse, sys
import theory, rows
from checklist import rows_to_checklist, checklist_to_rows, get_superior, \
  relation, get_children, set_children, get_synonyms, set_synonyms, \
  get_taxonomic_status, set_taxonomic_status, \
  get_source, get_outject, blurb, is_toplike
import newick
import checklist, match_records, workspace
from checklist import monitor
from rcc5 import *
from util import log
from theory import block_size

# Assign a parent to every node of AB except for the top.

# B checklist has priority
# Each node z has two possible parents p and q, one in each checklist.
# (If z is equivalent to z' in the other checklist then by
# construction the possible parents of the two are the same.)
# The 'correct' parent r would be the least one that is greater than both:
# r = p if r < p < q, r = q if r < q < p.
# If p = q then pick the one in the priority checklist.
# If p and q are incomparable (>< or !) then the parent is the one
# from the same checklist (r = p if z in A, r = q if z in B).

def span(AB):

  AB.top = None                 # Updated further down

  def traverse(C):              # C is AB or swap(AB)

    for y in checklist.preorder_records(C.B):    # starts with top

      # y could be A.top or B.top
      w = C.in_right(y)
      assert get_superior(w, None) == None

      # A.top is AB.B.top and has no superior
      if y == AB.B.top:
        assert not AB.top
        AB.top = w
        continue

      sup = None

      # Normalize non-dominant w to dominant equivalent
      if theory.isinA(AB, w):
        v_rel = theory.get_equivalent(AB, w)
        # N.b. w and v could be top
        if v_rel:    # in A
          # If w has an equivalent in B, synonymize the two
          sup = relation(EQ, v_rel.record,
                         note="treat equivalent as synonym")
          set_taxonomic_status(w, "equivalent") # Synonym
      if not sup:
        # p is one possible parent of w
        p_rel = cross_superior(C, w)

        # q is another possible parent of w
        qy_rel = get_superior(y, None)

        # If there are two possibilities, pick the smallest
        if qy_rel and p_rel:
          p = p_rel.record
          q = C.in_right(qy_rel.record)
          ship = theory.compare(C, p, q).relationship
          blot = "%s %s %s" % (blurb(p), rcc5_symbol(ship), blurb(q))
          if ship == GT or ship == GE:
            sup = relation(qy_rel.relationship, q,
                           note="%s, choosing right" % blot,
                           span=qy_rel.span)
          elif ship == LT or ship == LE:
            sup = relation(p_rel.relationship, p,
                           note="%s, choosing left" % blot,
                           span=p_rel.span)
          elif theory.isinB(AB, w):
            sup = relation(qy_rel.relationship, q,
                           note="%s, choosing dominant (on right)" % blot,
                           span=qy_rel.span)
          elif ship == EQ:
            sup = relation(p_rel.relationship, p,
                           note="%s, choosing dominant (on left)" % blot)
          else:
            # Happens a lot:
            # log("incomparable parent candidates: %s" % blot)
            # q is parent of w
            if theory.isinA(AB, p):
              sup = cross_superior(AB, p)
            else:
              sup = cross_superior(AB, q)
        elif qy_rel:
          sup = relation(qy_rel.relationship,
                         C.in_right(qy_rel.record),
                         note="root 1")
        elif p_rel:
          sup = relation(p_rel.relationship,
                         p_rel.record,
                         note="root 2")
        else:
          log("no parent candidates!? %s" % blurb(w))

        # Normalize sup to dominant equivalent (transfer children and
        # synonyms over to dominant node)
        if sup != None and sup.record != None:
          r = sup.record
          assert r, "foooo"
          s_rel = theory.get_equivalent(AB, r)
          if s_rel and theory.isinB(AB, s_rel.record):
            sup = relation(sup.relationship, s_rel.record,
                           note=sup.note + " normalized")

      #log("sup %s = %s" % (blurb(w), blurb(sup)))

      if sup:                           # not top
        checklist.link_superior(w, sup) # adds w to children/synonyms list

  traverse(AB)
  traverse(theory.swap(AB))
  assert AB.top
  AB.indexed = True
  return AB


# -----------------------------------------------------------------------------
# Given a record v in the A checklist, return the least node in B
# that's definitely greater than v.

# AB.B could be priority, or not

def cross_superior(AB, v):
  # Smallest in opposite checklist that's >= v
  rel = theory.get_estimate(v)     # candidate in AB
  assert rel
  w = rel.record
  while True:
    # Looking for least w such that v <= w0 < w
    rel = theory.compare(AB, v, w)
    ship = rel.relationship
    if ship == LT or ship == LE:
      break
    log("# foo %s %s" % (blurb(v), blurb(rel)))
    wsup = theory.local_sup(AB, w)
    if not wsup:
      assert is_toplike(w)
      log("#  no node in B is bigger than %s" % blurb(v))
      return None
    w = wsup.record

  assert ship == LT or ship == LE, \
    (blurb(v), rcc5_symbol(ship), blurb(w))

  return relation(ship, w, note=rel.note)

# -----------------------------------------------------------------------------

def test():
  def testit(m, n):
    log("\n-- Test: %s + %s --" % (m, n))
    a_iter = newick.parse_newick(m)
    b_iter = newick.parse_newick(n)
    AB = ingest_workspace(a_iter, "A",
                          b_iter, "B",
                          None)
    theory.theorize(AB)
    S = span(AB)
    print(newick.compose_newick(checklist.preorder_rows(S)))
  testit("a", "a")              # A + B
  testit("(c,d)a", "(c,e)b")
  testit("((a,b)e,c)d", "(a,(b,c)f)D")

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""
    """)
  parser.add_argument('--A', help="the A checklist, as path name or -",
                      default='-')
  parser.add_argument('--B', help="the B checklist, as path name or -",
                      default='-')
  parser.add_argument('--test', action='store_true', help="run smoke tests")
  args=parser.parse_args()
  if args.test:
    test()
  else:
    a_path = args.A
    b_path = args.B
    assert a_path != b_path
    with rows.open(a_path) as a_rows:
      with rows.open(b_path) as b_rows:
        # compute name matches afresh
        AB = ingest_workspace(a_rows.rows(), b_rows.rows())
        find_links(AB, matches_iter)
        find_links(AB)
        report = span(AB)       # without alignment!!??  no way.
        util.write_rows(report, sys.stdout)
