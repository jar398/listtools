#!/usr/bin/env python3

# Match specimen references
# Either in same checklist, or in different checklists
# Replaces typify.py

# By 'specimen' we're usually talking the type specimen of some taxon,
# but it would be nice eventually to generalize to broader classes of
# specimen references

# Probably the best way to do this would be to create an additional
# artificial taxon record for each non-taxon (non-type) specimen.
# The specimen would look like the type of the artificial taxon, and 
# all the rest of the machinery would just work.

import util
import property as prop
import parse
import simple

from util import log
from rcc5 import DISJOINT
from checklist import get_inferiors, get_rank, get_parts, blurb, blorb
from homotypy import compare_parts, MOTION

# An exemplar has fields (A_cluster, B_cluster, match_status, id)

left_cluster_prop = prop.declare_property("left_cluster")
right_cluster_prop = prop.declare_property("right_cluster")
match_status_prop = prop.declare_property("match_status")
exemplar_id_prop = prop.declare_property("exemplar_id")

left_cluster = prop.getter(left_cluster_prop)
right_cluster = prop.getter(right_cluster_prop)
match_status = prop.getter(match_status_prop)
exemplar_id = prop.getter(exemplar_id_prop)

make_exemplar = prop.constructor(left_cluster_prop,
                                 right_cluster_prop,
                                 match_status_prop,
                                 exemplar_id_prop)

# A cluster is a list of records (drawn from a single checklist)

# Returns a list of exemplars
# An exemplar is a pair of clusters, related by name match level.

def find_exemplars(A, B):
  A_clusters = find_clusters(A)
  B_clusters = find_clusters(B)
  by_epithet = index_clusters_by_epithet(A_clusters, B_clusters) 
  # by_epithet: ep -> [[acluster...], [bcluster...]]
  exemplars = []
  id = 1
  for sub in by_epithet.values():
    # sub = subproblem = [[acluster...], [bcluster...]]
    matches = match_in_subproblem(sub)
    # matches is a list of (i, j, score)
    # We can establish an exemplar for each match.
    # exemplar = matched pair of cluster lists
    for match in matches:
      (i, j, m) = match
      exemplars.append(make_exemplar(i, j, m, id))
      id += 1
  return exemplars

# Find clusters of same-holotype (homotypic) records in a checklist
# Returns: a list of clusters [x]

# Side effect: 'cluster' field of nodes in C get set.

def find_clusters(C):
  clusters = []
  def process(x, cluster):

    def starts_new_cluster(x, cluster):
      if cluster:
        l = get_leader(cluster)   # most rootward
        if simple.compare_per_checklist(x, l) == DISJOINT:
          return True
        elif compare_parts(x, l) < MOTION:
          return True
      else:
        if get_rank(x) == "species":
          return True
      return False

    if starts_new_cluster(x, cluster):
      clusters.append([x])
    elif cluster:
      cluster.append(x)

    for c in get_inferiors(x):
      process(c, cluster)
  process(C.top, None)
  return clusters

# The record that "leads" a cluster (most rootward)

def get_leader(cluster):
  return cluster[0]

# -----------------

# Returns: for each epithet, list of A-clusters and list of B-clusters.
# Could in principle examine all records for a variety of epithets?

def index_clusters_by_epithet(clusters1, clusters2):
  ix = {}                       # epithet -> (clusters, clusters)
  for i in clusters1:
    r1 = get_leader(i)
    ep = get_parts(r1).epithet
    pair = ix.get(ep)           # (i_s, js)
    if not pair:
      pair = ([], [])
      ix[ep] = pair
    pair[0].append(i)
  for j in clusters2:
    s1 = get_leader(j)
    ep = get_parts(s1).epithet
    pair = ix.get(ep)           # (i_s, js)
    if not pair:
      pair = ([], [])
      ix[ep] = pair
    pair[1].append(j)
  return ix

# Evaluate all candidate cluster matches in a subproblem.

def match_in_subproblem(sub):
  (i_s, js) = sub               # cluster lists
  cluster_matches = []
  for i in i_s:
    for j in js:
      m = match_clusters(i, j)  # match status
      if m >= MOTION:
        cluster_matches.append((i, j, m))
  # cluster_matches is a list of (A_cluster, B_cluster, score)
  cluster_matches.sort(key=lambda quux: -quux[2])
  # Remove redundant matches
  thinned = []
  A_matches_mep = prop.mep()
  B_matches_mep = prop.mep()
  # Sorted so that large m's are processed first
  for cm in cluster_matches:
    (i, j, m) = cm

    cm_before = prop.mep_get(A_matches_mep, get_leader(i), None)
    if cm_before:
      (_, j2, m2) = cm_before
      # ambiguity if match scores are the same
      if m == m2:
        log("# ambiguous %s -> %s, %s (%s)" %
            (blorb(get_leader(i)), blorb(get_leader(j2)),
             blorb(get_leader(j)), m))
      continue
    else: prop.mep_set(A_matches_mep, get_leader(i), cm)

    cm_before = prop.mep_get(B_matches_mep, get_leader(j), None)
    if cm_before:
      (i2, _, m2) = cm_before
      if m == m2:
        log("# ambiguous %s, %s -> %s (%s)" %
            (blorb(get_leader(i)), blorb(get_leader(i2)),
             blorb(get_leader(j)), m))
      continue
    else: prop.mep_set(B_matches_mep, get_leader(j), cm)

    thinned.append(cm)
  return thinned

# How well does a given cluster match one in the opposite checklist?

def match_clusters(i, j):
  mmax = -1
  for r in i:
    for s in j:
      m = compare_parts(r, s)
      if m > mmax:
        mmax = m
  return mmax
