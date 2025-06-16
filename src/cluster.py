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

from rcc5 import DISJOINT
from checklist import get_inferiors, get_rank
from homotypy import compare_parts

# A cluster has an id and a list of records (in checklist)

# Clusters
records_prop = prop.declare_property("records")
cluster_id_prop = prop.declare_property("cluster_id")
make_cluster = prop.constructor(records_prop, cluster_id_prop)

get_records = prop.getter(records_prop)
get_cluster_id = prop.getter(cluster_id_prop)

cluster_id = 1

def match_clusters(A, B):
  A_clusters = find_clusters(A)
  B_clusters = find_clusters(B)
  by_epithet = index_clusters_by_epithet(A_clusters, B_clusters) 
  # by_epithet: ep -> [[acluster...], [bcluster...]]
  choices = []
  for sub in by_epithet.values():
    # sub = subproblem = (i_s, js)
    matches = match_in_subproblem(sub)
    # matches is a list of (i, j, score)
    # We can establish an exemplar for each choice.
    # exemplar = matched pair of clusters
    for ch in check_cluster_matches(matches):
      choices.append(ch)
  return create_exemplars(choices)

# Find clusters of same-type (homotypic) records in a checklist, returned
# as a list of clusters

# Side effect: 'cluster' field of nodes in C get set.

def find_clusters(C):
  clusters = []
  def process(x, cluster):
    if not cluster and get_rank(x) == "species":
      cluster = new_cluster(x)
      clusters.append(cluster)

    if cluster:
      l = get_leader(cluster)   # most rootward
      if (simple.compare_per_checklist(x, l) == DISJOINT or
          compare_parts(x, l) < MOTION):
        cluster = new_cluster(x)
        clusters.append(cluster)
      else:
        get_records(cluster).append(x)

    for c in get_inferiors(x):
      process(c, cluster)
  process(C.top, None)
  return clusters

def new_cluster(leader):
  global cluster_id
  cluster = make_cluster([leader], id)
  cluster_id += 1
  return cluster

# The record that "leads" a cluster (most rootward)

def get_leader(cluster):
  return get_records(cluster)[0]

# -----------------

# Returns: for each epithet, list of A-clusters and list of B-clusters.
# Could in principle examine all records for a variety of epithets?

def index_clusters_by_epithet(clusters1, clusters2):
  ix = {}                       # epithet -> (clusters, clusters)
  for i in clusters1:
    r1 = get_records(get_leader(i))[0]
    ep = get_parts(r1).epithet
    pair = ix.get(ep)           # (i_s, js)
    if not pair:
      pair = ([], [])
      ix[ep] = pair
    pair[0].append(i)
  for j in clusters2:
    s1 = get_records(get_leader(j))[0]
    ep = get_parts(s1).epithet
    pair = ix.get(ep)           # (i_s, js)
    if not pair:
      pair = ([], [])
      ix[ep] = pair
    pair[1].append(j)
  return ix

# Evaluate all candidate cluster matches in a subproblem.

def match_in_subproblem(sub):
  (i_s, js) = sub
  matches = []
  for i in i_s:
    for j in js:
      m = match_clusters(i, j)
      if m >= MOTION:
        matches.append((i, j, m))
  return sorted(matches,
                key=lambda quux: -quux[2])

def match_clusters(i, j):
  mmax = -1
  for r in get_records(i):
    for s in get_records(j):
      m = compare_parts(r, s)
      if m > mmax: mmax = m
  return mmax

# Check for ambiguities within a single subproblem
# input is a list (i, j, score) sorted by score
# output is a list (i, j, score) of matches

def check_cluster_matches(matches):
  ii = {}                       # (i, j, score) keyed by cluster id
  jj = {}
  for (i, j, m) in matches:
    # i, j are clusters
    probe = ii.get(get_cluster_id(i), None)
    if probe:
      (_, _, m1) = probe
      assert m >= m1
      if m == m1:
        assert False, "ambiguous"
    else:
      ii[get_cluster_id(i)] = (i, j, m)

    probe = jj.get(get_cluster_id(j), None)
    if probe:
      (_, _, m2) = probe
      assert m >= m2
      if m == m2:
        assert False, "ambiguous"
    else:
      jj[get_cluster_id(j)] = (i, j, m)

  assert len(ii) == len(jj)

  for (_, j1, m1) in ii.values():
    jid = get_cluster_id(j1)
    (i2, j2, m2) = jj[jid]
    assert i is i2 and j2 is j1, "non reciprocal"

  return list(ii.values())

# matches is a list of (i, j, score) where i,j are clusters

def create_exemplars(matches):
  id = 1
  es = []
  for (i, j, m) in matches:
    e = make_exemplar(i, j, m, id)
    id += 1
    for r in get_records(i):
      set_exemplar(r, e)
    for r in get_records(j):
      set_exemplar(r, e)
    es.append(e)
  return es
