#!/usr/bin/env python3

import math
import property as prop
import checklist, workspace, simple
#import exemplar

from util import log, UnionFindable
from checklist import *
from workspace import *
from linkage import get_link
from linkage import pick_better_record
from typify import equate_typifications, find_typifications

# Exemplars (representing individual specimens or occurrences with
# known classification in BOTH checklists) are chosen heuristically
# via name matches.  Not every name match qualifies; only those that
# are 'tipward' do.  'Tipward' means that at least one of the two taxa
# has no tipward accepted descendant taxon.

# See AB.exemplar_ufs for map {xid: uf, ...}
# Each exemplar has a representative in A and one in B.
