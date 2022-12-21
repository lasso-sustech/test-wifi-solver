#!/usr/bin/env python3
import cvxpy as cp
import numpy as np
from utils import *


## Input
class NormalLink(LinkBase):
    LinkRate = 200*MB
    AC2 = [ ThruApp(min_thru=1*MB) ]

class NormalLinkAndDelay(LinkBase):
    LinkRate = 200*MB
    AC2 = [ ThruApp(min_thru=1*MB), DLApp(arrival=2.375*MB/PKT, max_qos=np.Inf) ]

class NormalLinkAndSidecar(LinkBase):
    LinkRate = 200*MB
    AC1 = [ RTApp(arrival=1.25*MB/PKT, max_qos=np.Inf) ]
    AC2 = [ ThruApp(min_thru=1*MB) ]

Links = [
    NormalLinkAndSidecar,
    NormalLinkAndDelay,
    NormalLink, NormalLink,
]


total_utility = 0
qos_function = 0
constraints = []
## Problems
for link in Links:
    link_utility = 0
    for aci, acq in enumerate(link.iter()):
        for app in acq:
            link_utility += app.utility / link.LinkRate
            qos_function += app.calc_qos()
            constraints.extend( app.constraints )
            pass
    total_utility += link_utility 
## bandwidth utility constraints
constraints.extend([
    total_utility >= 0.6,
    total_utility <= 1
])

## Solution
objective = cp.Minimize( qos_function )
problem = cp.Problem(objective, constraints)
problem.solve()

## Print
for link in Links:
    print(f'{link.name}: ')
    for aci, acq in enumerate(link.iter()):
        if acq:
            value = [x.variable.value/MB for x in acq]
            print(f'\tAC-{aci}: {value}')
