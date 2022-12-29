#!/usr/bin/env python3
import cvxpy as cp
import numpy as np
from utils import *


## Input
class NormalLink1(LinkBase):
    LinkRate = 200*MB
    AC2 = [ ThruApp(min_thru=1*MB, weight=1) ]

class NormalLink2(LinkBase):
    LinkRate = 200*MB
    AC2 = [ ThruApp(min_thru=1*MB, weight=1) ]

class NormalLinkAndDelay(LinkBase):
    LinkRate = 200*MB
    AC2 = [ ThruApp(min_thru=1*MB, weight=1), 
            DLApp(pkt_size=PKT, arrival=2.375*MB, max_qos=2000, weight=5000) ]

class NormalLinkAndSidecar(LinkBase):
    LinkRate = 200*MB
    AC1 = [ DLApp(pkt_size=PKT, arrival=1.25*MB, max_qos=2000, weight=5000) ]
    AC2 = [ ThruApp(min_thru=1*MB, weight=1) ]

AllLinks = [
    NormalLinkAndSidecar,
    NormalLinkAndDelay,
    NormalLink1, NormalLink2,
]



total_utility = 0
qos_function  = 0
constraints = []
## Problems
for link in AllLinks:
    link_utility = 0
    for aci, acq in enumerate(link.iter()):
        for app in acq:
            link_utility += app.utility / link.LinkRate
            qos_function += app.calc_qos(link, AllLinks)
            constraints.extend( app.constraints )
            pass
    constraints.append( link_utility <= 1.0 )
    total_utility += link_utility
## bandwidth utility constraints
constraints.extend([
    total_utility >= 0.6,
    total_utility <= 1.0
])

## Solution
objective = cp.Minimize( qos_function )
problem = cp.Problem(objective, constraints)
problem.solve()

## Print
for link in AllLinks:
    print(f'{link.name}: ')
    for aci, acq in enumerate(link.iter()):
        if acq:
            names = [app.name for app in acq if isinstance(app, ThruApp)]
            values = [app.variable.value/MB for app in acq if isinstance(app, ThruApp) ] # in unit MBps
            records = [f'{k}:{v:.3f}MBps' for k,v in zip(names, values)]
            print(f'\tAC-{aci}: {records}')
