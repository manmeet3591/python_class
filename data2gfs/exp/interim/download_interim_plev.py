#!/usr/bin/python
from ecmwfapi import ECMWFDataServer
import sys

server = ECMWFDataServer()

idate = sys.argv[1]
if len(idate) != 10:
    raise ValueError, "get_erainterim.py idate grid"
grid = sys.argv[2]
fname = sys.argv[3]

server.retrieve({
        'dataset' : "interim",
        'levtype' : 'pl',
        'date'    : idate[:8],
        'time'    : idate[-2:],
        'levelist': "1000/975/950/925/900/875/850/825/800/775/750/700/650/600/550/500/450/400/350/300/250/225/200/175/150/125/100/70/50/30/20/10/7/5/3/2/1",
        'type'    : "an",
        'param'   : 'z/u/v/t/q/o3/clwc/ciwc',
        'grid'    : '{0}/{0}'.format(grid),
        'area'    : "90/0/-90/360",
        'target'  : fname
        })

