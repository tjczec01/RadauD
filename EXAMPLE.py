# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 12:10:55 2020

Github: https://github.com/tjczec01

@author: Travis J Czechorski 

E-mail: tjczec01@gmail.com
"""

from __future__ import division, print_function, absolute_import
import os
dir_path = os.path.dirname(os.path.realpath(__file__))
exec(open(r'{}\__init__.py'.format(dir_path)).read())
from scipy.integrate import solve_ivp
from ivpd import solve_ivpd
from common import dm 

def flatten(lm):
              flatten = lambda l: [item for sublist in l for item in sublist] 
              return flatten(lm)

def RHS(t, y, args):
    pre = args
    
    try:
        yn = [dm(i, pre)**dm(3.0, pre) - dm(2.0, pre)*dm(i, pre)**dm(2.0, pre) - dm(8.0, pre)*dm(i, pre) + dm(1.0, pre) for i in y]
        return yn
    except:
        yn = dm(y, pre)**dm(3.0, pre) - dm(2.0, pre)*dm(y, pre)**dm(2.0, pre) - dm(8.0, pre)*dm(y, pre) + dm(1.0, pre)
        return [yn]


def jacob(t, y, args):
    pre = args
    return [[dm(3.0, pre)*(dm(y, pre)**dm(2.0, pre)) - dm(4.0, pre)*dm(y, pre) - dm(8.0, pre)]]

prec = [25]
tevs = [i/10 for i in range(0, 6, 1)]
sol1 = solve_ivp(RHS, [0.0, 0.5], [1.0], t_eval=tevs,  method="Radau", args=(prec), jac=jacob)
sol2 = solve_ivpd(RHS, [0.0, 0.5], [dm(1.0, prec)], t_eval=tevs,  method="RadauD", args=(prec), jac=jacob)
print(sol1.t.tolist())
print(sol1.y[0])
print(sol2.t)
print(flatten(sol2.y))
