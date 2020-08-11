#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 11 20:26:41 2019

@author: shoumikdc
"""

from random_matrix import *

#def find_local_ops():
''' decoupling protocol for the 3 mode case'''
S = generate_SP(3) 
idx = [0,3,1,4,2,5]
S = (((S.T)[idx]).T)[idx] # change to qpqpqp basis

v1 = S[0]
u1 = S.T[0]


""" 
def grobner_test_3():
    '''Main routine: generate 3-mode SP matrix and perform matrix operations. 
        We have S diag(X Y Z) S and want to decouple the three modes by setting 
        certain entries to zero.'''
    num_inconsistent = 0
    num_not_inconsistent = 0
    num_timeouts = 0
    total = 50
    
    for i in range(total):
        Sym = generate_SP(3) # generates 6 x 6 SP matrix with columns qqqppp
        assert Sym.shape == (6,6) # just to check
        
        idx = [0,3,1,4,2,5] # permutation of columns to get qpqpqp
        Sym = (((Sym.T)[idx]).T)[idx] # change of basis // permute cols then rows
        
        blocks = blockify(Sym, 2, 2) # split into subblocks
        A, B, C, D, E, F, G, H, I = (sp.Matrix(x) for x in blocks) # SymPy time! 
        
        # variable matrices
        variables = sp.symbols("a b c d e f g h k l m n")
        a, b, c, d, e, f, g, h, k, l, m, n = variables
        X = sp.Matrix([[a,b],[c,d]])
        Y = sp.Matrix([[e,f],[g,h]])
        Z = sp.Matrix([[k,l],[m,n]])
        
        R1 = A @ X @ C + B @ Y @ F + C @ Z @ I
        R2 = D @ X @ C + E @ Y @ F + F @ Z @ I
        
        eq1, eq2, eq3, eq4 = R1[0,0], R1[0,1], R1[1,0], R1[1,1]
        eq5, eq6, eq7, eq8 = R2[0,0], R2[0,1], R2[1,0], R2[1,1]
        eq9 = a*d - b*c - 1
        eq10 = e*h - f*g - 1
        eq11 = k*n - l*m - 1
        
        equations = [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11]
        
        # check if grobner basis is [1] or timeout 
        try:
            with timeout(1, exception=RuntimeError):
                G = sp.groebner(equations, variables)
                if (G == [1]):
                    num_inconsistent += 1
                else:
                    num_not_inconsistent += 1
        except RuntimeError:
            num_timeouts += 1
            pass
    
    print("Out of", total, "total trials:")
    print(num_inconsistent, "trials had inconsistent solutions with G == [1]")
    print(num_not_inconsistent, "trials had G != [1], so possibly consistent")
    print(num_timeouts, "trials took too long!")
"""