#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 21:19:47 2019

@author: shoumikdc
"""

import numpy as np
import sympy as sp

from interruptingcow import timeout
bound = 10


def haar_measure(n, matrix_type, int_type = "float"):
    '''A random matrix distributed with Haar measure; see arXiv:math-ph/0609050
    
        Note that np.multiply(a,b) <--> a*b while np.matmul(a,b) <--> a @ b
        We use the former below as it is faster than converting ph to a matrix
        and then performing matmul; but it only works because ph is diagonal
        
        Also keep in mind; random.randn samples from Gaussian G(0,1)
    '''
    
    if int_type == "float":
        rng = lambda n: np.random.randn(n,n)
    elif int_type == "int":
        rng = lambda n: np.random.randint(-bound, bound, size = (n,n))
    
    
    if matrix_type == "U":
        # U for "unitary"
        z = (rng(n) + 1j*rng(n))/np.sqrt(2.0)
    elif matrix_type == "O":
        # O for "orthogonal" i.e. real
        z = rng(n)
    else:
        raise ValueError('incorrect input')
        
    q,r = np.linalg.qr(z)
    d = np.diag(r)
    ph = d/np.absolute(d)
    q = np.multiply(q,ph,q) # the extra q is storage location
    return q 


def generate_diag_list(n, int_type = "float"):
    ''' Generate an list of n numbers with random entries -- must be used with 
    np.diag() to get an n x n diagonal matrix '''
    
    if int_type == "float":
        rng = lambda : np.random.uniform(-bound, bound)
    elif int_type == "int":
        rng = lambda : np.random.randint(-bound, bound)
    
    
    M = []
    for i in range(n):
        x = rng()
        while x == 0:
            x = rng()
        M.append(x)
        
    return M

    
def generate_squeezing(n, int_type = "float"):
    ''' Generate an np array D of length 2n, of the form [M, 1/M] where M is a 
        list of nonzero random numbers. Easily extended to diag[M, 1/M] where 
        M is an n x n diagonal matrix <--> D is a 2n x 2n matrix
    '''
    
    M = generate_diag_list(n, int_type)

    M_inv = [1/y for y in M]
    
    D_list = M + M_inv # concatenate two lists
    D_array = np.asarray(D_list)
    return D_array
            

def generate_SP(n, int_type = "float"):
    ''' Generate a random 2n x 2n symplectic matrix S using the Pre-Iwasawa
        decomposition S = K D L where D is a squeezing matrix diag[M, 1/M], K 
        has form K = [[Id, 0], [P, Id]] where P is symmetric. L is orthogonal
        and symplectic (i.e. L = [[X,Y],[-Y,X]] where X+iY is unitary)
    '''
    
    # Generate K
    O = haar_measure(n, "O", int_type) # Orthogonal change of basis
    A = np.diag(generate_diag_list(n, int_type)) # random eigenvalues: need to convert 
                                       # list to diagonal ndarray first! 
    P = O.T @ A @ O # P = O^T A O
    Id = np.identity(n)
    Z = np.zeros(shape=(n,n))
    K = np.block([[Id,Z], [P,Id]])
    
    
    # Generate symplectic orthogonal matrix L
    U2 = haar_measure(n, "U", int_type)
    X, Y = U2.real, U2.imag
    L = np.block([[X,Y], [-Y,X]]) 
    
    # Generate squeezing matrix D = diag[B, 1/B] 
    D = np.identity(2*n) #np.diag(generate_squeezing(n, int_type))
    
    # matrix multiplication S = KDL
    S = K @ D @ L
    return S

def J_mat(n, basis = "qqpp"):
    ''' Generate the canonical 2n x 2n symplectic matric J
        i.e. so that (S^t)J(S) = J for any S in SP(2N)'''
    
    Id = np.identity(n)
    if basis == "qqpp":   
        Z = np.zeros(shape=(n,n))
        J = np.block([[Z, Id], [-Id, Z]])
    elif basis == "qpqp":
        w = [[0, 1], [-1, 0]]
        J = np.kron(Id, w)
    return J


def blockify(array, nrows, ncols):
    ''' Split a matrix into blocks of dimensions nrows x ncols; only works when 
    nrows/ncols divide array rows/cols. See https://stackoverflow.com/questions
    /11105375/how-to-split-a-matrix-into-4-blocks-using-numpy
    '''
    r, h = array.shape
    return (array.reshape(h//nrows, nrows, -1, 
                          ncols).swapaxes(1, 2).reshape(-1, nrows, ncols))





       
    
    
    
    
    
    