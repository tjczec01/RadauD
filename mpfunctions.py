# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 19:39:06 2020

Github: https://github.com/tjczec01

@author: Travis J Czechorski 

E-mail: tjczec01@gmail.com
"""
from __future__ import division, print_function, absolute_import
import mpmath as mp
from mpmath import fdot, mpf, mpc
import os
import numpy as np

# For illustrative purposes.

__name__ = "mpfunctions"

ppath = os.path.dirname(os.path.abspath(__file__))
pth = r'{}'.format(ppath)

__all__ = ["MF", "mat"] 


def mv(v, ds):
        try:
            mp.dps = int(ds)
        except:
            mp.dps = int(ds[0])
        try:
               mv = mpf(v)
        except:
            try:
               mv = mpc(v.real, v.imag)
            except:
                try:
                    mv = mpf(v[0])
                except:
                    mv = mpc(v[0].real, v[0].imag)
        return mv
   
flatten = lambda l: [item for sublist in l for item in sublist]

def normd(A, prec):
       rows = len(A)
       cols = len(A[0])
       vals = []
       for i in range(rows):
              
              for j in range(cols):
                     vi = mv(abs(A[i][j]), prec)**mv(2, prec)
                     vals.append(vi)
                     
       vf = mv(sum(vals), prec)**mv((1/2), prec)
       return vf

class MF:
       
       def __init__(self, Ax, pr):
              self.Ax = Ax
              try:
                  self.pr = int(pr)
              except:
                  self.pr = int(pr[0])
              self.__name__ = "MF"
              
       
       flatten = lambda l: [item for sublist in l for item in sublist]
       
       def mmake(self, v):
             try:
                 mp.dps = int(self.pr)
             except:
                 mp.dps = int(self.pr[0])
             try:
                    mv = mpf(v)
             except:
                 try:
                    mv = mpc(v.real, v.imag)
                 except:
                     try:
                         mv = mpf(v[0])
                     except:
                         mv = mpc(v[0].real, v[0].imag)
             return mv     
       
       
       def mmakec(self, vs):
           mp.dps = int(self.pr)
           vr = []
           vc = []
           try:
                len(vs)
                VS = vs
           except:
                VS = [vs]
           for i in range(len(VS)):
                vr.append(VS[i].real)
                vc.append(VS[i].imag)
           try:
                 ivs = len(vr) 
                 mv = [mp.mpc(vr[i], vc[i]) for i in range(ivs)]
                 return flatten(mv)
           except:
                 mv = mp.mpc(vr[0], vc[0])
                 return mv 
       
       def TD(self):
              columns = len(self.Ax[0][:])
              rows = len(self.Ax)
              vf = self.mmake(0.0)
              tmat = [[vf for i in range(rows)] for j in range(columns)]
              for i in range(rows):
                     for j in range(columns):
                            vvf = self.mmake(self.Ax[i][j])
                            tmat[j][i] = vvf
              return tmat
                     
       def T(self, A):
              columns = len(A[0][:])
              rows = len(A)
              vf = self.mmake(0.0)
              tmat = [[vf for i in range(rows)] for j in range(columns)]
              for i in range(rows):
                     for j in range(columns):
                            vvf = self.mmake(A[i][j])
                            tmat[j][i] = vvf
              return tmat
       
       def dotdc(self, v1, v2):
              vv = fdot(v1, v2) 
              return vv
         
       def dotd(self, v1, v2):
              vv = fdot(v1, v2)
              return vv   
       
       def zeros_matrix(self, rows, cols):
           """
           Creates a matrix filled with zeros.
               :param rows: the number of rows the matrix should have
               :param cols: the number of columns the matrix should have
               :return: list of lists that form the matrix
           """
           M = []
           while len(M) < rows:
               M.append([])
               while len(M[-1]) < cols:
                   M[-1].append(self.mmake(0.0))
       
           return M
       
       def zerod(self, n):
              mm = []
              for ni in range(n):
                      vf = self.mmake(0.0)
                      mm.append([vf for mi in range(n)])
              
              for nn in range(n):
                      vfb = self.mmake(0.0)
                      mm[nn][nn] = vfb
              return mm 
       
       
       def LU_decompositiond(self):
           """Perform LU decomposition using the Doolittle factorisation."""
       
           N = len(self.Ax)
           L = self.zerod(N)
           U = self.zerod(N)
           
           
           def uvals(Um, k, n):
                   ulist = []
                   for i in range(k):
                         uu = Um[i][n]
                         ulist.append(uu)
                   return ulist
            
           def lvals(Lm, k, n):
                   lu = Lm[k]
                   lul = lu[0:k]
                   return lul
           
           def uvalsd(Um, k, n):
                   ulist = []
                   for i in range(k):
                         uu = Um[i][n]
                         ulist.append(uu)
                   return ulist
            
           def lvalsd(Lm, k, n):
                   llist = []
                   lu = Lm[k]
                   lul = lu[0:k]
                   for i in range(len(lul)):
                            val_ij = lul[i]
                            llist.append(val_ij)
                   return lul
              
           for k in range(N):
               try:
                    v1 = self.mmake(1.0)
                    L[k][k] = v1
                    v2 = self.mmake((self.mmake(self.Ax[k][k]) - self.mmake(self.dotd(lvalsd(L, k, k), uvalsd(U, k, k)))) / self.mmake(L[k][k]))
               except:
                   v1 = self.mmakec(complex(1.0, 0.0))
                   L[k][k] = v1
                   v2b = self.dotdc(lvalsd(L, k, k), uvalsd(U, k, k))
                   v2c = self.mmakec(self.Ax[k][k]) - self.mmakec(v2b)
                   v2 = self.mmakec(v2c) / self.mmakec(L[k][k])
               U[k][k] = v2
               try:
                    
                    for j in range(k+1, N):
                          val_i = self.mmake(self.mmake((self.mmake(self.Ax[k][j]) - self.mmake(self.dotd(lvalsd(L, k, k), uvalsd(U, j, j)))) / self.mmake(L[k][k])))
                          U[k][j] = val_i
                    for i in range(k+1, N):
                          val_ib = self.mmake(self.mmake((self.mmake(self.Ax[i][k]) - self.mmake(self.dotd(lvalsd(L, i, i), uvalsd(U, k, k)))) / self.mmake(U[k][k])))
                          L[i][k] = val_ib
               except:
                    for j in range(k+1, N):
                     val_i = self.mmakec(self.mmakec((self.mmakec(self.Ax[k][j]) - self.mmakec(self.dotdc(lvalsd(L, k, k), uvalsd(U, j, j)))) / self.mmakec(L[k][k])))
                     U[k][j] = val_i
                    for i in range(k+1, N):
                     val_ib = self.mmakec(self.mmakec((self.mmakec(self.Ax[i][k]) - self.mmakec(self.dotdc(lvalsd(L, i, i), uvalsd(U, k, k)))) / self.mmakec(U[k][k])))
                     L[i][k] = val_ib
       
           return L, U
       
       def LU_decompositionf(self):
           """Perform LU decomposition using the Doolittle factorisation."""
       
           N = len(self.Ax)
           L = self.zerod(N)
           U = self.zerod(N)
           
           
           def uvals(Um, k, n):
                   ulist = []
                   for i in range(k):
                         uu = Um[i][n]
                         ulist.append(uu)
                   return ulist
            
           def lvals(Lm, k, n):
                   lu = Lm[k]
                   lul = lu[0:k]
                   return lul
           
           def uvalsd(Um, k, n):
                   ulist = []
                   for i in range(k):
                         uu = Um[i][n]
                         ulist.append(uu)
                   return ulist
            
           def lvalsd(Lm, k, n):
                   llist = []
                   lu = Lm[k]
                   lul = lu[0:k]
                   for i in range(len(lul)):
                            val_ij = lul[i]
                            llist.append(val_ij)
                   return lul
              
           for k in range(N):
               v1 = 1.0
               L[k][k] = v1
               v2 = (self.Ax[k][k] - self.dotd(lvalsd(L, k, k), uvalsd(U, k, k))) / L[k][k]
               U[k][k] = v2
               for j in range(k+1, N):
                     val_i = (self.Ax[k][j] - self.dotd(lvalsd(L, k, k), uvalsd(U, j, j))) / L[k][k]
                     U[k][j] = val_i
               for i in range(k+1, N):
                     val_ib = (self.Ax[i][k] - self.dotd(lvalsd(L, i, i), uvalsd(U, k, k))) / U[k][k]
                     L[i][k] = val_ib
       
           return L, U 
       
       def forward_subd(self, L, b):
           """Given a lower triangular matrix L and right-side vector b,
           compute the solution vector y solving Ly = b."""
       
           y = []
           for i in range(len(b)):
               y.append(b[i])
               for j in range(i):
                     v1 = self.mmake(self.mmake(L[i][j])*self.mmake(y[j]))
                     v2 = self.mmake(y[i]) 
                     v3 = self.mmake(v2 - v1)
                     y[i]= v3
               y[i] = self.mmake(self.mmake(y[i])/self.mmake(L[i][i]))
       
           return y
       
       def backward_subd(self, U, y):
           """Given a lower triangular matrix U and right-side vector y,
           compute the solution vector x solving Ux = y."""
            
           x = [self.mmake(0.0) for ix in y]
           for i in range(len(x)-1, 0, -1):
               val_i = self.mmake((self.mmake(y[i-1]) - self.mmake(self.dotd(U[i-1][i:], x[i:]))) / self.mmake(U[i-1][i-1]))
               x[i-1] = self.mmake(val_i)
              
           return x
       
       def forward_sub(self, L, b):
           """Given a lower triangular matrix L and right-side vector b,
           compute the solution vector y solving Ly = b."""
           y = b
           try:
                
                for i in range(len(b)):
                     for j in range(i):
                          y[i] = self.mmake(self.mmake(y[i]) - self.mmake(self.mmake(L[i][j])*self.mmake(y[j])))
                y[i] = self.mmake(self.mmake(y[i])/self.mmake(L[i][i]))
                
           except:
              for i in range(len(b)):
                for j in range(i):
                     y[i] = self.mmakec(self.mmakec(y[i]) - self.mmakec(self.mmakec(L[i][j])*self.mmakec(y[j])))
                y[i] = self.mmakec(self.mmakec(y[i])/self.mmakec(L[i][i]))  
           return y
       
       def backward_sub(self, U, y):
           """Given a lower triangular matrix U and right-side vector y,
           compute the solution vector x solving Ux = y."""
            
           x = [self.mmake(0.0) for ix in y]
           XX = [self.mmakec(complex(0.0, 0.0)) for ix in y]
           for i in range(len(x)-1, -1, -1):
               try:
                        val_i = self.mmake(y[i-1] - self.dotd(U[i-1][i-1:], x[i-1:])) / self.mmake(U[i-1][i-1])
               except:
                    try:
                         try:
                              val_i = self.mmake(y[i-1] - self.dotd(U[i][i:], x[i:])) / self.mmake(U[i-1][i-1])
                         except:
                              try:
                                   val_i = self.mmake(y[i-1] - self.dotd(U[i][i-1:], x[i:])) / self.mmake(U[i-1][i-1])
                              except:
                                   val_i = self.mmake(y[i-1] - self.dotd(U[i][i:], x[i-1:])) / self.mmake(U[i-1][i-1])
                    except:
                         try:
                              val_i = (self.mmakec(y[i-1]) - self.dotdc(U[i-1][i-1:], XX[i-1:])) / self.mmakec(U[i-1][i-1])
                         except:
                              try:
                                   val_i = (self.mmakec(y[i-1]) - self.dotdc(U[i][i:], XX[i:])) / self.mmakec(U[i-1][i-1])
                              except:
                                   try:
                                        val_i = (self.mmakec(y[i-1]) - self.dotdc(U[i][i:], XX[i:])) / self.mmakec(U[i-1][i-1])
                                   except:
                                        val_i = (self.mmakec(y[i-1]) - self.dotdc(U[i][i:], XX[i:])) / self.mmakec(U[i-1][i-1])
               x[i-1] = val_i
           return x
       
       def lu_solved(self, L, U, b):
           # Step 1: Solve Uy = b using forward substitution
           # Step 2: Solve Lx = y using backward substitution
           y = self.forward_sub(L, b)
           x = self.backward_sub(U, y)
           return x
       
       def linear_solved(self, bd):
           Ld, Ud = self.LU_decompositionf()
           x = self.lu_solved(Ld, Ud, bd)
           return x
       
       def mfunc(self):
              rows = len(self.Ax)
              cols = len(self.Ax[0])
              AD = [[self.mmake(ij) for ij in self.Ax[ii]] for ii in range(len(self.Ax))]
              for i in range(rows):
                     for j in range(cols):
                            val_ij = self.mmake(self.Ax[i][j])
                            AD[i][j] = val_ij
              return AD
        
       def mfuncb(self, Ax):
              rows = len(Ax)
              cols = len(Ax[0])
              AD = [[self.mmake(ij) for ij in Ax[ii]] for ii in range(len(Ax))]
              for i in range(rows):
                     for j in range(cols):
                            val_ij = self.mmake(Ax[i][j])
                            AD[i][j] = val_ij
              return AD
       
       def mfuncl(self, Axl):
              vals = len(Axl)
              ADl = [self.mmake(il) for il in Axl]
              for i in range(vals):
                     val_i = self.mmake(Axl[i])
                     ADl[i] = val_i
              return ADl
       
       def matrix_multiplyd(self, A, B, pr):
           """
           Returns the product of the matrix A * B
               :param A: The first matrix - ORDER MATTERS!
               :param B: The second matrix
               :return: The product of the two matrices
           """
           # Section 1: Ensure A & B dimensions are correct for multiplication
           rowsA = len(A)
           colsA = len(A[0])
           rowsB = len(B)
           colsB = len(B[0])
           if colsA != rowsB:
               raise ArithmeticError(
                   'Number of A columns must equal number of B rows.')
       
           # Section 2: Store matrix multiplication in a new matrix
           C = self.zeros_matrix(rowsA, colsB)
           for i in range(rowsA):
               for j in range(colsB):
                   total = self.mmake(0)
                   for ii in range(colsA):
                       total += self.mmake(self.mmake(A[i][ii]) * self.mmake(B[ii][j]))
                   C[i][j] = self.mmake(total)
       
           return C
    
       def matrix_multiplydf(self, A, B, pr):
           """
           Returns the product of the matrix A * B
               :param A: The first matrix - ORDER MATTERS!
               :param B: The second matrix
               :return: The product of the two matrices
           """
           # Section 1: Ensure A & B dimensions are correct for multiplication
           rowsA = len(A)
           colsA = len(A[0])
           rowsB = len(B)
           colsB = len(B[0])
           if colsA != rowsB:
               raise ArithmeticError(
                   'Number of A columns must equal number of B rows.')
       
           # Section 2: Store matrix multiplication in a new matrix
           C = self.zeros_matrix(rowsA, colsB)
           for i in range(rowsA):
               for j in range(colsB):
                   total = 0
                   for ii in range(colsA):
                       total += float(A[i][ii]) * B[ii][j]
                   C[i][j] = total
       
           return C   
    
       def matrixadd(self, A, B):
                     Z = []
                     for i in range(len(A)):
                         row = []
                         for j in range(len(A[i])):
                             row.append(self.mmake(self.mmake(A[i][j]) + self.mmake(B[i][j])))
                         Z.append(row)
                         
                     return Z
       def matrixaddc(self, A, B):
                     Z = []
                     for i in range(len(A)):
                         row = []
                         for j in range(len(A[i])):
                             row.append(self.mmakec(A[i][j].real, A[i][j].imag) + self.mmakec(B[i][j].real, B[i][j].imag))
                         Z.append(row)
                         
                     return Z
              
       def matrixsub(self, A, B):
                     Z = []
                     for i in range(len(A)):
                         row = []
                         for j in range(len(A[i])):
                             row.append(self.mmake(self.mmake(A[i][j]) - self.mmake(B[i][j])))
                         Z.append(row)
                         
                     return Z
       
       
       def pivot_matrix(self, M):
           """Returns the pivoting matrix for M, used in Doolittle's method."""
           m = len(M)
           MM = m - 2
           # Create an identity matrix, with floating point values                                                                                                                                                                                            
           id_mat = [[float(i==j) for i in range(m)] for j in range(m)]
           ID = []
           r = [0]
           # Rearrange the identity matrix such that the largest element of                                                                                                                                                                                   
           # each column of M is placed on the diagonal of of M  
           for j in range(m):
                  
              rowv = max(range(j, m), key=lambda i: abs(M[i][j]))
              rv = max(r)
              if j == MM:
                     rowv = rowv + 1
              else:
                     pass
              if rowv > rv:
                     r.append(rowv)
              else:
                     pass
                   # Swap the rows                                                                                                                                                                                                                            
              id_mat[j], id_mat[max(r)] = id_mat[max(r)], id_mat[j]
              ID.append(id_mat[j])
              
           return ID
       
       def plud(self):
       
           """Performs an LU Decomposition of A (which must be square)                                                                                                                                                                                        
           into PA = LU. The function returns P, L and U."""
           n = len(self.Ax)
       
           # Create zero matrices for L and U                                                                                                                                                                                                                 
           L = [[self.mmake(0.0)] * n for i in range(n)]
           Lb = [[self.mmake(0.0)] * n for i in range(n)]
           U = [[self.mmake(0.0)] * n for i in range(n)]
           # Create the pivot matrix P and the multipled matrix PA 
           PP = self.pivot_matrix(self.Ax)                                                                                                                                                                                           
           P = self.decfuncb(PP)
           try:
                  PA = self.matrix_multiplyd(P, self.Ax, self.pr)
           except:
                  PA = self.matrix_multiplydf(P, self.Ax, self.pr)
       
           # Perform the LU Decomposition                                                                                                                                                                                                                     
           for j in range(n):
               # All diagonal entries of L are set to unity                                                                                                                                                                                                   
              L[j][j] = self.mmake(1.0)
              Lb[j][j] = self.mmake(0.0) 
       
              # LaTeX: u_{ij} = a_{ij} - \sum_{k=1}^{i-1} u_{kj} l_{ik}                                                                                                                                                                                      
              for i in range(j+1):
                     s1 = self.mmake(sum(self.mmake(U[k][j]) * self.mmake(L[i][k]) for k in range(i)))
                     U[i][j] = self.mmake(self.mmake(PA[i][j]) - self.mmake(s1))
       
              # LaTeX: l_{ij} = \frac{1}{u_{jj}} (a_{ij} - \sum_{k=1}^{j-1} u_{kj} l_{ik} )                                                                                                                                                                  
              for i in range(j, n):
                  s2 = self.mmake(sum(self.mmake(U[k][j]) * self.mmake(L[i][k]) for k in range(j)))
                  L[i][j] = self.mmake((self.mmake(PA[i][j]) - self.mmake(s2)) / self.mmake(U[j][j]))
                  Lb[i][j] = self.mmake((self.mmake(PA[i][j]) - self.mmake(s2)) / self.mmake(U[j][j]))
       
       
           for j in range(n):
               # All diagonal entries of L are set to unity                                                                                                                                                                                                   
              Lb[j][j] = self.mmake(0.0) 
           return (PP, P, L, U, Lb)
    
       
       def pluf(self):
       
           """Performs an LU Decomposition of A (which must be square)                                                                                                                                                                                        
           into PA = LU. The function returns P, L and U."""
           n = len(self.Ax)
       
           # Create zero matrices for L and U                                                                                                                                                                                                                 
           L = [[self.mmake(0.0)] * n for i in range(n)]
           Lb = [[self.mmake(0.0)] * n for i in range(n)]
           U = [[self.mmake(0.0)] * n for i in range(n)]
           # Create the pivot matrix P and the multipled matrix PA 
           PP = self.pivot_matrix(self.Ax)                                                                                                                                                                                           
           P = self.decfuncb(PP)
           try:
                  PA = self.matrix_multiplyd(P, self.Ax, self.pr)
           except:
                  PA = self.matrix_multiplydf(P, self.Ax, self.pr)
       
           # Perform the LU Decomposition                                                                                                                                                                                                                     
           for j in range(n):
               # All diagonal entries of L are set to unity                                                                                                                                                                                                   
              L[j][j] = self.mmake(1.0)
              Lb[j][j] = self.mmake(0.0) 
       
              # LaTeX: u_{ij} = a_{ij} - \sum_{k=1}^{i-1} u_{kj} l_{ik}                                                                                                                                                                                      
              for i in range(j+1):
                     s1 = sum(self.mmake(U[k][j]) * self.mmake(L[i][k]) for k in range(i))
                     U[i][j] = self.mmakec(PA[i][j].real, PA[i][j].imag) - self.mmake(s1)
       
              # LaTeX: l_{ij} = \frac{1}{u_{jj}} (a_{ij} - \sum_{k=1}^{j-1} u_{kj} l_{ik} )                                                                                                                                                                  
              for i in range(j, n):
                  s2 = sum(self.mmake(U[k][j]) * self.mmake(L[i][k]) for k in range(j))
                  L[i][j] = (self.mmakec(PA[i][j].real, PA[i][j].imag) - self.mmake(s2)) / self.mmakec(U[j][j].real, U[j][j].imag)
                  Lb[i][j] = (self.mmakec(PA[i][j].real, PA[i][j].imag) - self.mmake(s2)) / self.mmakec(U[j][j].real, U[j][j].imag)
       
       
           for j in range(n):
               # All diagonal entries of L are set to unity                                                                                                                                                                                                   
              Lb[j][j] = self.mmake(0.0) 
           return (PP, P, L, U, Lb)
       
       def LU_factor(self):
              PP, P, L, U, Lb = self.plud()
              Pp = [[self.mmake(P[i][j]) for i in range(len(PP[0]))] for j in range(len(PP))]
              LU = self.matrixadd(Lb, U)
              PIV = self.piv(Pp)
              return LU, PIV
       
       def LU_factorc(self):
              PP, P, L, U, Lb = self.pluf()
              Pp = [[self.mmake(P[i][j]) for i in range(len(PP[0]))] for j in range(len(PP))]
              LU = self.matrixaddc(Lb, U)
              PIV = self.piv(Pp)
              return LU, PIV
       
       def printm(self, Mm):
           """
           Print a matrix one row at a time
               :param M: The matrix to be printed
           """
           for row in Mm:
                try:
                    print([self.mmake(x) for x in row])
                except:
                     print([self.mmakec(x) for x in row])
    
       def printmf(self, Mm):
           """
           Print a matrix one row at a time
               :param M: The matrix to be printed
           """
           for row in Mm:
                try:
                    print([self.mmake(x) for x in row])
                except:
                     print([self.mmakec(x) for x in row])
               
       def csc_groups(self):
              A = self.Ax
              rows = len(A)
              cols = len(A[0])
              rows_l = []
              cols_l = []
              data_l = []
              for i in range(rows):
                     for j in range(cols):
                            val_ij = self.mmake(A[i][j])
                            if val_ij == 0 or val_ij == 0.0:
                                   pass
                            elif val_ij != 0 or val_ij != 0.0:
                                   rows_l.append(int(i))
                                   cols_l.append(int(j))
                                   data_l.append(self.mmake(A[i][j]))
              return rows_l, cols_l, data_l
       
       def csc_array(self, rows_l, cols_l, data_l):
              rs = max(rows_l)
              cs = max(cols_l)
              Z = self.zeros_matrix(rs+1, cs+1)
              k = 0
              for i, j in zip(rows_l, cols_l):
                     val = data_l[k]
                     k += 1
                     Z[i][j] = self.mmake(val)
                     
              return Z
       
       def piv(self, P):
              PTb = self.T(P)
              PT = PTb
              h = []
              Pr = len(PT)
              for i in range(Pr):
                     l = PT[i]
                     n = l.index(1.0)
                     h.append(int(n))
                     
              return h
       
       def solve(self, A, B):
              y = self.forward_sub(A, B)
              x = self.backward_sub(A, y)
              
              
              return x
              
       
       def factor_solve(self, B):
              LU, PIV = self.LU_factor()
              A = self.Ax
              AA = []
              for i in PIV:
                     AA.append(A[i])
              x = self.solve(AA, B)
              return x
       
       def factor_solvec(self, B):
              LU, PIV = self.LU_factorc()
              A = self.Ax
              AA = []
              for i in PIV:
                     AA.append(A[i])
              x = self.solve(AA, B)
              return x
       
       def ufsub(self, L, b):
           """ Unit row oriented forward substitution """
           for i in range(len(L[0])): 
               for j in range(i):
                   b[i] = self.mmake(b[i]) - self.mmake(L[i][j])*self.mmake(b[j])
           return b

       def bsub(self, U, y):
           """ Row oriented backward substitution """
           Us, US = np.array(U).shape
           for i in range(Us-1, 0, -1): 
               for j in range(i+1, US):
                   y[i] -= self.mmake(U[i][j])*self.mmake(y[j])
               y[i] = self.mmake(y[i])/self.mmake(U[i][i])
           return y
    
       def solves(self, LU, B):
              y = self.ufsub(LU, B)
              x = self.bsub(LU, y)
              
              return x
       
       def iden(self, n):
              mm = []
              for ni in range(n):
                     mm.append([self.mmake(0.0) for mi in range(n)])
              
              for nn in range(n):
                     mm[nn][nn] = self.mmake(1.0)
              return mm
       
       def idenf(self, n):
              mm = []
              for ni in range(n):
                     mm.append([self.mmake(0.0) for mi in range(n)])
              
              for nn in range(n):
                     mm[nn][nn] = self.mmake(1.0)
              return mm
       
       def idenv(self, n, v=1.0):
              mm = []
              for ni in range(n):
                     mm.append([self.mmake(0.0) for mi in range(n)])
              
              for nn in range(n):
                     mm[nn][nn] = self.mmake(1.0)*self.mmake(v)
              return mm
       
       def idenfv(self, n, v=1.0):
              mm = []
              for ni in range(n):
                     mm.append([self.mmake(0.0) for mi in range(n)])
              
              for nn in range(n):
                     mm[nn][nn] = self.mmake(1.0)*v
              return mm
         
       def idenvc(self, n, v=complex(1.0, 0.0)):
              mm = []
              for ni in range(n):
                     mm.append([self.mmakec(0.0) for mi in range(n)])
              
              for nn in range(n):
                     mm[nn][nn] = self.mmakec(1.0)*self.mmakec(v)
              return mm
              
                     
       
class mat:
    def __init__(self, M, pr):
        self.pr = pr
        self.M = M

    @property
    def T(self):
        try:
            return self._cached_T  # attempt to read cached attribute
        except AttributeError:
            self._cached_T = self._transpose()  # compute and cache
            return self._cached_T

    def _transpose(self):
        M = self.M
        rows = len(M)
        cols = len(M[0])
        if not isinstance(M[0], list):
            M = [M]
        rows = len(M)
        cols = len(M[0])
        def zeros_matrix(rows, cols):
                  M = []
                  while len(M) < rows:
                      M.append([])
                      while len(M[-1]) < cols:
                          M[-1].append(MF(self.M, self.pr).mmake(0.0))
                  return M
        MT = zeros_matrix(cols, rows)
        for i in range(rows):
            for j in range(cols):
                MT[j][i] = MF(self.M, self.pr).mmake(M[i][j])
        return MT 

MF.__name__ = "MF" 
mat.__name__ = "mat"  
