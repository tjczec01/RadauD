# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 19:39:06 2020

Github: https://github.com/tjczec01

@author: Travis J Czechorski 

E-mail: tjczec01@gmail.com
"""
from __future__ import division, print_function, absolute_import
import mpmath as mp
from mpmath import *
import os
import importlib.util
import numpy as np

# For illustrative purposes.

__name__ = "mpfunctions"

ppath = os.path.dirname(os.path.abspath(__file__))
pth = r'{}'.format(ppath)

__all__ = ["MF", "mat"] 

# De = decimal.Decimal

def mv(v, ds):
        mp.dps = int(ds)
        mv = mp.mpf(v)
        # print(mp.prec)
        return mv

# def dm(v, pr):
       # De = decimal.Decimal
       # val_i = De(v)
       # aa = '{}:.{}f{}'.format('{', pr,'}')
       # aa1 = '{}'.format(aa)
       # aa2 = str(aa1)
       # aa3 = str(aa2.format(val_i))
       # aa4 = De(aa3)
       # return aa4

# print(mv(4.555, 50))
# print(mp.prec)
# class ComplexDecimal(object):
#             def __init__(self, value, prec):
#                 self.real = dm(value.real, prec)
#                 self.imag = dm(value.imag, prec)
#                 self.vi = float(value.imag)
#                 self.prec = prec
#                 if self.vi >= 0 or self.vi >= 0.0:
#                        self.sign = '+'
#                 elif self.vi <= 0 or self.vi <= 0.0:
#                        self.sign = '-'
#             def __add__(self, other):
#                 result = ComplexDecimal(self)
#                 result.real += dm(other.real, self.prec)
#                 if other.imag >= 0 or other.imag >= 0.0:
#                        result.imag += dm(other.imag, self.prec)
#                 elif other.imag <= 0 or other.imag <= 0.0:
#                        result.imag -= dm(other.imag, self.prec)
#                 return result
        
#             __radd__ = __add__
#             def __str__(self):
#                       return f'({self.real} {self.sign} {abs(self.imag)}j)'
        
#             def sqrt(self):
#                 result = ComplexDecimal(self)
#                 if self.imag:
#                     raise NotImplementedError
#                 elif self.real > 0:
#                     result.real = self.real.sqrt()
#                     return result
#                 else:
#                     result.imag = (-self.real).sqrt()
#                     result.real = dm(0, self.prec)
#                     return result
             
#             def real(self):
#                 result = ComplexDecimal(self)
#                 result.real = self.real
#                 return result
         
#             def imag(self):
#                 result = ComplexDecimal(self)
#                 result.imag = self.imag
#                 return result

# class ComplexDecimalr(object):
#             def __init__(self, valuer, valuei, prec):
#                 self.real = dm(valuer, prec)
#                 self.imag = dm(valuei, prec)
#                 self.vi = float(valuei)
#                 self.prec = prec
#                 if self.vi >= 0 or self.vi >= 0.0:
#                        self.sign = '+'
#                 elif self.vi <= 0 or self.vi <= 0.0:
#                        self.sign = '-'
#             def __add__(self, other):
#                 result = ComplexDecimal(self)
#                 result.real += dm(self.real, self.prec)
#                 if self.vi >= 0 or self.vi >= 0.0:
#                        result.imag += dm(other.imag, self.prec)
#                 elif self.vi <= 0 or self.vi <= 0.0:
#                        result.imag -= dm(other.imag, self.prec)
#                 return result
        
#             __radd__ = __add__
#             def __str__(self):
#                       return f'({self.real} {self.sign} {abs(self.imag)}j)'
        
#             def sqrt(self):
#                 result = ComplexDecimal(self)
#                 if self.imag:
#                     raise NotImplementedError
#                 elif self.real > 0:
#                     result.real = self.real.sqrt()
#                     return result
#                 else:
#                     result.imag = (-self.real).sqrt()
#                     result.real = dm(0, self.prec)
#                     return result
             
#             def real(self):
#                 result = ComplexDecimal(self)
#                 result.real = self.real
#                 return result
         
#             def imag(self):
#                 result = ComplexDecimal(self)
#                 result.imag = self.imag
#                 return result

# def DCM(MUD, MUD2, prec):
#        M1R = ComplexDecimal(MUD, prec).real
#        M1C = ComplexDecimal(MUD, prec).imag
#        M2R = ComplexDecimal(MUD2, prec).real
#        M2C = ComplexDecimal(MUD2, prec).imag
#        vv = M1R*M2R
#        vv2 = M1C*M2R
#        vv3 = M1R*M2C
#        vv4 = M1C*M2C
#        vvi = vv - vv4
#        vvj = vv3 + vv2
#        # return ComplexDecimalr(vvi, vvj, prec)
#        return ComplexDecimal(complex(vvi, vvj), prec)

# def DCA(MUD, MUD2, prec):
#        M1R = ComplexDecimal(MUD, prec).real
#        M1C = ComplexDecimal(MUD, prec).imag
#        M2R = ComplexDecimal(MUD2, prec).real
#        M2C = ComplexDecimal(MUD2, prec).imag
#        vv = M1R + M2R
#        vv2 = M1C + M2C
#        return ComplexDecimal(complex(vv, vv2), prec)
#        # return ComplexDecimalr(vv, vv2, prec)

# def DCS(MUD, MUD2, prec):
#        M1R = ComplexDecimal(MUD, prec).real
#        M1C = ComplexDecimal(MUD, prec).imag
#        M2R = ComplexDecimal(MUD2, prec).real
#        M2C = ComplexDecimal(MUD2, prec).imag
#        vv = M1R - M2R
#        vv2 = M1C - M2C
#        # return ComplexDecimalr(vv , vv2, prec)
#        return ComplexDecimal(complex(vv, vv2), prec)

# def DCD(MUD, MUD2, prec):
#        M3R = ComplexDecimal(MUD2, prec).real
#        M3C = -ComplexDecimal(MUD2, prec).imag
#        bot = DCM(MUD2, complex(M3R, M3C), prec).real
#        top = DCM(MUD, complex(M3R, M3C), prec)
#        TR = top.real / bot
#        TC = top.imag / bot
#        return ComplexDecimal(complex(TR, TC), prec)

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
              
       # class ComplexDecimal(object):
       #      def __init__(self, value, prec):
       #          self.real = dm(value.real, prec)
       #          self.imag = dm(value.imag, prec)
       #          self.vi = float(value.imag)
       #          self.prec = prec
       #          if self.vi >= 0 or self.vi >= 0.0:
       #                 self.sign = '+'
       #          elif self.vi <= 0 or self.vi <= 0.0:
       #                 self.sign = '-'
       #      def __add__(self, other):
       #          result = ComplexDecimal(self)
       #          result.real += dm(other.real, self.prec)
       #          if other.imag >= 0 or other.imag >= 0.0:
       #                 result.imag += dm(other.imag, self.prec)
       #          elif other.imag <= 0 or other.imag <= 0.0:
       #                 result.imag -= dm(other.imag, self.prec)
       #          return result
        
       #      __radd__ = __add__
       #      # vi = value.imag
       #      # if vi >= 0 or vi >= 0.0:
       #      def __str__(self):
       #                return f'({self.real} {self.sign} {abs(self.imag)}j)'
       #      # elif vi <= 0 or vi <= 0.0:
       #      #            def __str__(self):
       #      #                   return f'({str(self.real)} - {str(self.imag)}j)'
        
       #      # def __str__(self):
       #      #     return f'({str(self.real)} + {str(self.imag)}j)'
        
       #      # def __strs__(self):
       #      #     return f'({str(self.real)} - {str(self.imag)}j)'
        
       #      def sqrt(self):
       #          result = ComplexDecimal(self)
       #          if self.imag:
       #              raise NotImplementedError
       #          elif self.real > 0:
       #              result.real = self.real.sqrt()
       #              return result
       #          else:
       #              result.imag = (-self.real).sqrt()
       #              result.real = dm(0, self.prec)
       #              return result
             
       #      def real(self):
       #          result = ComplexDecimal(self)
       #          result.real = self.real
       #          return result
         
       #      def imag(self):
       #          result = ComplexDecimal(self)
       #          result.imag = self.imag
       #          return result
              
       def mmake(self, v):
           mp.dps = int(self.pr)
           try:
               mv = mp.mpf(v)
           except:
               mv = mp.mpc(v.real, v.imag)
           return mv
       
       def mmakec(self, vr, vi):
           mp.dps = int(self.pr)
           mv = mp.mpc(vr, vi)
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
       
       def dotd(self, v1, v2):
              vv = self.mmake(sum([self.mmake(x)*self.mmake(y) for x, y in zip(v1, v2)]))
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
                         uu = self.mmake(Um[i][n])
                         ulist.append(uu)
                   return ulist
            
           def lvalsd(Lm, k, n):
                   llist = []
                   lu = Lm[k]
                   lul = lu[0:k]
                   for i in range(len(lul)):
                            val_ij = self.mmake(lul[i])
                            llist.append(val_ij)
                   return lul
              
           for k in range(N):
               v1 = self.mmake(1.0)
               L[k][k] = v1
               # print(type(self.Ax[k][k]))
               try:
                   v2 = self.mmake((self.mmake(self.Ax[k][k]) - self.mmake(self.dotd(lvalsd(L, k, k), uvalsd(U, k, k)))) / self.mmake(L[k][k]))
               except:
                   # v2r = self.mmake((self.mmake(self.Ax[k][k].real) - self.mmake(self.dotd(lvalsd(L, k, k), uvalsd(U, k, k)))) / self.mmake(L[k][k].real))
                    # print(v2r)
                   v2b = self.dotd(lvalsd(L, k, k), uvalsd(U, k, k))
                   v2c = self.mmakec(self.Ax[k][k].real, self.Ax[k][k].imag) - self.mmakec(v2b.real, v2b.imag)
                   v2 = self.mmakec(v2c.real, v2c.imag) / self.mmakec(L[k][k].real, L[k][k].imag)
                   # v2 = DCD(DCS(self.Ax[k][k], self.dotd(lvalsd(L, k, k), uvalsd(U, k, k)), self.pr), L[k][k], self.pr)
                   # print(v1i)
                   # v2i = self.mmake((self.mmake(self.Ax[k][k].imag) - self.mmake(self.dotd(lvalsd(L, k, k), uvalsd(U, k, k)))) / self.mmake(L[k][k].imag))
               U[k][k] = v2
               for j in range(k+1, N):
                     val_i = self.mmake(self.mmake((self.mmake(self.Ax[k][j]) - self.mmake(self.dotd(lvalsd(L, k, k), uvalsd(U, j, j)))) / self.mmake(L[k][k])))
                     U[k][j] = val_i
               for i in range(k+1, N):
                     val_ib = self.mmake(self.mmake((self.mmake(self.Ax[i][k]) - self.mmake(self.dotd(lvalsd(L, i, i), uvalsd(U, k, k)))) / self.mmake(U[k][k])))
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
               v2 = (self.Ax[k][k] - float(self.dotd(lvalsd(L, k, k), uvalsd(U, k, k)))) / L[k][k]
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
                     val_i = self.mmake(self.mmake(y[i])-self.mmake(self.mmake(L[i][j])*self.mmake(y[j])))
                     y[i]= v3
               y[i] = self.mmake(self.mmake(y[i])/self.mmake(L[i][i]))
       
           return y
       
       def backward_subd(self, U, y):
           """Given a lower triangular matrix U and right-side vector y,
           compute the solution vector x solving Ux = y."""
            
           x = [self.mmake(0.0) for ix in y]
           # Us = len(U[0])
           # US = len(U[1])
           # for i in range(Us-1, -1, -1): 
           #     for j in range(i+1, US):
           #         y[i] -= U[i][j]*y[j]
           #     x[i] = y[i]/U[i][i]
           # return x
           for i in range(len(x)-1, 0, -1):
               val_i = self.mmake((self.mmake(y[i-1]) - self.mmake(self.dotd(U[i-1][i:], x[i:]))) / self.mmake(U[i-1][i-1]))
               x[i-1] = self.mmake(val_i)
              
           return x
       
       def forward_sub(self, L, b):
           """Given a lower triangular matrix L and right-side vector b,
           compute the solution vector y solving Ly = b."""
       
           y = []
           for i in range(len(b)):
               y.append(b[i])
               for j in range(i):
                     v1 = L[i][j]*y[j]
                     v2 = y[i] 
                     v3 = v2 - v1
                     val_i = self.mmake(self.mmake(y[i])-self.mmake(self.mmake(L[i][j])*self.mmake(y[j])))
                     y[i]= val_i #v3
               try:
                   try:
                       y[i] = self.mmake(y[i][0])/self.mmake(L[i][i][0])
                   except:
                       try:
                           y[i] = self.mmake(y[i][0])/self.mmake(L[i][i])
                       except:
                           pass
               except:
                   y[i] = self.mmake(y[i])/self.mmake(L[i][i])
       
           return y
       
       def backward_sub(self, U, y):
           """Given a lower triangular matrix U and right-side vector y,
           compute the solution vector x solving Ux = y."""
            
           x = [self.mmake(0.0) for ix in y]
           XX = [self.mmakec(0.0, 0.0) for ix in y]
           # Us = len(U[0])
           # US = len(U[1])
           # for i in range(Us-1, -1, -1): 
           #     for j in range(i+1, US):
           #         y[i] -= U[i][j]*y[j]
           #     x[i] = y[i]/U[i][i]
           # return x
           for i in range(len(x), 0, -1):
               try:
                        val_i = (self.mmakec(y[i-1].real, y[i-1].imag) - np.dot(U[i-1][i:], XX[i:])) / self.mmakec(U[i-1][i-1].real, U[i-1][i-1].imag)
               except:
                        val_i = self.mmake(y[i-1] - np.dot(U[i-1][i:], x[i:])) / self.mmake(U[i-1][i-1][0])
               x[i-1] = val_i
              
           return x
       
       def lu_solved(self, L, U, b):
           # Step 1: Solve Uy = b using forward substitution
           # Step 2: Solve Lx = y using backward substitution
           y = self.forward_sub(L, b)
           x = self.backward_sub(U, y)
           # yv = []
           # for i in range(len(y)):
           #    val_yi = y[i]
           #    yv.append(val_yi)
           # x = self.backward_sub(U, y)
           return x
       
       def linear_solved(self, bd):
           Ld, Ud = self.LU_decompositionf()
           x = self.lu_solved(Ld, Ud, bd)
           return x
       
       def decfunc(self):
              rows = len(self.Ax)
              cols = len(self.Ax[0])
              AD = [[self.mmake(ij) for ij in self.Ax[ii]] for ii in range(len(self.Ax))]
              for i in range(rows):
                     for j in range(cols):
                            val_ij = self.mmake(self.Ax[i][j])
                            AD[i][j] = val_ij
              return AD
        
       def decfuncb(self, Ax):
              rows = len(Ax)
              cols = len(Ax[0])
              AD = [[self.mmake(ij) for ij in Ax[ii]] for ii in range(len(Ax))]
              for i in range(rows):
                     for j in range(cols):
                            val_ij = self.mmake(Ax[i][j])
                            AD[i][j] = val_ij
              return AD
       
       def decfuncl(self, Axl):
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
       
       # def matrixsubc(self, A, B):
       #               Z = []
       #               for i in range(len(A)):
       #                   row = []
       #                   for j in range(len(A[i])):
       #                       row.append(DCS(A[i][j], B[i][j], self.prec)) 
       #                   Z.append(row)
                         
       #               return Z
       
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
           jj = 0                                                                                                                                                                                             
           for j in range(m):
                  
              rowv = max(range(j, m), key=lambda i: abs(M[i][j]))
              rv = max(r)
              if j == MM:
                     # print(rowv)
                     rowv = rowv + 1
                      
              else:
                     pass
              if rowv > rv:
                     r.append(rowv)
                     row = max(r)
              else:
                     row = max(r)
              # print(row)
              # if j != max(r):
                   # Swap the rows                                                                                                                                                                                                                            
              id_mat[j], id_mat[max(r)] = id_mat[max(r)], id_mat[j]
              ID.append(id_mat[j])
              
              # ID[j-jj] = id_mat[max(r)]
           # print(ID)
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
           # print(P)
           # print(PA)
       
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
           # print(P)
           # print(PA)
       
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
               print([round(x, self.pr) for x in row])
    
       def printmf(self, Mm):
           """
           Print a matrix one row at a time
               :param M: The matrix to be printed
           """
           for row in Mm:
               print([round(self.mmake(x), self.pr) for x in row])
               
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
              # yv = []
              # for i in range(len(y)):
              #    val_yi = y[i]
              #    yv.append(val_yi)
              
              
              return x
              
       
       def factor_solve(self, B):
              LU, PIV = self.LU_factor()
              # print(PIV)
              A = self.Ax
              AA = []
              for i in PIV:
                     AA.append(A[i])
              # AA = A[PIV]
              x = self.solve(AA, B)
              return x
       
       def factor_solvec(self, B):
              LU, PIV = self.LU_factorc()
              # print(PIV)
              A = self.Ax
              AA = []
              for i in PIV:
                     AA.append(A[i])
              # AA = A[PIV]
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
           # print(Us, US)
           # US = len(U[1])
           for i in range(Us-1, 0, -1): 
               for j in range(i+1, US):
                   y[i] -= self.mmake(U[i][j])*self.mmake(y[j])
               y[i] = self.mmake(y[i])/self.mmake(U[i][i])
           return y
    
       def solves(self, LU, B):
              y = self.ufsub(LU, B)
              x = self.bsub(LU, y)
              # yv = []
              # for i in range(len(y)):
              #    val_yi = self.mmake(y[i])
              #    yv.append(val_yi)
              # x = self.bsub(LU, yv)
              
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
