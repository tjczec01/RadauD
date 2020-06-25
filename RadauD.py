from __future__ import division, print_function, absolute_import
import butchertableau as bt
import decimalfunctions as df
import numpy as np
from scipy.sparse import csc_matrix, issparse, eye
from scipy.optimize._numdiff import group_columns
from scipy.optimize import OptimizeResult
from common import validate_max_step, validate_tol, select_initial_step, norm, num_jac, EPS, warn_extraneous, validate_first_step, OdeSolution
from base import OdeSolver, DenseOutput
import decimal
import inspect
import sys
import os

plist = list(sys.path)
ppath = os.path.dirname(os.path.abspath(__file__))
pth = '{}'.format(ppath)
if pth in plist:
       pass
       
else:
       sys.path.append(r'{}'.format(ppath))
       
__name__ = "RadauD"


pth = r'{}'.format(ppath)

__all__ = ["RadauD", "dm", "Sf", "flatten", "solve_collocation_system",  "DCA", "DCS", "DCM", "DCD", "predict_factor", "RadauDenseOutput", "ComplexDecimal", "normd", "CP", "Tile", "Tfunc"]

def dm(v, pr):
       
       De = decimal.Decimal
       val_i = float(v)
       try:
           aa = '{}:.{}f{}'.format('{', int(pr),'}')
       except:
           aa = '{}:.{}f{}'.format('{', int(pr[0]),'}')
       aa1 = '{}'.format(aa)
       aa2 = str(aa1)
       aa3 = str(aa2.format(val_i))
       aa4 = De(aa3)
       return aa4
   
dm.__name__ = "dm"
flatten = lambda l: [item for sublist in l for item in sublist] 

def Sf(o):
    __name__ = "Sf"
    return int((o + 1)/2)

S6 = 6 ** 0.5
Eb = [-13 - 7 * S6, -13 + 7 * S6, -1] 
# E = [dm(ei, PREC[0])/dm(3.0, PREC[0]) for ei in Eb]
# Ef = df.DF(E, PREC[0]).decfuncl(E)
order = 5
# X = bt.butcher(order, PREC[0])
# S = X.stage()
S = Sf(order)
# A, B, C = X.radau() 
# Ainv = X.inv(A)        
# T, TI = X.Tmat(Ainv)  
# P = X.P(C)
# D = df.DF(A, PREC[0])
# EIGS = X.eigs(Ainv)
# print(EIGS)
MU_REAL = 3.637834252744502
MU_COMPLEX = 2.6810828736277488 - 3.050430199247423j

# print(T)
# print(TI)
# print(P)
# print(A)
# print(Ainv)
# print(B)
# print(C)

A = [[0.19681547722366058, -0.06553542585019845, 0.023770974348220134], 
      [0.39442431473908746, 0.29207341166522777, -0.04154875212599763], 
      [0.376403062700466, 0.512485826188424, 0.111111111111110]]

AI = [[3.22474487139158, 1.16784008469041, -0.253197264742181], 
      [-3.56784008469039, 0.775255128608406, 1.05319726474218], 
      [5.53197264742194, -7.53197264742187, 5.00000000000001]]

B = [0.376403062700466, 0.512485826188424, 0.111111111111110]

C = [0.15505102572168228, 0.6449489742783174, 1.0]

T = [[0.0944387624889745, -0.141255295020953, 0.0300291941051473], 
      [0.250213122965334, 0.204129352293798, -0.38294211275726], 
      [1.0, 1.0, 0.0]]

# TD = df.DF(T, PREC[0]).decfunc()

TI = [[4.178718591551935, 0.32768282076106514, 0.5233764454994487], 
      [-4.178718591551935, -0.3276828207610649, 0.476623554500551], 
      [0.5028726349458223, -2.571926949855616, 0.5960392048282263]]

# TID = df.DF(TI, PREC[0]).decfunc()
# print(TID)

P = [[10.048809399827414, -25.62959144707665, 15.580782047249254], 
      [-1.38214273316075, 10.29625811374331, -8.914115380582556], 
      [0.3333333333333328, -2.6666666666666616, 3.333333333333328]]

# E = [dm(ei, PREC[0])/dm(3.0, PREC[0]) for ei in Eb]
# Ef = df.DF(E, PREC[0]).decfuncl(E)
# TD = df.DF(T, PREC[0]).decfunc()
# TID = df.DF(TI, PREC[0]).decfunc()
# PD = df.DF(P, PREC[0]).decfunc()

NEWTON_MAXITER = 6  # Maximum number of Newton iterations.
MIN_FACTOR = 0.2  # Minimum allowed decrease in a step size.
MAX_FACTOR = 10  # Maximum allowed increase in a step size.



def Tfunc(TI):
    tl = len(TI)
    TI_REAL = TI[0]
    TC = []
    for i in range(1, tl, 2):
       tc = [TI[i][ii] + 1j * TI[i+1][ii] for ii in range(tl)]
       TC.append(tc)
       
    return TI_REAL, TC

# TI_REAL, TI_COMPLEX = Tfunc(TI)

TI_REAL = TI[0]
TI_COMPLEX = [TI[1][i] + 1j * TI[2][i] for i in range(len(T[0]))]

def CP(L):
       total = []
       current = 1
       for p in flatten(L):
              current *= p
              total.append(current)
       final = total[-1]
       return final

def Tile(l, x, y=None):
       til = []
       tf = []
       for i in range(x):
              for li in l:
                     til.append(li)
       if y is None or y == 0 or y == 0.0:
             return til 
       else:
              for yi in range(y):
                     tf.append(til)
              return tf

class ComplexDecimal(object):
    def __init__(self, value, prec):
        self.real = dm(value.real, prec)
        self.imag = dm(value.imag, prec)
        self.vi = dm(value.imag, prec)
        self.prec = prec
        if self.vi >= 0 or self.vi >= 0.0:
               self.sign = '+'
        elif self.vi <= 0 or self.vi <= 0.0:
               self.sign = '-'
    def __add__(self, other):
        result = ComplexDecimal(self)
        result.real += dm(other.real, self.prec)
        if other.imag >= 0 or other.imag >= 0.0:
               result.imag += dm(other.imag, self.prec)
        elif other.imag <= 0 or other.imag <= 0.0:
               result.imag -= dm(other.imag, self.prec)
        return result

    __radd__ = __add__
    def __str__(self):
              return f'({self.real} {self.sign} {abs(self.imag)}j)'

    def sqrt(self):
        result = ComplexDecimal(self)
        if self.imag:
            raise NotImplementedError
        elif self.real > 0:
            result.real = self.real.sqrt()
            return result
        else:
            result.imag = (-self.real).sqrt()
            result.real = dm(0, self.prec)
            return result
     
    def real(self):
        result = ComplexDecimal(self)
        result.real = self.real
        return result
 
    def imag(self):
        result = ComplexDecimal(self)
        result.imag = self.imag
        return result   

def DCM(MUD, MUD2, prec):
       __name__ = "DCM"
       M1R = ComplexDecimal(MUD, prec).real
       M1C = ComplexDecimal(MUD, prec).imag
       M2R = ComplexDecimal(MUD2, prec).real
       M2C = ComplexDecimal(MUD2, prec).imag
       vv = M1R*M2R
       vv2 = M1C*M2R
       vv3 = M1R*M2C
       vv4 = M1C*M2C
       vvi = vv - vv4
       vvj = vv3 + vv2
       return ComplexDecimal(complex(vvi, vvj), prec)

def DCA(MUD, MUD2, prec):
       __name__ = "DCA"
       M1R = ComplexDecimal(MUD, prec).real
       M1C = ComplexDecimal(MUD, prec).imag
       M2R = ComplexDecimal(MUD2, prec).real
       M2C = ComplexDecimal(MUD2, prec).imag
       vv = M1R + M2R
       vv2 = M1C + M2C
       return ComplexDecimal(complex(vv, vv2), prec)

def DCS(MUD, MUD2, prec):
       __name__ = "DCS"
       M1R = ComplexDecimal(MUD, prec).real
       M1C = ComplexDecimal(MUD, prec).imag
       M2R = ComplexDecimal(MUD2, prec).real
       M2C = ComplexDecimal(MUD2, prec).imag
       vv = M1R - M2R
       vv2 = M1C - M2C
       return ComplexDecimal(complex(vv, vv2), prec)

def DCD(MUD, MUD2, prec):
       __name__ = "DCD"
       M3R = ComplexDecimal(MUD2, prec).real
       M3C = -ComplexDecimal(MUD2, prec).imag
       bot = DCM(MUD2, complex(M3R, M3C), prec).real
       top = DCM(MUD, complex(M3R, M3C), prec)
       TR = top.real / bot
       TC = top.imag / bot
       return ComplexDecimal(complex(TR, TC), prec)

def normd(A, prec):
       __name__ = "normd"
       rows = len(A)
       cols = len(A[0])
       vals = []
       for i in range(rows):
              
              for j in range(cols):
                     vi = dm(abs(A[i][j]), prec)**dm(2, prec)
                     vals.append(vi)
                     
       vf = dm(sum(vals), prec)**dm(1/2, prec)
       return vf

# v1 = ComplexDecimal(TI_COMPLEX[0], PREC[0]) 
# v2 = ComplexDecimal(TI_COMPLEX[1], PREC[0]) 
# v3 = ComplexDecimal(TI_COMPLEX[2], PREC[0])                     
# TI_COMPLEXb = [v1, v2, v3]

# MU_REALd = dm(3.637834252744502, PREC[0])
# MU_COMPLEXd = ComplexDecimal(2.6810828736277488 - 3.050430199247423j, PREC[0])

def solve_collocation_system(fun, t, y, h, Z0, scale, tol,
                             LU_real, LU_complex, solve_lu, prec=9):
    """Solve the collocation system.
    Parameters
    ----------
    fun : callable
        Right-hand side of the system.
    t : float
        Current time.
    y : ndarray, shape (n,)
        Current state.
    h : float
        Step to try.
    Z0 : ndarray, shape (3, n)
        Initial guess for the solution. It determines new values of `y` at
        ``t + h * C`` as ``y + Z0``, where ``C`` is the Radau method constants.
    scale : float
        Problem tolerance scale, i.e. ``rtol * abs(y) + atol``.
    tol : float
        Tolerance to which solve the system. This value is compared with
        the normalized by `scale` error.
    LU_real, LU_complex
        LU decompositions of the system Jacobians.
    solve_lu : callable
        Callable which solves a linear system given a LU decomposition. The
        signature is ``solve_lu(LU, b)``.
    Returns
    -------
    converged : bool
        Whether iterations converged.
    n_iter : int
        Number of completed iterations.
    Z : ndarray, shape (3, n)
        Found solution.
    rate : float
        The rate of convergence.
    """
    TD = df.DF(T, prec).decfunc()
    TID = df.DF(TI, prec).decfunc()
    E = [dm(ei, prec)/dm(3.0, prec) for ei in Eb]
    Ef = df.DF(E, prec).decfuncl(E)
    n = len([y])
    M_real = MU_REAL / h
    M_complex = MU_COMPLEX / float(h) #DCD(MU_COMPLEX, h, prec) 
    ch = [dm(h*i, prec) for i in range(len(C))]
    ZP = np.zeros((S, n), dtype=float).tolist()
    for i in range(S):
            for j in range(n):
                   try:
                       ZP[i][j] = fun(dm(t + float(h)*float(C[i]), prec), dm(y, prec) + (dm(Z0[i][j], prec))).tolist()
                   except:
                       ZP[i][j] = [dm(fun(dm(t + float(h)*float(C[i]), prec), dm(y[0], prec) + (dm(Z0[i][j], prec))).tolist()[0], prec)]
    Z0P = np.array(flatten(ZP))
    Zd = df.DF(Z0P, prec).decfunc()
    try:
           W = np.array(np.array(TID).dot(Zd))
    except:
           W = np.array(np.array(TID).dot(Zd)) #, dtype=float
    Z = Z0
    F = np.empty((S, n))
    ch = [dm(h*i, prec) for i in range(len(C))]
    dW_norm_old = None
    dW = np.empty_like(W)
    dW2 = np.empty_like(W)
    converged = False
    rate = None
    for k in range(NEWTON_MAXITER):
        for i in range(S):
                   try:
                       F[i] = dm(fun(dm(t, prec) + ch[i], dm(y, prec) + dm(Z[i][0], prec)), prec)
                   except:
                       try:
                           F[i] = dm(fun(dm(t, prec) + ch[i], dm(y, prec) + dm(Z[i][0], prec))[0], prec)
                       except:
                           F[i] = dm(fun(dm(t, prec) + ch[i], dm(y[0], prec) + dm(Z[i][0], prec))[0], prec)

        if not np.all(np.isfinite(F)):
            break
        try:
            try:
                f_real = [dm(F.T.dot(TI_REAL), prec) - dm(M_real[0], prec) * dm(W[0][0], prec)]
            except:
                try:
                    f_real = [dm(F.T.dot(TI_REAL)[0], prec) - dm(M_real, prec) * dm(W[0][0], prec)]
                except:
                    f_real = [dm(F.T.dot(TI_REAL)[0], prec) - dm(M_real[0], prec) * dm(W[0], prec)]
        except:
            f_real = [dm(F.T.dot(TI_REAL), prec) - dm(M_real[0], prec) * dm(W[0][0], prec)]
        # fc1 = F.T.dot(TI_COMPLEX) - M_complex * TI_COMPLEX
        # fc2 = [M_complex * TI_COMPLEX[0][i] for i in range(len(TI_COMPLEX))]
        # fc3 = [i - j  for i, j in zip(fc1, fc2)]
        # f_complex = [F.dot(TI_COMPLEX)[0].tolist()[i] - [M_complex * TI_COMPLEX[0][i] for i in range(len(TI_COMPLEX))][i] for i in range(len(TI_COMPLEX))][0] #F.dot(TI_COMPLEX)[0] - [M_complex * TI_COMPLEX[0][i] for i in TI_COMPLEX]
        f_complex = F.T.dot(TI_COMPLEX) - M_complex *  complex(W[1] ,  W[2])
        dW_real = solve_lu(LU_real[0], LU_real[1], [f_real])
        dW_complex = solve_lu(LU_complex[0], LU_complex[1], [f_complex])
        dW[0] = dm(dW_real[0].real, prec)
        dW[1] = dm(dW_complex[0].real, prec)
        dW[2] = dm(dW_complex[0].imag, prec)
        dW2[0] = dW_real[0].real
        dW2[1] = dW_complex[0].real
        dW2[2] = dW_complex[0].imag
        dW_norm = normd(dW/scale, prec)
        if dW_norm_old is not None:
            rate = dW_norm / dW_norm_old

        if (rate is not None and (rate >= 1 or
                rate ** (NEWTON_MAXITER - k) / (1 - rate) * dW_norm > tol)):
            break
        W += dW
        Z = np.array(TD).dot(W)

        if (dW_norm == 0 or
                rate is not None and rate / (1 - rate) * dW_norm < tol):
            converged = True
            break

        dW_norm_old = dW_norm

    return converged, k + 1, Z, rate


def predict_factor(h_abs, h_abs_old, error_norm, error_norm_old, prec):
    """Predict by which factor to increase/decrease the step size.
    The algorithm is described in [1]_.
    Parameters
    ----------
    h_abs, h_abs_old : float
        Current and previous values of the step size, `h_abs_old` can be None
        (see Notes).
    error_norm, error_norm_old : float
        Current and previous values of the error norm, `error_norm_old` can
        be None (see Notes).
    Returns
    -------
    factor : float
        Predicted factor.
    Notes
    -----
    If `h_abs_old` and `error_norm_old` are both not None then a two-step
    algorithm is used, otherwise a one-step algorithm is used.
    References
    ----------
    .. [1] E. Hairer, S. P. Norsett G. Wanner, "Solving Ordinary Differential
           Equations II: Stiff and Differential-Algebraic Problems", Sec. IV.8.
    """
    if error_norm_old is None or h_abs_old is None or error_norm == 0:
        multiplier = 1
    else:
        multiplier = float(h_abs) / float(h_abs_old) * (float(error_norm_old) / float(error_norm)) ** 0.25

    with np.errstate(divide='ignore'):
        factor = min(1, multiplier) * error_norm ** -0.25

    return factor


class RadauD(OdeSolver):
    """Implicit Runge-Kutta method of Radau IIA family of order 5.
    The implementation follows [1]_. The error is controlled with a
    third-order accurate embedded formula. A cubic polynomial which satisfies
    the collocation conditions is used for the dense output.
    Parameters
    ----------
    fun : callable
        Right-hand side of the system. The calling signature is ``fun(t, y)``.
        Here ``t`` is a scalar, and there are two options for the ndarray ``y``:
        It can either have shape (n,); then ``fun`` must return array_like with
        shape (n,). Alternatively it can have shape (n, k); then ``fun``
        must return an array_like with shape (n, k), i.e. each column
        corresponds to a single column in ``y``. The choice between the two
        options is determined by `vectorized` argument (see below). The
        vectorized implementation allows a faster approximation of the Jacobian
        by finite differences (required for this solver).
    t0 : float
        Initial time.
    y0 : array_like, shape (n,)
        Initial state.
    t_bound : float
        Boundary time - the integration won't continue beyond it. It also
        determines the direction of the integration.
    first_step : float or None, optional
        Initial step size. Default is ``None`` which means that the algorithm
        should choose.
    max_step : float, optional
        Maximum allowed step size. Default is np.inf, i.e. the step size is not
        bounded and determined solely by the solver.
    rtol, atol : float and array_like, optional
        Relative and absolute tolerances. The solver keeps the local error
        estimates less than ``atol + rtol * abs(y)``. Here `rtol` controls a
        relative accuracy (number of correct digits). But if a component of `y`
        is approximately below `atol`, the error only needs to fall within
        the same `atol` threshold, and the number of correct digits is not
        guaranteed. If components of y have different scales, it might be
        beneficial to set different `atol` values for different components by
        passing array_like with shape (n,) for `atol`. Default values are
        1e-3 for `rtol` and 1e-6 for `atol`.
    jac : {None, array_like, sparse_matrix, callable}, optional
        Jacobian matrix of the right-hand side of the system with respect to
        y, required by this method. The Jacobian matrix has shape (n, n) and
        its element (i, j) is equal to ``d f_i / d y_j``.
        There are three ways to define the Jacobian:
            * If array_like or sparse_matrix, the Jacobian is assumed to
              be constant.
            * If callable, the Jacobian is assumed to depend on both
              t and y; it will be called as ``jac(t, y)`` as necessary.
              For the 'Radau' and 'BDF' methods, the return value might be a
              sparse matrix.
            * If None (default), the Jacobian will be approximated by
              finite differences.
        It is generally recommended to provide the Jacobian rather than
        relying on a finite-difference approximation.
    jac_sparsity : {None, array_like, sparse matrix}, optional
        Defines a sparsity structure of the Jacobian matrix for a
        finite-difference approximation. Its shape must be (n, n). This argument
        is ignored if `jac` is not `None`. If the Jacobian has only few non-zero
        elements in *each* row, providing the sparsity structure will greatly
        speed up the computations [2]_. A zero entry means that a corresponding
        element in the Jacobian is always zero. If None (default), the Jacobian
        is assumed to be dense.
    vectorized : bool, optional
        Whether `fun` is implemented in a vectorized fashion. Default is False.
    Attributes
    ----------
    n : int
        Number of equations.
    status : string
        Current status of the solver: 'running', 'finished' or 'failed'.
    t_bound : float
        Boundary time.
    direction : float
        Integration direction: +1 or -1.
    t : float
        Current time.
    y : ndarray
        Current state.
    t_old : float
        Previous time. None if no steps were made yet.
    step_size : float
        Size of the last successful step. None if no steps were made yet.
    nfev : int
        Number of evaluations of the right-hand side.
    njev : int
        Number of evaluations of the Jacobian.
    nlu : int
        Number of LU decompositions.
    References
    ----------
    .. [1] E. Hairer, G. Wanner, "Solving Ordinary Differential Equations II:
           Stiff and Differential-Algebraic Problems", Sec. IV.8.
    .. [2] A. Curtis, M. J. D. Powell, and J. Reid, "On the estimation of
           sparse Jacobian matrices", Journal of the Institute of Mathematics
           and its Applications, 13, pp. 117-120, 1974.
    """
    def __init__(self, fun, t0, y0, t_bound, max_step=np.inf,
                 rtol=1e-3, atol=1e-6, jac=None, jac_sparsity=None,
                 vectorized=False, first_step=None, prec=25, **extraneous):
        warn_extraneous(extraneous)
        super(RadauD, self).__init__(fun, t0, y0, t_bound, vectorized)
        self.y_old = None
        self.prec = prec
        self.max_step = validate_max_step(max_step)
        self.rtol, self.atol = validate_tol(rtol, atol, self.n)
        self.f = self.fun(self.t, self.y)
        E = [dm(ei, self.prec)/dm(3.0, self.prec) for ei in Eb]
        Ef = df.DF(E, self.prec).decfuncl(E)
        TD = df.DF(T, self.prec).decfunc()
        TID = df.DF(TI, self.prec).decfunc()
        PD = df.DF(P, self.prec).decfunc()
        # Select initial step assuming the same order which is used to control
        # the error.
        if first_step is None:
            self.h_abs = select_initial_step(self.fun, self.t, self.y, self.f, self.direction, S, self.rtol, self.atol)
        else:
            self.h_abs = validate_first_step(first_step, t0, t_bound)
        self.h_abs_old = None
        self.error_norm_old = None

        self.newton_tol = max(10  * EPS / rtol, min(0.03, rtol ** 0.5))
        self.sol = None

        self.jac_factor = None
        self.jac, self.J = self._validate_jac(jac, jac_sparsity)
        if issparse(self.J):
            def lu(A):
                self.nlu += 1
                A1 = df.DF(A, self.prec)
                L, U = A1.LU_decompositiond()
                return L, U

            def solve_lu(L, U, b):
                AA = df.DF([L, U], self.prec).lu_solved(L, U, b)
                return AA

            I = eye(self.n, format='csc')
        else:
            def lu(A):
                self.nlu += 1
                A1 = df.DF(A, self.prec)
                L, U = A1.LU_decompositiond()
                return L, U

            def solve_lu(L, U, b):
                AA = df.DF([L, U], self.prec).lu_solved(L, U, b)
                return AA

            I = np.identity(self.n)

        self.lu = lu
        self.solve_lu = solve_lu
        self.I = I

        self.current_jac = True
        self.LU_real = None
        self.LU_complex = None
        self.Z = None

    def _validate_jac(self, jac, sparsity):
        t0 = self.t
        y0 = self.y

        if jac is None:
            if sparsity is not None:
                if issparse(sparsity):
                    sparsity = csc_matrix(sparsity)
                groups = group_columns(sparsity)
                sparsity = (sparsity, groups)

            def jac_wrapped(t, y, f):
                self.njev += 1
                J, self.jac_factor = num_jac(self.fun_vectorized, t, y, f,
                                             self.atol, self.jac_factor,
                                             sparsity)
                return J
            J = jac_wrapped(t0, y0, self.f)
        elif callable(jac):
            J = jac(t0, y0)
            self.njev = 1
            if issparse(J):
                J = csc_matrix(J)

                def jac_wrapped(t, y, _=None):
                    self.njev += 1
                    return csc_matrix(jac(t, y), dtype=float)

            else:
                J = np.asarray(J, dtype=float)

                def jac_wrapped(t, y, _=None):
                    self.njev += 1
                    return np.asarray(jac(t, y), dtype=float)

            if J.shape != (self.n, self.n):
                raise ValueError("`jac` is expected to have shape {}, but "
                                 "actually has {}."
                                 .format((self.n, self.n), J.shape))
        else:
            if issparse(jac):
                J = csc_matrix(jac)
            else:
                J = np.asarray(jac, dtype=float)

            if J.shape != (self.n, self.n):
                raise ValueError("`jac` is expected to have shape {}, but "
                                 "actually has {}."
                                 .format((self.n, self.n), J.shape))
            jac_wrapped = None

        return jac_wrapped, J

    def _step_impl(self):
        t = self.t
        y = self.y
        f = self.f
        E = [dm(ei, self.prec)/dm(3.0, self.prec) for ei in Eb]
        Ef = df.DF(E, self.prec).decfuncl(E)
        # TD = df.DF(T, self.prec).decfunc()
        # TID = df.DF(TI, self.prec).decfunc()
        # PD = df.DF(P, self.prec).decfunc()
        max_step = self.max_step
        atol = self.atol
        rtol = self.rtol

        min_step = 10 * np.abs(np.nextafter(float(t), self.direction * np.inf) - float(t))        
        if self.h_abs > max_step:
            h_abs = max_step
            h_abs_old = None
            error_norm_old = None
        elif self.h_abs < min_step:
            h_abs = float(min_step)
            h_abs_old = None
            error_norm_old = None
        else:
            h_abs = float(self.h_abs)
            h_abs_old = self.h_abs_old
            error_norm_old = self.error_norm_old

        J = self.J
        LU_real = self.LU_real
        LU_complex = self.LU_complex

        current_jac = self.current_jac
        jac = self.jac

        rejected = False
        step_accepted = False
        message = None
        while not step_accepted:
            if h_abs < min_step:
                return False, self.TOO_SMALL_STEP
            h = float(h_abs) * self.direction
            t_new = float(t) + float(h)

            if self.direction * (float(t_new) - float(self.t_bound)) > 0:
                t_new = self.t_bound

            h = float(t_new) - float(t)
            h_abs = np.abs(h)

            if self.sol is None:
                Z0 = np.zeros((S, y.shape[0]))
            else:
                CH = np.array([t + h*i for i in C])
                Z1 = self.sol(CH).T
                Z2 = [dm(i[0], self.prec) for i in Z1.tolist()]
                Y = [dm(i, self.prec) for i in [y]]
                Z0  = np.array([[Z2[i] - Y[j] for j in range(len(Y))] for i in range(len(Z2))])
            try:
                yv = y[0]
            except:
                yv = y
            scale = dm(dm(float(atol), self.prec) + dm(yv, self.prec) * dm(rtol, self.prec), self.prec)
            
            def matrix_subtraction(A, B):
                  """
                  Subtracts matrix B from matrix A and returns difference
                      :param A: The first matrix
                      :param B: The second matrix
                      :return: Matrix difference
                  """
                  # Section 1: Ensure dimensions are valid for matrix subtraction
                  rowsA = len(A)
                  colsA = len(A[0])
                  rowsB = len(B)
                  colsB = len(B[0])
                  if rowsA != rowsB or colsA != colsB:
                      raise ArithmeticError('Matrices are NOT the same size.')
              
                  # Section 2: Create a new matrix for the matrix difference
                  C = df.DF(A, self.prec).zeros_matrix(rowsA, colsB)
              
                  # Section 3: Perform element by element subtraction
                  for i in range(rowsA):
                      for j in range(colsB):
                          C[i][j] = A[i][j] - B[i][j]
              
                  return C
           
            def matrixsubc(A, B):
                     Z = []
                     for i in range(len(A)):
                         row = []
                         for j in range(len(A[i])):
                             VC = DCS(complex(A[i][j].real, A[i][j].imag), complex(B[i][j].real, B[i][j].imag), self.prec)
                             row.append(complex(VC.real, VC.imag)) 
                         Z.append(row)
                         
                     return Z
            converged = False
            while not converged:
                if LU_real is None or LU_complex is None:
                    # JJ = df.DF([0.0], self.prec)
                    # vr = MU_REAL / float(h)
                    # vc = MU_COMPLEX / complex(h, 0.0)
                    # IR = df.DF([0.0], self.prec).idenv(len(J), vr)
                    # IC = df.DF([0.0], self.prec).idenfv(len(J), complex(vc.real, vc.imag))
                    Lr, Ur = self.lu(MU_REAL / h * self.I - J)
                    LU_real = [Lr, Ur]
                    Lc, Uc = self.lu(MU_COMPLEX / h * self.I - J)
                    LU_complex = [Lc, Uc]

                converged, n_iter, Z, rate = solve_collocation_system(
                    self.fun, t, y, h, Z0, scale, self.newton_tol,
                    LU_real, LU_complex, self.solve_lu, self.prec)

                if not converged:
                    if current_jac:
                        break

                    J = self.jac(t, y, f)
                    current_jac = True
                    LU_real = None
                    LU_complex = None

            if not converged:
                h_abs *= 0.5
                LU_real = None
                LU_complex = None
                continue
            try:
                y_new = dm(y[0], self.prec) + dm(Z[-1][0], self.prec)
            except:
                y_new = dm(y, self.prec) + dm(Z[-1][0], self.prec)
            ZE = np.array(Z).T.dot(np.array(Ef)) / dm(h, self.prec)
            
            error = self.solve_lu(LU_real[0], LU_real[1], [dm(f[0], self.prec) + dm(ZE[0], self.prec)])
            scale = atol + float(np.maximum(np.abs(y), np.abs(y_new))) * rtol
            error_norm = norm(error[0].real / scale)
            safety = 0.9 * (2 * NEWTON_MAXITER + 1) / (2 * NEWTON_MAXITER + n_iter)

            if rejected and error_norm > 1:
                error = df.DF(LU_real, self.prec).solve(LU_real, self.fun(t, float(y[0]) + float(error[0].real)) + float(ZE))
                error_norm = norm(float(error[0]) / scale)

            if error_norm > 1:
                factor = predict_factor(h_abs, h_abs_old,
                                        error_norm, error_norm_old, self.prec)
                h_abs *= max(MIN_FACTOR, safety * factor)

                LU_real = None
                LU_complex = None
                rejected = True
            else:
                step_accepted = True

        recompute_jac = jac is not None and n_iter > 2 and rate > 1e-3

        factor = predict_factor(h_abs, h_abs_old, error_norm, error_norm_old, self.prec)
        factor = min(MAX_FACTOR, safety * factor)

        if not recompute_jac and factor < 1.2:
            factor = 1
        else:
            LU_real = None
            LU_complex = None

        f_new = self.fun(t_new, y_new)
        if recompute_jac:
            J = jac(t_new, y_new, f_new)
            current_jac = True
        elif jac is not None:
            current_jac = False

        self.h_abs_old = self.h_abs
        self.error_norm_old = error_norm

        self.h_abs = h_abs * factor

        self.y_old = y

        self.t = t_new
        self.y = y_new
        self.f = f_new

        self.Z = Z

        self.LU_real = LU_real
        self.LU_complex = LU_complex
        self.current_jac = current_jac
        self.J = J

        self.t_old = t
        self.sol = self._compute_dense_output()

        return step_accepted, message

    def _compute_dense_output(self):
        Qb = df.DF(self.Z, self.prec)
        Q = Qb.matrix_multiplyd(df.DF(self.Z, self.prec).TD(), df.DF(P, self.prec).decfunc(), self.prec)
        return RadauDenseOutput(self.t_old, self.t, self.y_old, Q, self.prec)

    def _dense_output_impl(self):
        return self.sol


class RadauDenseOutput(DenseOutput):
    def __init__(self, t_old, t, y_old, Q, prec):
        super(RadauDenseOutput, self).__init__(t_old, t)
        self.h = t - float(t_old)
        self.Q = Q
        self.order =  np.array(Q).shape[1] - 1
        self.y_old = y_old
        self.prec = prec

    def _call_impl(self, t):
        x = (np.array(t, dtype=float) - np.array(self.t_old, dtype=float)) / float(self.h)
        xd = [dm(xi, self.prec) for xi in x.tolist()]
        if t.ndim == 0:
            pb = np.tile(np.array(xd), self.order + 1, 0)
            p = np.cumprod(pb)
        else:
            pb = np.tile(np.array(xd), (self.order + 1, 1))
            p = np.cumprod(pb, axis=0)
        # Here we don't multiply by h, not a mistake.
        # try:
        #     pl = [[dm(pii, self.prec) for pii in range(len(p.tolist()[0]))] for pi in p.tolist()]
        # except:
        #     pl = dm(p[0], self.prec)
        y = np.dot(self.Q, p)
        yv = y.tolist()[0]
        if y.ndim == 2:
            try:
                try:
                    yvb = [dm(i, self.prec) + dm(self.y_old[0], self.prec) for i in yv]
                except:
                    try:
                        yvb = [dm(i[0], self.prec) + dm(self.y_old, self.prec) for i in yv]
                    except:
                        try:
                            yvb = [dm(i, self.prec) + dm(self.y_old, self.prec) for i in yv]
                        except:
                            yvb = [dm(i[0], self.prec) + dm(self.y_old[0], self.prec) for i in yv]
            except:
                yvb = [dm(i, self.prec) + dm(self.y_old, self.prec) for i in yv]
        else:
            y += self.y_old

        return np.array([yvb])
 
METHODS = {'RadauD' : RadauD}

MESSAGES = {0: "The solver successfully reached the end of the integration interval.",
            1: "A termination event occurred."}


class OdeResult(OptimizeResult):
    pass


def prepare_events(events):
    """Standardize event functions and extract is_terminal and direction."""
    if callable(events):
        events = (events,)

    if events is not None:
        is_terminal = np.empty(len(events), dtype=bool)
        direction = np.empty(len(events))
        for i, event in enumerate(events):
            try:
                is_terminal[i] = event.terminal
            except AttributeError:
                is_terminal[i] = False

            try:
                direction[i] = event.direction
            except AttributeError:
                direction[i] = 0
    else:
        is_terminal = None
        direction = None

    return events, is_terminal, direction


def solve_event_equation(event, sol, t_old, t):
    """Solve an equation corresponding to an ODE event.
    The equation is ``event(t, y(t)) = 0``, here ``y(t)`` is known from an
    ODE solver using some sort of interpolation. It is solved by
    `scipy.optimize.brentq` with xtol=atol=4*EPS.
    Parameters
    ----------
    event : callable
        Function ``event(t, y)``.
    sol : callable
        Function ``sol(t)`` which evaluates an ODE solution between `t_old`
        and  `t`.
    t_old, t : float
        Previous and new values of time. They will be used as a bracketing
        interval.
    Returns
    -------
    root : float
        Found solution.
    """
    from scipy.optimize import brentq
    return brentq(lambda t: event(t, sol(t)), t_old, t,
                  xtol=4 * EPS, rtol=4 * EPS)


def handle_events(sol, events, active_events, is_terminal, t_old, t):
    """Helper function to handle events.
    Parameters
    ----------
    sol : DenseOutput
        Function ``sol(t)`` which evaluates an ODE solution between `t_old`
        and  `t`.
    events : list of callables, length n_events
        Event functions with signatures ``event(t, y)``.
    active_events : ndarray
        Indices of events which occurred.
    is_terminal : ndarray, shape (n_events,)
        Which events are terminal.
    t_old, t : float
        Previous and new values of time.
    Returns
    -------
    root_indices : ndarray
        Indices of events which take zero between `t_old` and `t` and before
        a possible termination.
    roots : ndarray
        Values of t at which events occurred.
    terminate : bool
        Whether a terminal event occurred.
    """
    roots = [solve_event_equation(events[event_index], sol, t_old, t)
             for event_index in active_events]

    roots = np.asarray(roots)

    if np.any(is_terminal[active_events]):
        if t > t_old:
            order = np.argsort(roots)
        else:
            order = np.argsort(-roots)
        active_events = active_events[order]
        roots = roots[order]
        t = np.nonzero(is_terminal[active_events])[0][0]
        active_events = active_events[:t + 1]
        roots = roots[:t + 1]
        terminate = True
    else:
        terminate = False

    return active_events, roots, terminate


def find_active_events(g, g_new, direction):
    """Find which event occurred during an integration step.
    Parameters
    ----------
    g, g_new : array_like, shape (n_events,)
        Values of event functions at a current and next points.
    direction : ndarray, shape (n_events,)
        Event "direction" according to the definition in `solve_ivp`.
    Returns
    -------
    active_events : ndarray
        Indices of events which occurred during the step.
    """
    g, g_new = np.asarray(g), np.asarray(g_new)
    up = (g <= 0) & (g_new >= 0)
    down = (g >= 0) & (g_new <= 0)
    either = up | down
    mask = (up & (direction > 0) |
            down & (direction < 0) |
            either & (direction == 0))

    return np.nonzero(mask)[0]


def solve_ivpd(fun, t_span, y0, method='RK45', prec=25, t_eval=None, dense_output=False,
              events=None, vectorized=False, args=None, **options):
    """Solve an initial value problem for a system of ODEs.
    This function numerically integrates a system of ordinary differential
    equations given an initial value::
        dy / dt = f(t, y)
        y(t0) = y0
    Here t is a one-dimensional independent variable (time), y(t) is an
    n-dimensional vector-valued function (state), and an n-dimensional
    vector-valued function f(t, y) determines the differential equations.
    The goal is to find y(t) approximately satisfying the differential
    equations, given an initial value y(t0)=y0.
    Some of the solvers support integration in the complex domain, but note
    that for stiff ODE solvers, the right-hand side must be
    complex-differentiable (satisfy Cauchy-Riemann equations [11]_).
    To solve a problem in the complex domain, pass y0 with a complex data type.
    Another option always available is to rewrite your problem for real and
    imaginary parts separately.
    Parameters
    ----------
    fun : callable
        Right-hand side of the system. The calling signature is ``fun(t, y)``.
        Here `t` is a scalar, and there are two options for the ndarray `y`:
        It can either have shape (n,); then `fun` must return array_like with
        shape (n,). Alternatively it can have shape (n, k); then `fun`
        must return an array_like with shape (n, k), i.e. each column
        corresponds to a single column in `y`. The choice between the two
        options is determined by `vectorized` argument (see below). The
        vectorized implementation allows a faster approximation of the Jacobian
        by finite differences (required for stiff solvers).
    t_span : 2-tuple of floats
        Interval of integration (t0, tf). The solver starts with t=t0 and
        integrates until it reaches t=tf.
    y0 : array_like, shape (n,)
        Initial state. For problems in the complex domain, pass `y0` with a
        complex data type (even if the initial value is purely real).
    method : string or `OdeSolver`, optional
        Integration method to use:
            * 'RK45' (default): Explicit Runge-Kutta method of order 5(4) [1]_.
              The error is controlled assuming accuracy of the fourth-order
              method, but steps are taken using the fifth-order accurate
              formula (local extrapolation is done). A quartic interpolation
              polynomial is used for the dense output [2]_. Can be applied in
              the complex domain.
            * 'RK23': Explicit Runge-Kutta method of order 3(2) [3]_. The error
              is controlled assuming accuracy of the second-order method, but
              steps are taken using the third-order accurate formula (local
              extrapolation is done). A cubic Hermite polynomial is used for the
              dense output. Can be applied in the complex domain.
            * 'DOP853': Explicit Runge-Kutta method of order 8 [13]_.
              Python implementation of the "DOP853" algorithm originally
              written in Fortran [14]_. A 7-th order interpolation polynomial
              accurate to 7-th order is used for the dense output.
              Can be applied in the complex domain.
            * 'Radau': Implicit Runge-Kutta method of the Radau IIA family of
              order 5 [4]_. The error is controlled with a third-order accurate
              embedded formula. A cubic polynomial which satisfies the
              collocation conditions is used for the dense output.
            * 'BDF': Implicit multi-step variable-order (1 to 5) method based
              on a backward differentiation formula for the derivative
              approximation [5]_. The implementation follows the one described
              in [6]_. A quasi-constant step scheme is used and accuracy is
              enhanced using the NDF modification. Can be applied in the
              complex domain.
            * 'LSODA': Adams/BDF method with automatic stiffness detection and
              switching [7]_, [8]_. This is a wrapper of the Fortran solver
              from ODEPACK.
        Explicit Runge-Kutta methods ('RK23', 'RK45', 'DOP853') should be used
        for non-stiff problems and implicit methods ('Radau', 'BDF') for
        stiff problems [9]_. Among Runge-Kutta methods, 'DOP853' is recommended
        for solving with high precision (low values of `rtol` and `atol`).
        If not sure, first try to run 'RK45'. If it makes unusually many
        iterations, diverges, or fails, your problem is likely to be stiff and
        you should use 'Radau' or 'BDF'. 'LSODA' can also be a good universal
        choice, but it might be somewhat less convenient to work with as it
        wraps old Fortran code.
        You can also pass an arbitrary class derived from `OdeSolver` which
        implements the solver.
    t_eval : array_like or None, optional
        Times at which to store the computed solution, must be sorted and lie
        within `t_span`. If None (default), use points selected by the solver.
    dense_output : bool, optional
        Whether to compute a continuous solution. Default is False.
    events : callable, or list of callables, optional
        Events to track. If None (default), no events will be tracked.
        Each event occurs at the zeros of a continuous function of time and
        state. Each function must have the signature ``event(t, y)`` and return
        a float. The solver will find an accurate value of `t` at which
        ``event(t, y(t)) = 0`` using a root-finding algorithm. By default, all
        zeros will be found. The solver looks for a sign change over each step,
        so if multiple zero crossings occur within one step, events may be
        missed. Additionally each `event` function might have the following
        attributes:
            terminal: bool, optional
                Whether to terminate integration if this event occurs.
                Implicitly False if not assigned.
            direction: float, optional
                Direction of a zero crossing. If `direction` is positive,
                `event` will only trigger when going from negative to positive,
                and vice versa if `direction` is negative. If 0, then either
                direction will trigger event. Implicitly 0 if not assigned.
        You can assign attributes like ``event.terminal = True`` to any
        function in Python. 
    vectorized : bool, optional
        Whether `fun` is implemented in a vectorized fashion. Default is False.
    args : tuple, optional
        Additional arguments to pass to the user-defined functions.  If given,
        the additional arguments are passed to all user-defined functions.
        So if, for example, `fun` has the signature ``fun(t, y, a, b, c)``,
        then `jac` (if given) and any event functions must have the same
        signature, and `args` must be a tuple of length 3.
    options
        Options passed to a chosen solver. All options available for already
        implemented solvers are listed below.
    first_step : float or None, optional
        Initial step size. Default is `None` which means that the algorithm
        should choose.
    max_step : float, optional
        Maximum allowed step size. Default is np.inf, i.e. the step size is not
        bounded and determined solely by the solver.
    rtol, atol : float or array_like, optional
        Relative and absolute tolerances. The solver keeps the local error
        estimates less than ``atol + rtol * abs(y)``. Here `rtol` controls a
        relative accuracy (number of correct digits). But if a component of `y`
        is approximately below `atol`, the error only needs to fall within
        the same `atol` threshold, and the number of correct digits is not
        guaranteed. If components of y have different scales, it might be
        beneficial to set different `atol` values for different components by
        passing array_like with shape (n,) for `atol`. Default values are
        1e-3 for `rtol` and 1e-6 for `atol`.
    jac : array_like, sparse_matrix, callable or None, optional
        Jacobian matrix of the right-hand side of the system with respect
        to y, required by the 'Radau', 'BDF' and 'LSODA' method. The
        Jacobian matrix has shape (n, n) and its element (i, j) is equal to
        ``d f_i / d y_j``.  There are three ways to define the Jacobian:
            * If array_like or sparse_matrix, the Jacobian is assumed to
              be constant. Not supported by 'LSODA'.
            * If callable, the Jacobian is assumed to depend on both
              t and y; it will be called as ``jac(t, y)`` as necessary.
              For 'Radau' and 'BDF' methods, the return value might be a
              sparse matrix.
            * If None (default), the Jacobian will be approximated by
              finite differences.
        It is generally recommended to provide the Jacobian rather than
        relying on a finite-difference approximation.
    jac_sparsity : array_like, sparse matrix or None, optional
        Defines a sparsity structure of the Jacobian matrix for a finite-
        difference approximation. Its shape must be (n, n). This argument
        is ignored if `jac` is not `None`. If the Jacobian has only few
        non-zero elements in *each* row, providing the sparsity structure
        will greatly speed up the computations [10]_. A zero entry means that
        a corresponding element in the Jacobian is always zero. If None
        (default), the Jacobian is assumed to be dense.
        Not supported by 'LSODA', see `lband` and `uband` instead.
    lband, uband : int or None, optional
        Parameters defining the bandwidth of the Jacobian for the 'LSODA'
        method, i.e., ``jac[i, j] != 0 only for i - lband <= j <= i + uband``.
        Default is None. Setting these requires your jac routine to return the
        Jacobian in the packed format: the returned array must have ``n``
        columns and ``uband + lband + 1`` rows in which Jacobian diagonals are
        written. Specifically ``jac_packed[uband + i - j , j] = jac[i, j]``.
        The same format is used in `scipy.linalg.solve_banded` (check for an
        illustration).  These parameters can be also used with ``jac=None`` to
        reduce the number of Jacobian elements estimated by finite differences.
    min_step : float, optional
        The minimum allowed step size for 'LSODA' method. 
        By default `min_step` is zero.
    Returns
    -------
    Bunch object with the following fields defined:
    t : ndarray, shape (n_points,)
        Time points.
    y : ndarray, shape (n, n_points)
        Values of the solution at `t`.
    sol : `OdeSolution` or None
        Found solution as `OdeSolution` instance; None if `dense_output` was
        set to False.
    t_events : list of ndarray or None
        Contains for each event type a list of arrays at which an event of
        that type event was detected. None if `events` was None.
    y_events : list of ndarray or None
        For each value of `t_events`, the corresponding value of the solution.
        None if `events` was None.
    nfev : int
        Number of evaluations of the right-hand side.
    njev : int
        Number of evaluations of the Jacobian.
    nlu : int
        Number of LU decompositions.
    status : int
        Reason for algorithm termination:
            * -1: Integration step failed.
            *  0: The solver successfully reached the end of `tspan`.
            *  1: A termination event occurred.
    message : string
        Human-readable description of the termination reason.
    success : bool
        True if the solver reached the interval end or a termination event
        occurred (``status >= 0``).
    References
    ----------
    .. [1] J. R. Dormand, P. J. Prince, "A family of embedded Runge-Kutta
           formulae", Journal of Computational and Applied Mathematics, Vol. 6,
           No. 1, pp. 19-26, 1980.
    .. [2] L. W. Shampine, "Some Practical Runge-Kutta Formulas", Mathematics
           of Computation,, Vol. 46, No. 173, pp. 135-150, 1986.
    .. [3] P. Bogacki, L.F. Shampine, "A 3(2) Pair of Runge-Kutta Formulas",
           Appl. Math. Lett. Vol. 2, No. 4. pp. 321-325, 1989.
    .. [4] E. Hairer, G. Wanner, "Solving Ordinary Differential Equations II:
           Stiff and Differential-Algebraic Problems", Sec. IV.8.
    .. [5] `Backward Differentiation Formula
            <https://en.wikipedia.org/wiki/Backward_differentiation_formula>`_
            on Wikipedia.
    .. [6] L. F. Shampine, M. W. Reichelt, "THE MATLAB ODE SUITE", SIAM J. SCI.
           COMPUTE., Vol. 18, No. 1, pp. 1-22, January 1997.
    .. [7] A. C. Hindmarsh, "ODEPACK, A Systematized Collection of ODE
           Solvers," IMACS Transactions on Scientific Computation, Vol 1.,
           pp. 55-64, 1983.
    .. [8] L. Petzold, "Automatic selection of methods for solving stiff and
           nonstiff systems of ordinary differential equations", SIAM Journal
           on Scientific and Statistical Computing, Vol. 4, No. 1, pp. 136-148,
           1983.
    .. [9] `Stiff equation <https://en.wikipedia.org/wiki/Stiff_equation>`_ on
           Wikipedia.
    .. [10] A. Curtis, M. J. D. Powell, and J. Reid, "On the estimation of
            sparse Jacobian matrices", Journal of the Institute of Mathematics
            and its Applications, 13, pp. 117-120, 1974.
    .. [11] `Cauchy-Riemann equations
             <https://en.wikipedia.org/wiki/Cauchy-Riemann_equations>`_ on
             Wikipedia.
    .. [12] `Lotka-Volterra equations
            <https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations>`_
            on Wikipedia.
    .. [13] E. Hairer, S. P. Norsett G. Wanner, "Solving Ordinary Differential
            Equations I: Nonstiff Problems", Sec. II.
    .. [14] `Page with original Fortran code of DOP853
            <http://www.unige.ch/~hairer/software.html>`_.
    Examples
    --------
    Basic exponential decay showing automatically chosen time points.
    >>> from scipy.integrate import solve_ivp
    >>> def exponential_decay(t, y): return -0.5 * y
    >>> sol = solve_ivp(exponential_decay, [0, 10], [2, 4, 8])
    >>> print(sol.t)
    [ 0.          0.11487653  1.26364188  3.06061781  4.81611105  6.57445806
      8.33328988 10.        ]
    >>> print(sol.y)
    [[2.         1.88836035 1.06327177 0.43319312 0.18017253 0.07483045
      0.03107158 0.01350781]
     [4.         3.7767207  2.12654355 0.86638624 0.36034507 0.14966091
      0.06214316 0.02701561]
     [8.         7.5534414  4.25308709 1.73277247 0.72069014 0.29932181
      0.12428631 0.05403123]]
    Specifying points where the solution is desired.
    >>> sol = solve_ivp(exponential_decay, [0, 10], [2, 4, 8],
    ...                 t_eval=[0, 1, 2, 4, 10])
    >>> print(sol.t)
    [ 0  1  2  4 10]
    >>> print(sol.y)
    [[2.         1.21305369 0.73534021 0.27066736 0.01350938]
     [4.         2.42610739 1.47068043 0.54133472 0.02701876]
     [8.         4.85221478 2.94136085 1.08266944 0.05403753]]
    Cannon fired upward with terminal event upon impact. The ``terminal`` and
    ``direction`` fields of an event are applied by monkey patching a function.
    Here ``y[0]`` is position and ``y[1]`` is velocity. The projectile starts
    at position 0 with velocity +10. Note that the integration never reaches
    t=100 because the event is terminal.
    >>> def upward_cannon(t, y): return [y[1], -0.5]
    >>> def hit_ground(t, y): return y[0]
    >>> hit_ground.terminal = True
    >>> hit_ground.direction = -1
    >>> sol = solve_ivp(upward_cannon, [0, 100], [0, 10], events=hit_ground)
    >>> print(sol.t_events)
    [array([40.])]
    >>> print(sol.t)
    [0.00000000e+00 9.99900010e-05 1.09989001e-03 1.10988901e-02
     1.11088891e-01 1.11098890e+00 1.11099890e+01 4.00000000e+01]
    Use `dense_output` and `events` to find position, which is 100, at the apex
    of the cannonball's trajectory. Apex is not defined as terminal, so both
    apex and hit_ground are found. There is no information at t=20, so the sol
    attribute is used to evaluate the solution. The sol attribute is returned
    by setting ``dense_output=True``. Alternatively, the `y_events` attribute
    can be used to access the solution at the time of the event.
    >>> def apex(t, y): return y[1]
    >>> sol = solve_ivp(upward_cannon, [0, 100], [0, 10], 
    ...                 events=(hit_ground, apex), dense_output=True)
    >>> print(sol.t_events)
    [array([40.]), array([20.])]
    >>> print(sol.t)
    [0.00000000e+00 9.99900010e-05 1.09989001e-03 1.10988901e-02
     1.11088891e-01 1.11098890e+00 1.11099890e+01 4.00000000e+01]
    >>> print(sol.sol(sol.t_events[1][0]))
    [100.   0.]
    >>> print(sol.y_events)
    [array([[-5.68434189e-14, -1.00000000e+01]]), array([[1.00000000e+02, 1.77635684e-15]])]
    As an example of a system with additional parameters, we'll implement
    the Lotka-Volterra equations [12]_.
    >>> def lotkavolterra(t, z, a, b, c, d):
    ...     x, y = z
    ...     return [a*x - b*x*y, -c*y + d*x*y]
    ...
    We pass in the parameter values a=1.5, b=1, c=3 and d=1 with the `args`
    argument.
    >>> sol = solve_ivp(lotkavolterra, [0, 15], [10, 5], args=(1.5, 1, 3, 1),
    ...                 dense_output=True)
    Compute a dense solution and plot it.
    >>> t = np.linspace(0, 15, 300)
    >>> z = sol.sol(t)
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(t, z.T)
    >>> plt.xlabel('t')
    >>> plt.legend(['x', 'y'], shadow=True)
    >>> plt.title('Lotka-Volterra System')
    >>> plt.show()
    """
    __name__ = "solve_ivpd"
    if method not in METHODS and not (
            inspect.isclass(method) and issubclass(method, OdeSolver)):
        raise ValueError("`method` must be one of {} or OdeSolver class."
                         .format(METHODS))

    t0, tf = float(t_span[0]), float(t_span[1])

    if args is not None:
        # Wrap the user's fun (and jac, if given) in lambdas to hide the
        # additional parameters.  Pass in the original fun as a keyword
        # argument to keep it in the scope of the lambda.
        fun = lambda t, x, fun=fun: fun(t, x, *args)
        jac = options.get('jac')
        if callable(jac):
            options['jac'] = lambda t, x: jac(t, x, *args)

    if t_eval is not None:
        t_eval = np.asarray(t_eval)
        if t_eval.ndim != 1:
            raise ValueError("`t_eval` must be 1-dimensional.")

        if np.any(t_eval < min(t0, tf)) or np.any(t_eval > max(t0, tf)):
            raise ValueError("Values in `t_eval` are not within `t_span`.")

        d = np.diff(t_eval)
        if tf > t0 and np.any(d <= 0) or tf < t0 and np.any(d >= 0):
            raise ValueError("Values in `t_eval` are not properly sorted.")

        if tf > t0:
            t_eval_i = 0
        else:
            # Make order of t_eval decreasing to use np.searchsorted.
            t_eval = t_eval[::-1]
            # This will be an upper bound for slices.
            t_eval_i = t_eval.shape[0]

    if method in METHODS:
        method = METHODS[method]

    solver = method(fun, t0, y0, tf, vectorized=vectorized, **options)

    if t_eval is None:
        # print(y0)
        ts = [t0]
        ys = [dm(y0[0], prec)]
    elif t_eval is not None and dense_output:
        ts = []
        ti = [t0]
        ys = []
    else:
        ts = []
        ys = []

    interpolants = []

    events, is_terminal, event_dir = prepare_events(events)

    if events is not None:
        if args is not None:
            # Wrap user functions in lambdas to hide the additional parameters.
            # The original event function is passed as a keyword argument to the
            # lambda to keep the original function in scope (i.e. avoid the
            # late binding closure "gotcha").
            events = [lambda t, x, event=event: event(t, x, *args)
                      for event in events]
        g = [event(t0, y0) for event in events]
        t_events = [[] for _ in range(len(events))]
        y_events = [[] for _ in range(len(events))]
    else:
        t_events = None
        y_events = None

    status = None
    while status is None:
        message = solver.step()

        if solver.status == 'finished':
            status = 0
        elif solver.status == 'failed':
            status = -1
            break

        t_old = solver.t_old
        t = solver.t
        y = solver.y

        if dense_output:
            sol = solver.dense_output()
            interpolants.append(sol)
        else:
            sol = None

        if events is not None:
            g_new = [event(t, y) for event in events]
            active_events = find_active_events(g, g_new, event_dir)
            if active_events.size > 0:
                if sol is None:
                    sol = solver.dense_output()

                root_indices, roots, terminate = handle_events(
                    sol, events, active_events, is_terminal, t_old, t)

                for e, te in zip(root_indices, roots):
                    t_events[e].append(te)
                    y_events[e].append(sol(te))

                if terminate:
                    status = 1
                    t = roots[-1]
                    y = sol(t)

            g = g_new

        if t_eval is None:
            ts.append(t)
            ys.append(dm(y, prec))
        else:
            # The value in t_eval equal to t will be included.
            if solver.direction > 0:
                t_eval_i_new = np.searchsorted(t_eval, t, side='right')
                t_eval_step = t_eval[t_eval_i:t_eval_i_new]
            else:
                t_eval_i_new = np.searchsorted(t_eval, t, side='left')
                # It has to be done with two slice operations, because
                # you can't slice to 0-th element inclusive using backward
                # slicing.
                t_eval_step = t_eval[t_eval_i_new:t_eval_i][::-1]

            if t_eval_step.size > 0:
                if sol is None:
                    sol = solver.dense_output()
                ts.append(t_eval_step)
                ys.append(sol(t_eval_step).tolist()[-1])
                t_eval_i = t_eval_i_new
        
        if t_eval is not None and dense_output:
            ti.append(t)

    message = MESSAGES.get(status, message)

    if t_events is not None:
        t_events = [np.asarray(te) for te in t_events]
        y_events = [np.asarray(ye) for ye in y_events]

    if t_eval is None:
        ts = np.array(ts)
        ys = df.DF([ys], prec).TD()
    else:
        ts = np.hstack(ts)
        ys = df.DF([np.hstack(ys)], prec).decfunc()

    if dense_output:
        if t_eval is None:
            sol = OdeSolution(ts, interpolants)
        else:
            sol = OdeSolution(ti, interpolants)
    else:
        sol = None

    return OdeResult(t=ts, y=ys, sol=sol, t_events=t_events, y_events=y_events,
                     nfev=solver.nfev, njev=solver.njev, nlu=solver.nlu,
                     status=status, message=message, success=status >= 0)
       
def exponential_decay(t, y, args): 
       pre = args
       try:
              yn = [dm(-0.5, pre)*dm(i, pre) for i in y]
              return yn
              
       except:
              yn = dm(-0.5, pre)*dm(y, pre)
              return [yn]




