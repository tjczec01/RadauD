# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 11:46:52 2020

Github: https://github.com/tjczec01

@author: Travis J Czechorski

E-mail: tjczec01@gmail.com
"""
from __future__ import division, print_function, absolute_import
import os
import sys
import inspect
import numpy as np
from scipy.optimize import OptimizeResult, brentq
from RadauM import RadauM
from base import OdeSolver
from common import OdeSolution

EPS = 2.220446049250313e-16

plist = list(sys.path)
ppath = os.path.dirname(os.path.abspath(__file__))
PTH = '{}'.format(ppath)

if PTH in plist:
    pass
else:
    sys.path.append(r'{}'.format(ppath))

# __all__ = ["OdeResult", "prepare_events", "solve_event_equation",
#            "handle_events", "find_active_events", "solve_ivpm"]

METHODS = {'RadauM' : RadauM}

MESSAGES = {0: "The solver successfully reached the end of the integration interval.",
            1: "A termination event occurred."}

class OdeResult(OptimizeResult):
    """Continuous ODE solution.
    It is organized as a collection of `DenseOutput` objects which represent
    local interpolants. It provides an algorithm to select a right interpolant
    for each given point.
    The interpolants cover the range between `t_min` and `t_max` (see
    Attributes below). Evaluation outside this interval is not forbidden, but
    the accuracy is not guaranteed.
    When evaluating at a breakpoint (one of the values in `ts`) a segment with
    the lower index is selected.
    Parameters
    ----------
    ts : array_like, shape (n_segments + 1,)
        Time instants between which local interpolants are defined. Must
        be strictly increasing or decreasing (zero segment with two points is
        also allowed).
    interpolants : list of DenseOutput with n_segments elements
        Local interpolants. An i-th interpolant is assumed to be defined
        between ``ts[i]`` and ``ts[i + 1]``.
    Attributes
    ----------
    t_min, t_max : float
        Time range of the interpolation.
    """
    __name__ = "OdeResult"
    # pass


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


def solve_event_equation(event, sol, t_old, t_v):
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

    return brentq(lambda t_v: event(t_v, sol(t_v)), t_old, t_v,
                  xtol=4.0 * EPS, rtol=4.0 * EPS)


def handle_events(sol, events, active_events, is_terminal,
                  t_l):
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
    t_old = t_l[0]
    t_v = t_l[1]
    roots = [solve_event_equation(events[event_index], sol, t_old, t_v)
             for event_index in active_events]
    roots = np.asarray(roots)
    if np.any(is_terminal[active_events]):
        if t_v > t_old:
            order = np.argsort(roots)
        else:
            order = np.argsort(-roots)
        active_events = active_events[order]
        roots = roots[order]
        t_v = np.nonzero(is_terminal[active_events])[0][0]
        active_events = active_events[:t_v + 1]
        roots = roots[:t_v + 1]
        terminate = True
    else:
        terminate = False
    return active_events, roots, terminate


def find_active_events(g_0, g_new, direction):
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
    g_0, g_new = np.asarray(g_0), np.asarray(g_new)
    u_p = (g_0 <= 0) & (g_new >= 0)
    down = (g_0 >= 0) & (g_new <= 0)
    either = u_p | down
    mask = (u_p & (direction > 0) |
            down & (direction < 0) |
            either & (direction == 0))

    return np.nonzero(mask)[0]


# def solve_ivpm(fun, t_span, y0, method='RK45', t_eval=None, dense_output=False,
#     solve_ivpm.__name__ = "solve_ivpm"
#     prec = args
#     if method not in METHODS and not (
#             inspect.isclass(method) and issubclass(method, OdeSolver)):
#         raise ValueError("`method` must be one of {} or OdeSolver class."
#                           .format(METHODS))

#     t0, tf = float(t_span[0]), float(t_span[1])

#     if args is not None:
#         # Wrap the user's fun (and jac, if given) in lambdas to hide the
#         # additional parameters.  Pass in the original fun as a keyword
#         # argument to keep it in the scope of the lambda.
#         fun = lambda t, x, fun=fun: fun(t, x, *args)
#         jac = options.get('jac')
#         if callable(jac):
#             options['jac'] = lambda t, x: jac(t, x, *args)

#     if t_eval is not None:
#         t_eval = np.asarray(t_eval)
#         if t_eval.ndim != 1:
#             raise ValueError("`t_eval` must be 1-dimensional.")

#         if np.any(t_eval < min(t0, tf)) or np.any(t_eval > max(t0, tf)):
#             raise ValueError("Values in `t_eval` are not within `t_span`.")

#         d = np.diff(t_eval)
#         if tf > t0 and np.any(d <= 0) or tf < t0 and np.any(d >= 0):
#             raise ValueError("Values in `t_eval` are not properly sorted.")

#         if tf > t0:
#             t_eval_i = 0
#         else:
#             # Make order of t_eval decreasing to use np.searchsorted.
#             t_eval = t_eval[::-1]
#             # This will be an upper bound for slices.
#             t_eval_i = t_eval.shape[0]

#     if method in METHODS:
#         method = METHODS[method]

#     solver = method(fun, t0, y0, tf, vectorized=vectorized, **options)

#     if t_eval is None:
#         # print(y0)
#         ts = [t0]
#         y_s = [y0]
#     elif t_eval is not None and dense_output:
#         ts = []
#         ti = [t0]
#         y_s = []
#     else:
#         ts = []
#         y_s = []

#     interpolants = []

#     events, is_terminal, event_dir = prepare_events(events)

#     if events is not None:
#         if args is not None:
#             # Wrap user functions in lambdas to hide the additional parameters.
#             # The original event function is passed as a keyword argument to the
#             # lambda to keep the original function in scope (i.e. avoid the
#             # late binding closure "gotcha").
#             events = [lambda t, x, event=event: event(t, x, *args)
#                       for event in events]
#         g = [event(t0, y0) for event in events]
#         t_events = [[] for _ in range(len(events))]
#         y_events = [[] for _ in range(len(events))]
#     else:
#         t_events = None
#         y_events = None

#     status = None
#     while status is None:
#         message = solver.step()

#         if solver.status == 'finished':
#             status = 0
#         elif solver.status == 'failed':
#             status = -1
#             break

#         t_old = solver.t_old
#         t = solver.t
#         y = solver.y

#         if dense_output:
#             sol = solver.dense_output()
#             interpolants.append(sol)
#         else:
#             sol = None

#         if events is not None:
#             g_new = [event(t, y) for event in events]
#             active_events = find_active_events(g, g_new, event_dir)
#             if active_events.size > 0:
#                 if sol is None:
#                     sol = solver.dense_output()

#                 root_indices, roots, terminate = handle_events(
#                     sol, events, active_events, is_terminal, t_old, t)

#                 for e, te in zip(root_indices, roots):
#                     t_events[e].append(te)
#                     y_events[e].append(sol(te))

#                 if terminate:
#                     status = 1
#                     t = roots[-1]
#                     y = sol(t)

#             g = g_new

#         if t_eval is None:
#             ts.append(t)
#             y_s.append(y)
#         else:
#             # The value in t_eval equal to t will be included.
#             if solver.direction > 0:
#                 t_eval_i_new = np.searchsorted(t_eval, t, side='right')
#                 t_eval_step = t_eval[t_eval_i:t_eval_i_new]
#             else:
#                 t_eval_i_new = np.searchsorted(t_eval, t, side='left')
#                 # It has to be done with two slice operations, because
#                 # you can't slice to 0-th element inclusive using backward
#                 # slicing.
#                 t_eval_step = t_eval[t_eval_i_new:t_eval_i][::-1]

#             if t_eval_step.size > 0:
#                 if sol is None:
#                     sol = solver.dense_output()
#                 ts.append(t_eval_step)
#                 y_s.append(sol(t_eval_step))
#                 t_eval_i = t_eval_i_new

#         if t_eval is not None and dense_output:
#             ti.append(t)

#     message = MESSAGES.get(status, message)

#     if t_events is not None:
#         t_events = [np.asarray(te) for te in t_events]
#         y_events = [np.asarray(ye) for ye in y_events]

#     if t_eval is None:
#         ts = np.array(ts)
#         y_s = np.vstack(y_s).T
#     else:
#         ts = np.hstack(ts)
#         y_s = np.hstack(y_s)

#     if dense_output:
#         if t_eval is None:
#             sol = OdeSolution(ts, interpolants)
#         else:
#             sol = OdeSolution(ti, interpolants)
#     else:
#         sol = None

#     return OdeResult(t=ts, y=y_s, sol=sol, t_events=t_events, y_events=y_events,
#                       nfev=solver.nfev, njev=solver.njev, nlu=solver.nlu,
#                       status=status, message=message, success=status >= 0)

def solve_ivpm(fun, t_span, y_0, method='RadauM', t_eval=None, dense_output=False,
               events=None, vectorized=False, args=None, **options):
    """Solve an initial value problem for a sy_stem of ODEs.
    This function numerically integrates a sy_stem of ordinary differential
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
    Another option alway_s available is to rewrite your problem for real and
    imaginary parts separately.
    Parameters
    ----------
    fun : callable
        Right-hand side of the sy_stem. The calling signature is ``fun(t, y)``.
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
        Jacobian matrix of the right-hand side of the sy_stem with respect
        to y, required by the 'Radau', 'BDF' and 'LSODA' method. The
        Jacobian matrix has shape (n, n) and its element (i, j) is equal to
        ``d f_i / d y_j``.  There are three way_s to define the Jacobian:
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
        a corresponding element in the Jacobian is alway_s zero. If None
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
        Contains for each event type a list of array_s at which an event of
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
    .. [7] A. C. Hindmarsh, "ODEPACK, A Sy_stematized Collection of ODE
            Solvers," IMACS Transactions on Scientific Computation, Vol 1.,
            pp. 55-64, 1983.
    .. [8] L. Petzold, "Automatic selection of methods for solving stiff and
            nonstiff sy_stems of ordinary differential equations", SIAM Journal
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
    As an example of a sy_stem with additional parameters, we'll implement
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
    >>> plt.title('Lotka-Volterra Sy_stem')
    >>> plt.show()
    """
    if method not in METHODS and not (inspect.isclass(method) and issubclass(method, OdeSolver)):
        raise ValueError("`method` must be one of {} or OdeSolver class.".format(METHODS))

    t_0, t_f = float(t_span[0]), float(t_span[1])

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

        if np.any(t_eval < min(t_0, t_f)) or np.any(t_eval > max(t_0, t_f)):
            raise ValueError("Values in `t_eval` are not within `t_span`.")

        d_p = np.diff(t_eval)
        if t_f > t_0 and np.any(d_p <= 0) or t_f < t_0 and np.any(d_p >= 0):
            raise ValueError("Values in `t_eval` are not properly sorted.")

        if t_f > t_0:
            t_eval_i = 0
        else:
            # Make order of t_eval decreasing to use np.searchsorted.
            t_eval = t_eval[::-1]
            # This will be an upper bound for slices.
            t_eval_i = t_eval.shape[0]

    if method in METHODS:
        method = METHODS.get('RadauM', [RadauM])

    solver = method(fun, t_0, y_0, t_f, vectorized=vectorized, **options)

    if t_eval is None:
        # print(y0)
        t_s = [t_0]
        y_s = [y_0]
    elif t_eval is not None and dense_output:
        t_s = []
        t_i = [t_0]
        y_s = []
    else:
        t_s = []
        y_s = []

    interpolants = []

    events, is_terminal, event_dir = prepare_events(events)

    if events is not None:
        if args is not None:
            # Wrap user functions in lambdas to hide the additional parameters.
            # The original event function is passed as a keyword argument to the
            # lambda to keep the original function in scope (i.e. avoid the
            # late binding closure "gotcha").
            events = [lambda t_v, x, event=event: event(t_v, x, *args)
                      for event in events]
        g_0 = [event(t_0, y_0) for event in events]
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
        t_v = solver.t
        y_v = solver.y

        if dense_output:
            sol = solver.dense_output()
            interpolants.append(sol)
        else:
            sol = None

        if events is not None:
            g_new = [event(t_v, y_v) for event in events]
            active_events = find_active_events(g_0, g_new, event_dir)
            if active_events.size > 0:
                if sol is None:
                    sol = solver.dense_output()
                t_l = [t_old, t_v]
                root_indices, roots, terminate = handle_events(
                    sol, events, active_events, is_terminal, t_l)

                for e_v, t_e in zip(root_indices, roots):
                    t_events[e_v].append(t_e)
                    y_events[e_v].append(sol(t_e))

                if terminate:
                    status = 1
                    t_v = roots[-1]
                    y_v = sol(t_v)

            g_0 = g_new

        if t_eval is None:
            t_s.append(t_v)
            y_s.append(y_v)
        else:
            # The value in t_eval equal to t will be included.
            if solver.direction > 0:
                t_eval_i_new = np.searchsorted(t_eval, t_v, side='right')
                t_eval_step = t_eval[t_eval_i:t_eval_i_new]
            else:
                t_eval_i_new = np.searchsorted(t_eval, t_v, side='left')
                # It has to be done with two slice operations, because
                # you can't slice to 0-th element inclusive using backward
                # slicing.
                t_eval_step = t_eval[t_eval_i_new:t_eval_i][::-1]

            if t_eval_step.size > 0:
                if sol is None:
                    sol = solver.dense_output()
                t_s.append(t_eval_step)
                y_s.append(sol(t_eval_step))
                t_eval_i = t_eval_i_new

        if t_eval is not None and dense_output:
            t_i.append(t_v)

    message = MESSAGES.get(status, message)

    if t_events is not None:
        t_events = [np.asarray(t_e) for t_e in t_events]
        y_events = [np.asarray(ye) for ye in y_events]

    if t_eval is None:
        t_s = np.array(t_s)
        y_s = np.vstack(y_s).T
    else:
        t_s = np.hstack(t_s)
        y_s = np.hstack(y_s)

    if dense_output:
        if t_eval is None:
            sol = OdeSolution(t_s, interpolants)
        else:
            sol = OdeSolution(t_i, interpolants)
    else:
        sol = None

    return OdeResult(t=t_s, y=y_s, sol=sol, t_events=t_events, y_events=y_events,
                     nfev=solver.nfev, njev=solver.njev, nlu=solver.nlu,
                     status=status, message=message, success=status >= 0)
