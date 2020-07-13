# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 12:10:55 2020

Github: https://github.com/tjczec01

@author: Travis J Czechorski

E-mail: tjczec01@gmail.com
"""
import os
dir_path = os.path.dirname(os.path.realpath(__file__))
exec(open(r'{}\__init__.py'.format(dir_path)).read())
from scipy.integrate import solve_ivp
from ivpd import solve_ivpd
from ivpm import solve_ivpm
from common import dm
import mpmath as mp
import numpy as np
import matplotlib.pyplot as plt
import math
from tqdm import tqdm
import sympy as sp
from sympy import Heaviside, DiracDelta

def mv(v, pr):
        mp.dps = int(pr)
        try:
               mv = mp.mpf(v)
        except:
               mv = mp.mpc(v.real, v.imag)
        return mv

def get_change(current, previous):
    if current == previous:
        return mv(100.0, pi[0])
    try:
        return mv(mv(abs(current - previous) / mv(current, pi[0]) , pi[0]) , pi[0])
    except ZeroDivisionError:
        return 0

def flatten(lm):
              flatten = lambda l: [item for sublist in l for item in sublist]
              return flatten(lm)

def RHS(t, y, args):
    eq1 = 2.0*y[0]**3.0
    eq2 = 3.0*y[1]**2.0
    return [eq1, eq2]

def jacob(t, y, args):
    return [[6.0*y[0]**2.0, 0.0], [0.0, 6.0*y[1]]]

def RHSd(t, y, args):
    pre = args
    pre = args
    Y0 = dm(y[0], pre)
    Y1 = dm(y[1], pre)
    eq1 = dm(2.0, pre)*Y0**(dm(3.0, pre))
    eq2 = dm(3.0, pre)*Y1**(dm(2.0, pre))
    return [eq1, eq2]


def jacobd(t, y, args):
    pre = args
    Y0 = dm(y[0], pre)
    Y1 = dm(y[1], pre)
    return [[dm(6.0, pre)*Y0**dm(2.0, pre), dm(0.0, pre)], [dm(0.0, pre), dm(6.0, pre)*Y1]]



def RHSm(t, y, args):
    pre = args
    Y0 = mv(y[0], pre)
    Y1 = mv(y[1], pre)
    eq1 = mv(2.0, pre)*Y0**(mv(3.0, pre))
    eq2 = mv(3.0, pre)*Y1**(mv(2.0, pre))
    return [eq1, eq2]

def jacobm(t, y, args):
    pre = args
    Y0 = mv(y[0], pre)
    Y1 = mv(y[1], pre)
    return [[mv(6.0, pre)*Y0**mv(2.0, pre), mv(0.0, pre)], [mv(0.0, pre), mv(6.0, pre)*Y1]]


pi = [50]
mp.dps = pi[0]
init = 1.00000
tevs = [i/100 for i in range(0, 51, 1)]
sol1 = solve_ivp(RHS, [0.0, 0.5], [1.0, 1.5], t_eval=tevs, method="Radau", args=(pi), jac=jacob, rtol=1E-9, atol=1E-9)
# sol2 = solve_ivpd(RHSd, [0.0, 0.5], [dm(1.0, pi[0]), dm(1.5, pi[0])], t_eval=tevs, method="RadauD", prec=pi[0], args=(pi), jac=jacobd, rtol=1E-9, atol=1E-9) # V
# print(sol2.y)
# print(sol2.y[1])
sol3 = solve_ivpm(RHSm, [0.0, 0.5], [mv(1.0, pi[0]), mv(1.5, pi[0])],  method="RadauM", t_eval=tevs, prec=pi[0], args=(pi), jac=jacobm, rtol=1E-9, atol=1E-9)
# print(sol1.t.tolist())
# time_l = sol1.t.tolist()

# s2 = flatten(sol2.y)
# s3 = sol3.y
s1a = sol1.y[0].tolist()
s1b = sol1.y[-1].tolist()
# print(sol3.y)
print(s1a)
# print(sol2.y[0])
for i in range(len(sol3.y[0])):
    print(mpf(sol3.y[0][i]))
print(sol3.y[0])
print(s1b)
# print(sol2.y[1])
print(sol3.y[1])
# print(sol2.y)

# print(sol2.y[1])

# print(sol1.y[0].tolist())
# print(sol2.t.tolist())
# print(flatten(sol2.y))
# print(sol3.t.tolist())
# mp.nprint(flatten(sol3.y), pi[0])
# print(flatten(sol3.y))
# print("")
# ch = [get_change(mv(str(i), pi[0]), j) for i, j in zip(s2, s3)]
# ch2 = [mv(get_change(mv(str(i), pi[0]), j), pi[0]) for i, j in zip(s1, s3)]
# ch3 = [mv(get_change(mv(str(i), pi[0]), mv(str(j), pi[0])), pi[0]) for i, j in zip(s1, s2)]

# ch4 = [i - j for i, j in zip(ch2, ch3)]
# ch5 = [i - j  for i, j in zip(ch, ch2)]
# ch6 = [i - j  for i, j in zip(ch, ch3)]
# ch7 = [abs((mv(str(i), pi[0])/j) - mv(1.0, pi[0]))*mv(100.0, pi[0]) for i, j in zip(s2, s3)]
# ch8 = [abs((mv(str(i), pi[0])/j) - mv(1.0, pi[0]))*mv(100.0, pi[0]) for i, j in zip(s3, s2)]
# ch9 = [abs((mv(str(i), pi[0])/mv(str(j), pi[0])) - mv(1.0, pi[0]))*mv(100.0, pi[0]) for i, j in zip(s1, s2)]
# ch10 = [abs((i/j) - mv(1.0, pi[0]))*mv(100.0, pi[0]) for i, j in zip(s1, s3)]

# fig = plt.figure()
# plt.plot(time_l, ch7, 'k-', label='RadauD/RadauM')
# plt.plot(time_l, ch8, 'r--', label='RadauD/RadauM')
# plt.legend([r'$\frac{Radau-D}{Radau-M}$', r'$\frac{Radau-M}{Radau-D}$'], loc="best")
# plt.xlabel('Time')
# plt.ylabel('Value [%]')
# plt.title('Decimal vs. Mpmath')
# plt.grid()
# plt.show()

# fig = plt.figure()
# plt.plot(time_l, ch9, 'g-', label='Radau/RadauD')
# plt.plot(time_l, ch10, 'b--', label='RadauD/RadauM')
# plt.legend([r'$\frac{Radau}{Radau-D}$', r'$\frac{Radau}{Radau-M}$'], loc="best")
# plt.xlabel('Time')
# plt.ylabel('Value [%]')
# plt.title('Decimal vs. Mpmath')
# plt.grid()
# plt.show()

# for i, j in zip(s2, s3):
#       print("Decimal Length {} : Mpf Length {}".format(len(str(j)), len(str(i))))
#       print("{} : {}".format(j, i))

clear = lambda: os.system('cls')
cwd = os.getcwd()
dir_path = os.path.dirname(os.path.realpath(__file__))
path_fol = "{}\Modified SEIRS Model Validation".format(dir_path)

try:
       os.mkdir(path_fol)
except:
       pass

def roundup(x, places):
    return int(math.ceil(x / int(places))) * int(places)

def R0(α, β, μ, γ):
       R_0 = (α/(μ + α))*(β/(μ + γ))
       return R_0

def R_t(σ, R_0):
       return (1 - σ)*R_0

def R_s(t, *args):
       τ_ξ, ξ, times, R = args
       te = t - τ_ξ
       teindex = min(range(len(times)), key=lambda i: abs(times[i] - te)) - 1
       if teindex == len(R):
              teindex -= 1
       elif teindex >= len(R):
              while teindex >= len(R):
                     teindex -= 1
       return ξ*float(R[teindex])

def R_c(t, τ_ξ, ξ):
       tt = t - τ_ξ
       return ξ*tt

def τ_σf(τ_pre, τ_post):
       return τ_pre + τ_post

def τ_in(α):
       return α**-1

def τ_re(γ):
       return γ**-1

def τ_ho(σ):
       return σ**-1

def alpha(τ_inc):
       return τ_inc**-1

def beta(R_0, γ):
       return R_0*γ

def betaf(R_0, α, γ, μ):
       val1 = R_0*(μ + α)*(μ + γ)
       val2 = val1/α
       return val2

def gamma(τ_rec):
       return τ_rec**-1

def R_sm(t, *args):
       τ_ξ, ξ, times, R = args
       te = t - float(τ_ξ)
       teindex = min(range(len(times)), key=lambda i: abs(times[i] - te)) - 1
       if teindex == len(R):
              teindex -= 1
       elif teindex >= len(R):
              while teindex >= len(R):
                     teindex -= 1
       return ξ*R[teindex]

# print(betaf(2.8, 0.48, 0.0205, 3**-1))
# Λ is the birth rate in the overall population
# µ is the death rate due to conditions other than the COVID-19
# β is the rate of transmission per S-I contact
# α is the rate of which an exposed person becomes infected
# γ is the recovery rate
# σ is the efficiency of the control action
# ξ is the percentage of the recovered population who are resusceptible
# κ represents the percentage of non-elderly who recovered
# κ_old represents the percentage of non-elderly and elderly who recovered
# τ_(pre - σ) represents the time to initiate the control action after the first confirmed case at t = 0
# τ_(post - σ) represents the time after the control action has been initiated but before the effects are evidenced in the outputs of the system
# τ_inc incubation time
# τ_rec recovery time
# τ_hos represents the time spent hospitalised
# τ_ξ represents the duration of temporary immunity
# N represents the stock population
# N_old represents the percentage of elderly population (above 65 years of age)
# R_0 represents the basic reproduction number
# t represents the time

def SEIRS(t, y, *args):
       σ, β, γ, α, Λ, μ, ξ, κ, κ_old, τ_ξ, τ_σ, N, N_old, time, Is, Ss, Rs  = args
       args2 = [τ_ξ, ξ, time, Rs]
       S = y[0]
       E = y[1]
       I = y[2]
       R = y[3]
       D = y[4]
       if t >= τ_σ and t >= τ_ξ:
              tindex = min(range(len(time)), key=lambda i: abs(time[i]- t)) - 1  # time.index(tt)
              if tindex == len(Is):
                     tindex -= 1
              elif tindex >= len(Is):
                     tindex = len(Is) - 1
              It = Is[tindex]
              St = Ss[tindex]
              dsdt = Λ - μ*S - (β*I*S)/N + (σ*β*It*St)/N + R_s(t, *args2)
              drdt = (γ + μ)*I - μ*R - R_s(t, *args2)
       if t >= τ_σ and t < τ_ξ:
               tindex = min(range(len(time)), key=lambda i: abs(time[i]- t)) - 1 # time.index(tt)
               if tindex == len(Is):
                     tindex -= 1
               elif tindex >= len(Is):
                      tindex = len(Is) - 1
               It = Is[tindex]
               St = Ss[tindex]
               dsdt = Λ - μ*S - (β*I*S)/N + (σ*β*It*St)/N
               drdt = (γ + μ)*I - μ*R
       elif t < τ_σ:
              It = 0
              St = 0
              dsdt = Λ - μ*S - (β*I*S)/N
              drdt = (γ + μ)*I - μ*R
       dedt = (β*I*S)/N - (σ*β*It*St)/N - (μ + α)*E
       didt = (μ + α)*E - (γ + μ)*I - γ*((1 - κ_old)*N_old + (1 - κ)*(1 - N_old))*I
       dDdt = γ*((1 - κ_old)*N_old + (1 - κ)*(1 - N_old))*I
       return [dsdt, dedt, didt, drdt, dDdt]



# def SEIRSD(t, y, *args):
#        σ, β, γ, α, Λ, μ, ξ, κ, κ_old, τ_ξ, τ_σ, N, N_old, time, Is, Ss, Rs, pre  = args
#        args2 = [τ_ξ, ξ, time, Rs]
#        print(t)
#        print(y)
#        S = dm(y[0], pre)
#        E = dm(y[1], pre)
#        I = dm(y[2], pre)
#        R = dm(y[3], pre)
#        D = dm(y[4], pre)
#        if t >= τ_σ and t >= τ_ξ:
#               tindex = min(range(len(time)), key=lambda i: abs(time[i]- t)) - 1  # time.index(tt)
#               if tindex == len(Is):
#                      tindex -= 1
#               elif tindex >= len(Is):
#                      tindex = len(Is) - 1
#               It = Is[tindex]
#               St = Ss[tindex]
#               dsdt = dm(Λ, pre) - dm(μ, pre)*dm(S, pre) - (dm(β, pre)*dm(I, pre)*dm(S, pre))/dm(N, pre) + (dm(σ, pre)*dm(β, pre)*dm(It, pre)*dm(St, pre))/dm(N, pre) + dm(R_s(t, *args2), pre)
#               drdt = (dm(γ, pre) + dm(μ, pre))*dm(I, pre) - dm(μ, pre)*dm(R, pre) - dm(R_s(t, *args2), pre)
#        if t >= τ_σ and t < τ_ξ:
#                tindex = min(range(len(time)), key=lambda i: abs(time[i]- t)) - 1 # time.index(tt)
#                if tindex == len(Is):
#                      tindex -= 1
#                elif tindex >= len(Is):
#                       tindex = len(Is) - 1
#                It = Is[tindex]
#                St = Ss[tindex]
#                dsdt = dm(Λ, pre) - dm(μ, pre)*dm(S, pre) - (dm(β, pre)*dm(I, pre)*dm(S, pre))/dm(N, pre) + (dm(σ, pre)*dm(β, pre)*dm(It, pre)*dm(St, pre))/dm(N, pre) #Λ - μ*S - (β*I*S)/N + (σ*β*It*St)/N
#                drdt = dm((γ + μ), pre)*dm(I, pre) - dm(μ, pre)*dm(R, pre) #(γ + μ)*I - μ*R
#        elif t < τ_σ:
#               It = 0
#               St = 0
#               dsdt = dm(Λ, pre) - dm(μ, pre)*dm(S, pre) - (dm(β, pre)*dm(I, pre)*dm(S, pre))/dm(N, pre) + (dm(σ, pre)*dm(β, pre)*dm(It, pre)*dm(St, pre))/dm(N, pre) #Λ - μ*S - (β*I*S)/N  #+ R_s(t, *args2)
#               drdt = dm((γ + μ), pre)*dm(I, pre) - dm(μ, pre)*dm(R, pre) #(γ + μ)*I - μ*R
#        dedt = (dm(β, pre)*dm(I, pre)*dm(S, pre))/dm(N, pre) - (dm(σ, pre)*dm(β, pre)*dm(It, pre)*dm(St, pre))/dm(N, pre) - dm((μ + α), pre)*dm(E, pre)
#        didt = dm((μ + α), pre)*dm(E, pre) - dm((γ + μ), pre)*dm(I, pre) - dm(γ, pre)*dm(((dm(1, pre) - dm(κ_old, pre))*dm(N_old, pre) + (dm(1, pre) - dm(κ, pre))*(dm(1, pre) - dm(N_old, pre))), pre)*dm(I, pre)
#        dDdt = dm(γ, pre)*dm(((dm(1, pre) - dm(κ_old, pre))*dm(N_old, pre) + (dm(1, pre) - dm(κ, pre))*(dm(1, pre) - dm(N_old, pre))), pre)*dm(I, pre)

#        return [dm(dsdt, pre), dm(dedt, pre), dm(didt, pre), dm(drdt, pre), dm(dDdt, pre)]

def SEIRSM(tm, ym, *args):
       σm, βm, γm, αm, Λm, μm, ξm, κm, κ_oldm, τ_ξm, τ_σm, Nm, N_oldm, time, ts, Ism, Ssm, Rsm, pre = args
       args2 = [float(τ_ξm), float(ξm), time, Rsm]
       S = mv(ym[0], pre)
       E = mv(ym[1], pre)
       I = mv(ym[2], pre)
       R = mv(ym[3], pre)
       D = mv(ym[4], pre)
       M1 = mv(1.0, pre)
       tll = [abs(time[i] - tm)  for i in range(len(time))]
       if tm >= τ_σ and tm >= τ_ξ:
              tindex = min(range(len(time)), key=lambda i: tll[i]) - 1
              if tindex == len(Ism):
                     tindex -= 1
              elif tindex >= len(Ism):
                     tindex = len(Ism) - 1
              It = Ism[tindex]
              St = Ssm[tindex]
              dsdt = Λm - μm*S - (βm*I*S)/Nm + (σm*βm*It*St)/Nm + mv(R_s(tm, *args2), pre)
              drdt = (γm + μm)*I - μm*R - mv(R_s(tm, *args2), pre)
       if tm >= τ_σ and tm < τ_ξ:
               tindex = min(range(len(time)), key=lambda i: tll[i]) - 1
               if tindex == len(Ism):
                     tindex -= 1
               elif tindex >= len(Ism):
                      tindex = len(Ism) - 1
               It = Ism[tindex]
               St = Ssm[tindex]
               dsdt = Λm - μm*S - (βm*I*S)/Nm + (σm*βm*It*St)/Nm
               drdt = (γm + μm)*I - μm*R
       elif tm < τ_σ:
              It = mv(0.0, pre)
              St = mv(0.0, pre)
              dsdt = Λm - μm*S - (βm*I*S)/Nm + (σm*βm*It*St)/Nm
              drdt = (γm + μm)*I - μm*R
       dedt = (βm*I*S)/Nm - (σm*βm*It*St)/Nm - (μm + αm)*E
       didt = (μm + αm)*E - (γm + μm)*I - γm*((M1 - κ_oldm)*N_oldm + (M1 - κm)*(M1 - N_oldm))*I
       dDdt = γm*((M1 - κ_oldm)*N_oldm + (M1 - κm)*(M1 - N_oldm))*I
       return [dsdt, dedt, didt, drdt, dDdt]

def JacM(t, y, *args):
       σ, β, γ, α, Λ, μ, ξ, κ, κ_old, τ_ξ, τ_σ, N, N_old, time, ts, Is, Ss, Rs, pre  = args
       xx, yy = sp.symbols('xx yy')
       S = mv(y[0], pre)
       E = mv(y[1], pre)
       I = mv(y[2], pre)
       R = mv(y[3], pre)
       D = mv(y[4], pre)
       Dv = mv(0.0, pre) #mv(sp.integrate(DiracDelta(xx), (xx, 0.0, R - τ_ξ)), pre)
       M1 = mv(1.0, pre)
       jacc = [[I*β*σ/N - I*β/N - μ,      mv(0.0, pre),                                      S*β*σ/N - S*β/N,      ξ*Dv, mv(0.0, pre)],
               [   -I*β*σ/N + I*β/N, -α - μ,                                     -S*β*σ/N + S*β/N,                          mv(0.0, pre), mv(0.0, pre)],
               [                  mv(0.0, pre),  α + μ, -γ*(N_old*(M1 - κ_old) + (M1 - N_old)*(M1 - κ)) - γ - μ,                          mv(0.0, pre), mv(0.0, pre)],
               [                  mv(0.0, pre),      mv(0.0, pre),                                                γ + μ, -μ - ξ*Dv, mv(0.0, pre)],
               [                  mv(0.0, pre),      mv(0.0, pre),          γ*(N_old*(M1 - κ_old) + (M1 - N_old)*(M1 - κ)),                          mv(0.0, pre), mv(0.0, pre)]]
       return jacc
       # return [[I*β*σ/N - I*β/N - μ,      0.0,                                      S*β*σ/N - S*β/N,      ξ*Dv, 0.0],
       #         [   -I*β*σ/N + I*β/N, -α - μ,                                     -S*β*σ/N + S*β/N,                          0.0, 0.0],
       #         [                  0.0,  α + μ, -γ*(N_old*(1 - κ_old) + (1 - N_old)*(1 - κ)) - γ - μ,                          0.0, 0.0],
       #         [                  0.0,      0.0,                                                γ + μ, -μ - ξ*Dv, 0.0],
       #         [                  0.0,      0.0,          γ*(N_old*(1 - κ_old) + (1 - N_old)*(1 - κ)),                          0.0, 0.0]]

pib = [50]
Init_inf = 3
days = 1200
intval = 1000
daysi = int(days/10)
tint = days/intval
time_list = [float(i*tint) for i in range(intval+1)]
zhi_list = [0, 30, 90, 360]
τ_inc = 5.1
τ_rec = 18.8
R_0i = 5.2

for i in range(len(zhi_list)):
        σ = 0.0
        β = beta(R_0i, 1.0/18.8)
        γ = gamma(τ_rec)
        α = alpha(τ_inc)
        Λ = 0 # Birth rate
        μ = 0 # Death rate
        ξ = 0.01
        κ = 0.98
        κ_old = 0.96
        τ_ξ = zhi_list[i]
        τ_pre = 0
        τ_post = 0
        τ_σ = τ_σf(τ_pre, τ_post)
        N = 51.5 * 10**6
        N_old = 0.15
        S = [N]
        E = [20*Init_inf]
        I = [Init_inf]
        R = [0]
        D = [0]
        Sm = [mv(N, pib[0])]
        Em = [mv(20.0, pib[0])*mv(Init_inf, pib[0])]
        Im = [mv(float(Init_inf), pib[0])]
        Rm = [mv(0.0, pib[0])]
        Dm = [mv(0.0, pib[0])]
        for ii in tqdm(range(intval)):

              t_start = time_list[ii]
              t_end = time_list[ii+1]
              tevss = [ii/intval for ii in range(int(t_start)*intval, int(t_end)*intval, intval)]
              Y0 = [S[-1], E[-1], I[-1], R[-1], D[-1]]
              # Y0D = [dm(i, pi[0]) for i in Y0]
              Y0M = [Sm[-1], Em[-1], Im[-1], Rm[-1], Dm[-1]]
              α = alpha(τ_inc)
              γ = gamma(τ_rec)
              tevss = [t_start, t_end]
              argsl = [mv(σ, pib[0]), mv(β, pib[0]), mv(γ, pib[0]), mv(α, pib[0]), mv(Λ, pib[0]), mv(μ, pib[0]), mv(ξ, pib[0]), mv(κ, pib[0]), mv(κ_old, pib[0]), mv(τ_ξ, pib[0]), mv(τ_σ, pib[0]), mv(N, pib[0]), mv(N_old, pib[0]), time_list, [t_start, t_end], Im[:], Sm[:], Rm[:], pib[0]]
              argslf = [σ, β, γ, α, Λ, μ, ξ, κ, κ_old, τ_ξ, τ_σ, N, N_old, time_list, [t_start, t_end], Im[:], Sm[:], Rm[:], pib[0]]
              answer = solve_ivp(SEIRS, [float(t_start), float(t_end)], Y0,  t_eval=tevss,  method = 'Radau', args=(σ, β, γ, α, Λ, μ, ξ, κ, κ_old, τ_ξ, τ_σ, N, N_old, time_list, I[:], S[:], R[:]))
              # answer2 = solve_ivpd(SEIRSD, [t_start, t_end], Y0D, method = 'RadauD', args=(σ, β, γ, α, Λ, μ, ξ, κ, κ_old, τ_ξ, τ_σ, N, N_old, time_list, I[:], S[:], R[:], pi[0]))
              answer3 = solve_ivpm(SEIRSM, [t_start, t_end], Y0M, method='RadauM', t_eval=tevss, dense_output=False, prec=pib[0], jac=JacM, args=(argsl)) #mv(σ, pib[0]), mv(β, pib[0]), mv(γ, pib[0]), mv(α, pib[0]), mv(Λ, pib[0]), mv(μ, pib[0]), mv(ξ, pib[0]), mv(κ, pib[0]), mv(κ_old, pib[0]), mv(τ_ξ, pib[0]), mv(τ_σ, pib[0]), mv(N, pib[0]), mv(N_old, pib[0]), time_list, Im[:], Sm[:], Rm[:])) #  jac=JacM,  rtol=1E-1, atol=1E-2
              Sn = answer.y[0][-1]
              En = answer.y[1][-1]
              In = answer.y[2][-1]
              Rn = answer.y[3][-1]
              Dn = answer.y[4][-1]
              S.append(Sn)
              E.append(En)
              I.append(In)
              R.append(Rn)
              D.append(Dn)
              Snm = answer3.y[0][-1]
              Enm = answer3.y[1][-1]
              Inm = answer3.y[2][-1]
              Rnm = answer3.y[3][-1]
              Dnm = answer3.y[4][-1]
              Sm.append(Snm)
              Em.append(Enm)
              Im.append(Inm)
              Rm.append(Rnm)
              Dm.append(Dnm)

        Sp = [(i/N)*100.0 for i in S]
        Ep = [(i/N)*100.0 for i in E]
        Ip = [(i/N)*100.0 for i in I]
        Rp = [(i/N)*100.0 for i in R]
        Dp = [(i/N)*100.0 for i in D]

        Ip1 = (I).index(max(I))
        peakn = int(days*(Ip1/intval))
        Ip2 = (Ip).index(max(Ip))
        peakn2 = int(days*(Ip2/intval))
        Imax = max(I)*1.05
        fig = plt.figure()
        plt.plot(time_list, I, 'b-', label=r'$\it{Infected}$')
        plt.plot(time_list, D, '-', color='orange', label=r'$\it{Dead}$')
        plt.plot(time_list, Im, 'g--', label=r'$\it{Infected-M}$')
        plt.plot(time_list, Dm, 'k--', label=r'$\it{Dead-M}$')
        plt.legend([r'$\it{Infected}$', r'$\it{Dead}$', r'$\it{Infected-M}$', r'$\it{Dead-M}$'], loc="best", fontsize=15)
        plt.xlim((0,days))
        plt.ylim((0,Imax))
        plt.yticks([roundup(i*(Imax/10), intval) for i in range(11)])
        plt.xticks([int(i*daysi) for i in range(11)])
        plt.gca().get_yaxis().get_major_formatter().set_scientific(False)
        plt.gca().set_yticklabels([r'{:,}'.format(int(x)) for x in plt.gca().get_yticks()])
        plt.xlabel(r'$\bf{Time \ [Days]}$', fontsize=15)
        plt.ylabel(r'$\bf{Number \ of \ people}$', fontsize=15)
        plt.title(r'$\bf{SEIRS \ Method  \ for \ Spread \ of \ Disease}$', fontsize=18)
        plt.grid()
        fig.savefig(r"{}\SEIRS-{} Dead vs Infected.pdf".format(path_fol, τ_ξ), bbox_inches='tight')
        fig.savefig(r"{}\SEIRS-{} Dead vs Infected.svg".format(path_fol, τ_ξ), bbox_inches='tight')
        plt.show()
        plt.draw()
