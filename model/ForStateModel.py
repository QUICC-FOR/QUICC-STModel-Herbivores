# -*- coding: utf-8 -*-

#----------------------------------------------------------------------------
# THREE STATES MODEL
#----------------------------------------------------------------------------
from sympy import *

# state variables (G= grasslands, S=seedlings, T=trees)
G, S, T = symbols('G S T')

# parameters
# alpha = successsion rate (1/ageMaturity)
# delta = mortality rate (fire frequency * intensity + other disturbances)
# f = fecundity (growth rate)
# c = colonisation rate (recruitment in grasslands or empty patches)
alpha, delta, f, c = symbols('alpha delta f c')


# equations
G = 1-T-S
dT  = alpha*S - delta*T
dS  = f*T*c*(G+delta*T) - alpha*S

# resolution equilibre
equations = [Eq(dT, 0), Eq(dS, 0)]
equilibrium = solve(equations, T,S)
# one solution is null and do not interest us
Tstar = equilibrium[1][0]
Sstar = equilibrium[1][1]

Gstar = 1 - Tstar - Sstar
simplify(Gstar)

# one equilibrium state
systemEq = Matrix(equations)
stateVars = Matrix([T,S])
J = systemEq.jacobian(stateVars)

J

#subs
#Tstar.subs([(b, 0.1), (nu, 10), (a, 0.05), (c, 0)])


#----------------------------------------------------------------------------
# FOUR STATES MODEL
#----------------------------------------------------------------------------
from sympy import *

# state variables (G= grasslands, S=seedlings, C=boreal, D=temperate)
G, S, C, D = symbols('G S C D')

# parameters
# alpha = successsion rate (1/ageMaturity)
# delta = mortality rate (fire frequency * intensity + other disturbances)
# f = fecundity (growth rate)
# c = colonisation rate (recruitment in grasslands or empty patches)
alphac, alphad, deltac, deltad, fc, c, fd, k = symbols('alphac alphad deltac deltad fc fd c k')


# equations
G = 1-C-D-S
dC  = alphac*k*S - deltac*C
dD = alphad*(1-k)*S - deltad*D
dS  = fc*C*c*(G+deltac*C) + fd*D*c*(G+deltad*D) - (alphac*k + alphad*(1-k))*S


# Jacobian

systemEq = Matrix([dC,dD,dS])
stateVars = Matrix([C,D,S])
systemEq.jacobian(stateVars)


# resolution equilibre
equations = [Eq(dC, 0), Eq(dD, 0), Eq(dS, 0)]
#equilibrium = solve(equations, C, D,S)
## one solution is null and do not interest us
#Cstar = equilibrium[1][0]
#Dstar = equilibrium[1][1]
#Sstar = equilibrium[1][2]
#
#Gstar = 1 - Cstar - Dstar - Sstar
#simplify(Gstar)


# INVASIBILITY ANALYSIS

# Compute the Jacobian in absence of D and M
#
#f = systemEq
#v = stateVar
#J = jacobian(f,v)
#
#D = 0
#M = 0
#C = (aC-e)/aC
#
#J = [ aC*(aD*(D + M) - 1)*(C + D + M - 1) - e + aC*(aD*(D + M) - 1)*(C + M),           aC*(aD*(D + M) - 1)*(C + M) + aC*aD*(C + M)*(C + D + M - 1), sC + aC*(aD*(D + M) - 1)*(C + D + M - 1) + aC*(aD*(D + M) - 1)*(C + M) + aC*aD*(C + M)*(C + D + M - 1);           aD*(aC*(C + M) - 1)*(D + M) + aC*aD*(D + M)*(C + D + M - 1), aD*(aC*(C + M) - 1)*(C + D + M - 1) - e + aD*(aC*(C + M) - 1)*(D + M), sD + aD*(aC*(C + M) - 1)*(C + D + M - 1) + aD*(aC*(C + M) - 1)*(D + M) + aC*aD*(D + M)*(C + D + M - 1);               - aC*aD*(C + M)*(D + M) - aC*aD*(D + M)*(C + D + M - 1),               - aC*aD*(C + M)*(D + M) - aC*aD*(C + M)*(C + D + M - 1),  - e - sC - sD - aC*aD*(C + M)*(D + M) - aC*aD*(C + M)*(C + D + M - 1) - aC*aD*(D + M)*(C + D + M - 1)]
#
#l1 = e - aC
#
#l2 = (aD*e - (- 4*aC^2*aD*e*sC + aC^2*sC^2 + 2*aC^2*sC*sD + aC^2*sD^2 + 4*aC*aD*e^2*sC + 2*aC*aD*e*sC + 2*aC*aD*e*sD + aD^2*e^2)^(1/2))/(2*aC) - sC/2 - sD/2 - e
#
#l3 =  (aD*e + (- 4*aC^2*aD*e*sC + aC^2*sC^2 + 2*aC^2*sC*sD + aC^2*sD^2 + 4*aC*aD*e^2*sC + 2*aC*aD*e*sC + 2*aC*aD*e*sD + aD^2*e^2)^(1/2))/(2*aC) - sC/2 - sD/2 - e
#
#% THE TWO EIGENVALUES ARE THE SAME
#% NOW COMPUTE THE CRITICAL aD TO HAVE SPECIES D INVADING
#solve(l2,aD)
#aD_crit_D = (aC*(e + sC + sD))/(e + sD + sC*(1 - aC + e))
#
#
#
#% Compute the Jacobian in absence of C and M
#
#
#f = [dCdt;dDdt;dMdt]
#v = [C;D;M]
#J = jacobian(f,v)
#
#C = 0
#M = 0
#D = (aD-e)/aD
#
#J = [ aC*(aD*(D + M) - 1)*(C + D + M - 1) - e + aC*(aD*(D + M) - 1)*(C + M),           aC*(aD*(D + M) - 1)*(C + M) + aC*aD*(C + M)*(C + D + M - 1), sC + aC*(aD*(D + M) - 1)*(C + D + M - 1) + aC*(aD*(D + M) - 1)*(C + M) + aC*aD*(C + M)*(C + D + M - 1);           aD*(aC*(C + M) - 1)*(D + M) + aC*aD*(D + M)*(C + D + M - 1), aD*(aC*(C + M) - 1)*(C + D + M - 1) - e + aD*(aC*(C + M) - 1)*(D + M), sD + aD*(aC*(C + M) - 1)*(C + D + M - 1) + aD*(aC*(C + M) - 1)*(D + M) + aC*aD*(D + M)*(C + D + M - 1);               - aC*aD*(C + M)*(D + M) - aC*aD*(D + M)*(C + D + M - 1),               - aC*aD*(C + M)*(D + M) - aC*aD*(C + M)*(C + D + M - 1),  - e - sC - sD - aC*aD*(C + M)*(D + M) - aC*aD*(C + M)*(C + D + M - 1) - aC*aD*(D + M)*(C + D + M - 1)]
#
#
#l1 = e - aC
#
#l2 = -(2*aD*e - aC*e + aD*sC + aD*sD + (aC^2*e^2 - 4*aC*aD^2*e*sD + 4*aC*aD*e^2*sD + 2*aC*aD*e*sC + 2*aC*aD*e*sD + aD^2*sC^2 + 2*aD^2*sC*sD + aD^2*sD^2)^(1/2))/(2*aD)
#
#% THE TWO EIGENVALUES ARE THE SAME
#% NOW COMPUTE THE CRITICAL aD TO HAVE SPECIES D INVADING
#solve(l2,aD)
#aD_crit_C = aC*(e + sC + sD + e*sD)/(e + sC + sD*(1 + aC))
#
#
#
#% ---> MAIN RESULT: THERE ARE ALTERNATIVE STABLE STATES, BUT NO WAY TO GET COEXISTENCE. FOR SOME REGION OF THE PARAMETER SPACE WE HAVE THE FIRST SPECIES TO GET THERE TO WIN AND FOR OTHER REGIONS WE HAVE COMPETITIVE EXCLUSION

#----------------------------------------------------------------------------
#  Carrying capacity Herbivore Model
#----------------------------------------------------------------------------
from sympy import *
from scipy import linalg
import numpy as np

# carrying capacity
I, m, p, z, c = symbols('I m p z c')
K = (c*I - (p+q*exp(-z*c*I)))*(1-m) - m
Ik = solve(Eq(K, 0), I)
simplify(Ik)

theta, v, T, nu, H, A = symbols('theta v T nu H A')
K = theta*v*T/(nu*H + v*T) -A
fT = solve(Eq(K, 0), T)
simplify(fT)

#- fonction selection (gaussienne)
mu, sigma, phi, m = symbols('mu, sigma, phi, m')
eq1 = exp(-(1/2)*((mu)/sigma)**2) - phi
eq2 = exp(-(1/2)*((m - mu)/sigma)**2) - phi/2
res = solve([Eq(eq1,0),Eq(eq2,0)], mu, sigma)

#- fonction selection (logistique)
A, B, phi, p = symbols('A, B, phi, p')
eq1 = phi + (A-phi)/(1+B*exp(-p)) - 1
eq2 = phi + (A-phi)/(1+B*exp(0)) - (1-phi)/2
res = solve([Eq(eq1,0),Eq(eq2,0)], A, B)

simplify((phi - 1)*(phi + 1)/(-2*phi + (3*phi - 1)*exp(m) + 2) + 1 -phi)

#----------------------------------------------------------------------------
#  Herbivore coexistence
#----------------------------------------------------------------------------
from sympy import *
from sympy import init_printing
init_printing()
init_session()
from scipy.special import lambertw


thetab,thetag, taub, taug, nub, nug, mub, mug, u, v, w, A, B, T, S, G, qb, qg, pb, pg, mb, mg, zb, zg, rb, rg, k, Ikb, Ikg, Hg, Hb= symbols('thetab thetag taub taug nub nug mub mug u v w A B T S G qb qg pb pg mb mg zb zg rb rg k Ikb Ikg Hg Hb')


eqb = taub*k*u*S/(Hb*mub + k*u*S)  + (thetab*v*T/(Hb*nub + v*T))*( 1/ (1+ exp(rb*( (taub*k*u*S/(Hb*mub + k*u*S)) -pb-mb ))) ) -Ikb
Hbstar = roots(eqb)


B = (1+exp(rg*( taug*w*G/(mug*Hg+w*G) -pg-mg) ))*( Ikg - taug*w*G/(mug*Hg+w*G) )
eqg = - B *nug*Hg/ ((1-k)*u*(B-thetag))

iso = 1- G - S - eqb

coexistence = solve(Eq(iso,0), G)
symplify(coexistence)



#coexistence = solve([Eq(eqb, 0 ), Eq(eqg, 0), Eq(T + G + S , 1)], T, S, G)



