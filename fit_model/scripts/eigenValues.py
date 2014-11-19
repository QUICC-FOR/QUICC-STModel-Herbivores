from sympy import *
from scipy import linalg
import numpy as np
from sympy import init_printing

T, B, M, R = symbols('T B M R')
at, ab, bt, bb, tt, tb, e = symbols('at ab bt bb tt tb e')

# Empty state
R = 1 - T - B - M


# Differential equations describing the dynamics of the state variables
dTdt = at*(T+M)*(1-ab*(B+M))*R + tt*M - bb*(B+M)*T - e*T
dBdt = ab*(B+M)*(1-at*(T+M))*R + tb*M - bt*(T+M)*B - e*B
dMdt = at*(T+M)*ab*(B+M)*R + bb*(B+M)*T + bt*(T+M)*B - tt*M - tb*M -e*M

# Compute the Jacobian
systemEq = Matrix([dTdt, dBdt, dMdt])
stateVar = Matrix([T, B, M])

Jacob = systemEq.jacobian(stateVar)


# Solve the model for the case where T and M are at 0
# and B = 1 - e/ab (ie a l'équilibre si T=0 et M=0)
J_B = Jacob.subs([(T,0),(B,1-e/ab) , (M,0)])
eigs_J_B = J_B.eigenvals().keys()
invT = eigs_J_B

# Solve the model for the case where B and M are at 0
# and T = 1 - e/at (ie a l'équilibre si B=0 et M=0)
J_T = Jacob.subs([(T,1-e/at), (B,0), (M,0)])
eigs_J_T = J_T.eigenvals().keys()
invB = eigs_J_T
