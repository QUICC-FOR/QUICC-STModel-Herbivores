# -*- coding: utf-8 -*-
"""
Created on Sat Apr 12 22:01:14 2014
resolution equation diff
@author: isabelle
"""

from __future__ import division
from scipy import *
from pylab import *
from scipy.integrate import odeint      # Module de résolution des équations différentielles

tau1 = 1
tau2 = 0.5

def deriv(syst, t):
    NBi = syst[0]                       # Variable1 NBi(t)
    NPo = syst[1]                       # Variable2 NPo(t)
    NPb = syst[2]                       # Variable3 NPb(t)
    dNBidt=-NBi/tau1                    # Equation différentielle 1
    dNPodt=NBi/tau1-NPo/tau2            # Equation différentielle 2
    dNPbdt=NPo/tau2                     # Equation différentielle 3
    return [dNBidt,dNPodt,dNPbdt]       # Dérivées des variables

# Paramètres d'intégration
start = 0
end = 5
numsteps = 50
t = linspace(start,end,numsteps)

# Conditions initiales et résolution
NBi0=1
NPo0=0
NPb0=0
syst_CI=array([NBi0,NPo0,NPb0])         # Tableau des CI
Sols=odeint(deriv,syst_CI,t)            # Résolution numérique des équations différentielles

# Récupération des solutions
NBi = Sols[:, 0]
NPo = Sols[:, 1]
NPb = Sols[:, 2]

# Graphiques des solutions
# Solutions analytiques
plot(t, exp(-t/tau1), '-b', lw=1.5, label=u"Sol.An. Bi")                                                                 
plot(t, (tau2*NBi0/(tau1-tau2))*(exp(-t/tau1)-exp(-t/tau2))+NPo0*exp(-t/tau2), '-r', lw=1.5, label=u"Sol.An. Po")         
plot(t, NBi0*(1-exp(-t/tau1)-(tau2/(tau1-tau2))*(exp(-t/tau1)-exp(-t/tau2))) + NPo0*(1-exp(-t/tau2)) + NPb0, '-g', lw=1.5, label=u"Sol.An. Pb")
# Solutions numériques
plot(t, NBi, 'o', ms=6, mfc='w', mec='b',label=u"Sol.Num. Bi")
plot(t, NPo, 'o', ms=6, mfc='w', mec='r',label=u"Sol.Num. Po")
plot(t, NPb, 'o', ms=6, mfc='w', mec='g',label=u"Sol.Num. Pb")


xlabel(ur"$t/\tau1 \, (\varnothing)$", fontsize=16)     # Label de l'axe des abscisses
ylabel(ur"$N/N_{0}(Bi) \, (\varnothing)$", fontsize=16) # Label de l'axe des ordonnées

legend()                                                # Appel de la légende
show()