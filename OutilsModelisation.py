"""Optimisation Données"""

import numpy as np
import scipy.integrate as integr
import matplotlib.pyplot as plt

"""Constantes de travail"""

"""Fixées"""

m=0.580    # Masse Tour
h=0.372    # Frottements tiges
g=9.81     # Champ de pesenteur terrestre
k=39.46    # Constante de rappel du ressort équivalent

"""Paramétrables"""
"""Valeurs de la maquette"""

w=7.57E-4  # Frottement pendule
L=0.24     # Longueur TMD
md=0.058   # Masse TMD

x= lambda t:0
xp=lambda t:0

###############
# Utilitaires #
###############

def systeme_spp_opp(Y,t):
    """Renvoie les dérivées secondes à l'instant t"""
    #Pour les notations, se référer au document système matriciel
    o,op,s,sp=Y
    
    #Matrice A
    Aa=md*L
    Ab=np.cos(o)*md
    Ac=md*L*np.cos(o)
    Ad=m+md
    detA=Aa*Ad-Ab*Ac

    assert detA!=0 , "Determinant de A nul !"
    detAi=1/detA #L'inverse de detA qui nous sera utile
    
    #Matrice B
    Ba=-md*g*np.sin(o)-w*op
    Bb=k*(x(t)-s)+h*(xp(t)-sp)+md*L*np.sin(o)*op**2-w*op*np.cos(o)
    
    X=[]
    X.append(detAi*(Ad*Ba-Ab*Bb))
    X.append(detAi*(-Ac*Ba+Aa*Bb))
    return(X)

def fCI(o=0,op=0,s=0,sp=0):
    """Renvoie les conditions initiales"""
    return(o,op,s,sp)
    
    
################
# Méthode Scipy 
################

def f(Y,t):
    X=systeme_spp_opp(Y,t)
    return(np.array([Y[1],X[0],Y[3],X[1]]))

def Methode_Scipy(CI,duree,precision):
    T = np.arange(0, duree, precision)
    Y = integr.odeint(f, np.array(CI), T)
    
    # On veut renvoyer aussi les courbes des dérivées secondes
    OPP=[]
    SPP=[]
    for i in range(len(Y[:])):
        X=systeme_spp_opp(Y[i],T[i])
        OPP.append(X[0])
        SPP.append(X[1])
    V=[Y[:,0],Y[:,1],OPP,Y[:,2],Y[:,3],SPP]
    return(T,V)

##############
# Affichages #
##############

def modification_para(xx=lambda t:0,xpp=lambda t:0,ww=7.57E-4,mdd=.058,LL=.24):
    global x,xp,w,md,L
    
    x,xp,w,md,L=xx,xpp,ww,mdd,LL
    
def aff(couleur,temps=15,extr=1):
    T,V=Methode_Scipy(fCI(),temps,1e-3)
    
    choix=[5]
    textelegende=["Theta","Dérivée de Theta","Dérivée seconde de Theta","S","Dérivée de S","Dérivée seconde de S"]
    for c in choix:
        plt.plot(T,V[c],label=textelegende[c],color=couleur)
        
        if extr:
            plt.plot([0,temps],[max(V[c]),max(V[c])],color=couleur)
            plt.plot([0,temps],[min(V[c]),min(V[c])],color=couleur)
        
    plt.legend()   
    