import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from OutilsModelisation import *

##########
# Entrée #
##########

vm=.1
tm=.5
dt=.05
def xx(t):
    if t<= tm: 
        return(t**2*vm/(2*tm))
    elif t<= tm+dt:
        return(vm/dt*(t*(tm+dt-t/2)-tm/2*(tm+dt)))
    else:
        return xx(tm+dt)
def xpp(t):
    if t<= tm:
        return(vm/tm*t)
    elif t<= tm+dt:
        return(vm/dt*(tm+dt-t))
    else:
        return(0)
"""
xx=lambda t : 0.01*np.sin(8.2*t)
xpp=lambda t : 8.2*0.01*np.cos(8.2*t)
"""
#########################################
# Fonction d'optimisation du frottement #
#########################################

def MoyenneQuadratique(Y):
    """Renvoie la moyenne quadratique de la simulation"""
    tot=0
    for i in Y:
        tot+=i
    moy=tot/len(Y)
    
    totq=0
    for i in Y:
        totq+=(i-moy)**2
    return np.sqrt(totq/len(Y))

def OptimisationW(m_,l_,wmin,wmax,nbW,dureeScipy,precisionScipy,meth):
    """Renvoie l'accélération extreme ressentie minimale obtenue en faisant varier les frottements"""
    
    essais = np.linspace(wmin,wmax,nbW)
    listeAccExtreme = []
    
    for w_ in essais: 
        
        modification_para(xx=xx,xpp=xpp,ww=10**w_,mdd=m_,LL=l_)
        sim = Methode_Scipy(fCI(),dureeScipy,precisionScipy)[1][5]
        
        if meth:
            ma=abs(max(sim))
            mi=abs(min(sim))
            listeAccExtreme.append(max(ma,mi))
    
        else : listeAccExtreme.append(MoyenneQuadratique(sim))

    accExtremeMinimale = min( listeAccExtreme )
    p=listeAccExtreme.index(accExtremeMinimale)
    woptimal= essais[p]
    
    return (accExtremeMinimale,woptimal)


##########################
# Surface d'optimisation #
##########################      

""" Constantes """

nb  = 10
m_m = 0.06
l_m = .5

""""""

# Plan M,L :

Axes = Axes3D(plt.figure())

def Plan(Mm,Lm,largeurGrille):
    M_ = np.linspace(0.05, Mm, largeurGrille)
    L_ = np.linspace(1E-4, Lm, largeurGrille)
    return(M_,L_)
MM,LL=Plan(m_m,l_m,nb)
L,M = np.meshgrid(LL,MM)

# Cotes :

Z = np.zeros([nb,nb],float) 
W = np.zeros([nb,nb],float)

""" Constantes """

Cwmin            = -7
Cwmax            = 0.5
CnbW             = 10
CdureeScipy      = 2
CprecisionScipy  = 1e-2

methode          = True # True pour extreme, False pour quadratique
""""""
def Cotes():
    c=nb**2
    for i in range(nb):
        for j in range(nb):
            print(c)
            c-=1
            
            Lz,Lw=OptimisationW(MM[i],LL[j],Cwmin,Cwmax,CnbW,CdureeScipy,CprecisionScipy,methode)
            Z[i,j]=Lz
            W[i,j]=Lw
            
# Affichage 3D :

Cotes()

Axes.plot_surface(M,L,Z)
Axes.set_xlabel("Masse pendule (kg)",fontdict={'color':'darkred',
                  'weight': 'bold',
                    'size': 10})

Axes.set_ylabel("Longueur pendule (m)",fontdict={'color':'darkred',
                  'weight': 'bold',
                    'size': 10})
Axes.set_zlabel("Accélération du sommet de la tour (m.s\u207B\u00B2)",fontdict={'color':'darkred',
                  'weight': 'bold',
                    'size': 10})



plt.show()



    
    
    
    
    
    
    
    
    
    
    