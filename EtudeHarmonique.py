"""Optimisation Données"""

import numpy as np
import scipy.integrate as integr
import matplotlib.pyplot as plt

"""Constantes de travail"""

"""Fixées"""
m=0.580 #Masse Tour
h=0.372 #Frottement tiges
g=9.81 #Champ de pesenteur terrestre
k=39.46 #Constante de rappel du ressort équivalent

"""Paramétrables"""
w=7.57E-4 #Frottement pendule
L=0.24 #Longueur TMD
md=m/10 #Masse TMD

"""Excitation"""
x= lambda t: 0  #Excitation
xp= lambda t: 0 #Dérivée de l'exciation

##############
# Utilitaires
##############

def systeme_spp_opp(s,sp,o,op,t):
    """Renvoie les dérivées secondes à l'instant t"""
    #Pour les notations, se référer au document système matriciel
    
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
    return(o,op,s,sp)
    
    
################
# Méthode Scipy 
################

def f(Y,t):
    X=systeme_spp_opp(Y[2],Y[3],Y[0],Y[1],t)
    return(np.array([Y[1],X[0],Y[3],X[1]]))

def MethodeScipy(CI,duree):
    T = np.arange(0, duree, 10**-2)
    Y = integr.odeint(f, np.array(CI), T)
    #Acceleration:
    opp=[]
    spp=[]
    for i in range(len(Y[:,0])):
        X=systeme_spp_opp(Y[i,2],Y[i,3],Y[i,0],Y[i,1],T[i])
        opp.append(X[0])
        spp.append(X[1])
    
    return(max(spp))
                       
#############
# Resonance #
#############

def resonance():
    global x,xp
    essais= 10**np.linspace(0,2,100)
    reponses = []
    for j in range(len(essais)) :
        print(100-j)
        i=essais[j]
        x= lambda t: np.sin(i *t)*0.01
        xp= lambda t: i* np.cos(i *t)*0.01
        a=MethodeScipy(fCI(),10)
        reponses.append(a)
        
    plt.plot(essais,reponses,label="Valeur maximale atteinte par l'accélération du sommet de la tour" ,color='darkorange')
    
    plt.legend()
    plt.grid(True,which="both",linestyle='--')
    plt.xscale('log')
    plt.xlabel("Pulsation de l'excitation (rad.s\u207B\u00B9)",fontdict={'color':'darkred',
                  'weight': 'normal',
                    'size': 12})
    plt.ylabel("Accélération du sommet (m.s\u207B\u00B2)",fontdict={'color':'darkred',
                  'weight': 'normal',
                    'size': 12})
    plt.title("Système soumis à une excitation sinusoïdale d'amplitude\n1 cm et de pulsation variable",
              fontdict={'color':'darkred',
                  'weight': 'bold',
                    'size': 16})
    plt.show()
    print(essais[list(reponses).index(max(reponses))])
    
        
    
    
