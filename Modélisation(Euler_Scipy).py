
"""Modélisation : Comparaison Euler/Scipy"""

import numpy as np
import scipy.integrate as integr
import matplotlib.pyplot as plt

"""Constantes de travail"""

g=9.81     # Champ de pesenteur terrestre

m=.58      # Masse Tour
k=39.46    # Constante de rappel du ressort équivalent au bâtiment
h=0.372    # Coefficient du frottement dans les tiges

md=0.058    # Masse TMD
L=.237     # Longueur pendule (TMD)
w=7.527E-4 # Frottement pendule : h sur les calucls de Nico


##############
# Utilitaires
##############

def systeme_spp_opp(Y):
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
    Bb=-k*s-h*sp+md*L*np.sin(o)*op**2-w*op*np.cos(o)
    
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
    X=systeme_spp_opp(Y)
    return(np.array([Y[1],X[0],Y[3],X[1]]))

def Methode_Scipy(CI,duree,dt):
    T = np.arange(0, duree, 10**-dt)
    Y = integr.odeint(f, np.array(CI), T)
    
    # On veut renvoyer aussi les courbes des dérivées secondes
    OPP=[]
    SPP=[]
    for i in range(len(Y[:])):
        X=systeme_spp_opp(Y[i])
        OPP.append(X[0])
        SPP.append(X[1])
    
    V=[Y[:,0],Y[:,1],OPP,Y[:,2],Y[:,3],SPP]
    return(T,V,dt," par Scipy ")
    
################
# Méthode Euler 
################

def Methode_Euler(CI,duree,dt):
    l=[0,1,3,4]#Les grandeurs auxquelles on applique Euler
    
    T = np.arange(0, duree, 10**-dt)
    V=[[CI[0]],[CI[1]],[],[CI[2]],[CI[3]],[]]#CI : notaion standard
    for t in T:
        Y=[V[0][-1],V[1][-1],V[3][-1],V[4][-1]]
        X=systeme_spp_opp(Y)#Le vecteur dérivée seconde
    
        V[2].append(X[0])#opp
        V[5].append(X[1])#spp
        #Mise en place Euler
        for i in l:
            V[i].append(V[i][-1]+V[i+1][-1]*10**-dt)
    for i in l:#Les grandeurs non dérivées seconde ont une valeur en trop
        V[i].pop()
    
    return(T,V,dt," par Euler ")

############
# Affichage 
############

def affichage_courbes(donnees,choix):
    """Choix des courbes : Notation standard"""
    T,V,precision,nom_methode=donnees
    
    
    textelegende=["Theta","Dérivée de Theta","Dérivée seconde de Theta","S","Dérivée de S","Dérivée seconde de S"]
    for c in choix:
        plt.plot(T, V[c],label=textelegende[c]+nom_methode+" dt = 1E-"+str(precision))
    plt.legend()

def comparaison(comp):
    #Comp contient le code de la méthode : 0 pour Euler, 1 pour Scipy
    # et la précision pour chaque métode
    
    for c in comp:
        if c[0]==0:
            affichage_courbes(Methode_Euler(fCI(s=0.1),20,c[1]),[5])
        else:
            affichage_courbes(Methode_Scipy(fCI(s=0.1),20,c[1]),[5])
    plt.title("Comparaison des méthodes d'Euler et Scipy")
    plt.show()
    