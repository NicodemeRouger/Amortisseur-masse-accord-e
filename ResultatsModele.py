import numpy as np
import matplotlib.pyplot as plt

from OutilsModelisation import *

######################################################

def MoyenneQuadratique(Y):
    tot=0
    for i in Y:
        tot+=i
    moy=tot/len(Y)
    
    totq=0
    for i in Y:
        totq+=(i-moy)**2
    return np.sqrt(totq/len(Y))


##########################
# Définition des entrées #
##########################

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

def EchelonPosition(couleur,ww=7.57E-4,mdd=.058,LL=0.24):
    modification_para(xx=xx,xpp=xpp,ww=ww,mdd=mdd,LL=LL)
    aff(couleur)

#EchelonPosition("indianred")
EchelonPosition("purple",mdd=20,LL=0.10424583,ww=10**-1.666666666666666)
EchelonPosition("lightblue",ww=10**-9)
EchelonPosition("lightgreen",mdd=.041745833,LL=0.5,ww=10**-1.666666666666666)
plt.show()

############################
# Optimisation frottements #
############################

def VariationWmin(mdd=0.058,LL=0.24):
   
    modification_para(xx,xpp,1E10,mdd,LL)
    mplus=min(Methode_Scipy(fCI(),5,1e-3)[1][5])
    
    varW=np.linspace(-2,2,100)
    deltaMax=[]
    for j in range(len(varW)):
        print(100-j)
        
        modification_para(xx,xpp,10**varW[j],mdd,LL)
        mmoins=min(Methode_Scipy(fCI(),5,1e-3)[1][5])
                   
        deltaMax.append(abs(mplus-mmoins))
    
    plt.plot(varW,deltaMax)
    print(varW[deltaMax.index(max(deltaMax))])    
    
    
        
def VariationWmax(mdd=0.058,LL=0.24):
   
    modification_para(xx,xpp,1E10,mdd,LL)
    mplus=max(Methode_Scipy(fCI(),5,1e-3)[1][5])
    
    varW=np.linspace(-2,2,100)
    deltaMax=[]
    for j in range(len(varW)):
        print(100-j)
        
        modification_para(xx,xpp,10**varW[j],mdd,LL)
        mmoins=max(Methode_Scipy(fCI(),5,1e-3)[1][5])
                   
        deltaMax.append(abs(mplus-mmoins))
          
    
    
    plt.plot(varW,deltaMax)

def VariationW_quad(mdd=0.058,LL=0.24):
    varW=np.linspace(-6,2,100)
    deltaMax=[]
    for j in range(len(varW)):
        print(100-j)
        
        modification_para(xx,xpp,10**varW[j],mdd,LL)
        deltaMax.append(MoyenneQuadratique(Methode_Scipy(fCI(),5,1e-3)[1][5]))
                   
    print(varW[deltaMax.index(min(deltaMax))])    
    
    plt.plot(varW,deltaMax)
    plt.show()
    
def VariationWext(mdd=0.058,LL=0.24):
    varW=np.linspace(-7,3,200)
    deltaMax=[]
    for j in range(len(varW)):
        print(100-j)
        
        modification_para(xx,xpp,10**varW[j],mdd,LL)
        Y=Methode_Scipy(fCI(),10,1e-2)[1][5]
        deltaMax.append(max(max(Y),abs(min(Y))))
                   
    print(varW[deltaMax.index(min(deltaMax))])    
    
    plt.plot(varW,deltaMax)
    plt.show()

######################################################
    
def VariationLext(mdd=0.058,ww=10**-1.3030303030303):
    varL=np.linspace(.2,.4,100)
    deltaMax=[]
    for j in range(len(varL)):
        
        
        modification_para(xx,xpp,ww,mdd,varL[j])
        Y=Methode_Scipy(fCI(),5,1e-3)[1][5]
        deltaMax.append(max(max(Y),abs(min(Y))))
                   
    print(varL[deltaMax.index(min(deltaMax))])    
    
    plt.plot(varL,deltaMax)
    plt.show()

def VariationLQuad(mdd=0.058,ww=10**-.636363636363):
    varL=np.linspace(0.1,0.5,100)
    deltaMax=[]
    for j in range(len(varL)):
        modification_para(xx,xpp,ww,mdd,varL[j])
        Y=Methode_Scipy(fCI(),5,1e-3)[1][5]
        deltaMax.append(MoyenneQuadratique(Y))
                   
    print(varL[deltaMax.index(min(deltaMax))])    
    
    plt.plot(varL,deltaMax)
    plt.show()
    
######################################################    

def VariationMext(LL=0.24,ww=7.57E-4):
    varM=np.linspace(.0001,1,100)
    deltaMax=[]
    for j in range(len(varM)):
        
        
        modification_para(xx,xpp,ww,varM[j],LL)
        Y=Methode_Scipy(fCI(),5,1e-3)[1][5]
        deltaMax.append(max(max(Y),abs(min(Y))))
                   
    print(varM[deltaMax.index(min(deltaMax))])    
    
    plt.plot(varM,deltaMax)
    plt.show()


def VariationMQuad(LL=0.29,ww=10**-1.3030303030303):
    varM=np.linspace(0.01,1,100)
    deltaMax=[]
    for j in range(len(varM)):
        modification_para(xx,xpp,ww,varM[j],LL)
        Y=Methode_Scipy(fCI(),5,1e-3)[1][5]
        deltaMax.append(MoyenneQuadratique(Y))
                   
    print(varM[deltaMax.index(min(deltaMax))])    
    
    plt.plot(varM,deltaMax)
    plt.show()

######################################################
    
def OptimisationExt(m_,l_,wmin,wmax,nb,tps):
    """w dans [10**wmin, 10**wmax]"""
    varW=np.linspace(wmin,wmax,nb)
    deltaMax=[]
    for j in range(nb):
        modification_para(xx,xpp,10**varW[j],m_,l_)
        Y=Methode_Scipy(fCI(),tps,1e-3)[1][5]
        deltaMax.append(max(abs(max(Y)),abs(min(Y))))
    plt.plot(varW,deltaMax)
    plt.show()
    a=min(deltaMax)
    return(varW[deltaMax.index(a)],a)
