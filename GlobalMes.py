import matplotlib.pyplot as mp
"""Sur un unique axe"""

def mise_en_forme():
    L=[]

    with open("//Users//maxencerichard//Documents//PSI*//TIPE//ArduinoSim//Tendeurs_Masse_1",'r') as file:
        txt=file.read()

    tmot=""
    for l in txt:
        if l not in ("\t","\n"):
            tmot+=l
        else:
            L.append(int(tmot))
            tmot=""
    return(L)

def tracer(L):
    abs=[50*i for i in range(len(L))]
    mp.plot(abs,L)
    mp.show()

def traitement(L):
    pass

def extremes(L):
    i=0
    while L[i]==L[0]:
        i+=1
        
    ex=[L[0],L[i]]
    y=[0,i]
    m= 0 if L[i]<L[0] else 1

    for j in range(i,len(L)):
        if m==0 and L[j]>L[j-1]:
            ex.append(L[j-1])
            y.append((j-1)/40)
            m=1
        elif m==1 and L[j]<L[j-1]:
            ex.append(L[j-1])
            y.append((j-1)/40)
            m=0


    return(y,ex)


def affichage(L,L2):
    abs=[50*i for i in range(len(L2))]
    mp.plot(abs,L2)

    mp.plot(L[1],L[0],".")
    mp.show()
