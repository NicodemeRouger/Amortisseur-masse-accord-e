"""Sur un unique axe"""

def Arduino1axe(Serie):
    L=[]
    doc="/Users/maxencerichard/Documents/PSI*/TIPE/CapturesCoolTerm/Capture"+str(Serie)
    with open(doc,'r') as file:
        lignes=file.readlines()
    for i in range(3,len(lignes)-1):
        #La dernière donnée est souvent intraitable car tronquée lors de son écriture
        #Les trois lignes sont des lignes d'informations sur la série de données
        L.append(float((int(lignes[i][:-1])-330)*9.81/80))#Le dernier caractère est un retour à la ligne
    
    absc=[i/40 for i in range(len(L))]
        
    return(absc,L)

