import sys
from scipy.spatial import distance


##CONSTANTES
QQ=0.084
F=332

try:
    nom=sys.argv[1]
    chain=sys.argv[2]
except:
    print("Ce programme nécessite en premier argument le nom du fichier PDB à traiter, contenant les Hydrogènes, et en deuxième argument la lettre en majuscule de la chaine concernée.\n\nExemple : python3 FG_dssp.py 1btaFH.pdb A")
    sys.exit()
    
    
pdb=open(nom,"r")    #Fichier PDB comportant les H
coord=[]
attribut=[]                  #Rendu final vierge de structure
verif=[]


for line in pdb:            #Récupère les coordonnées de tout les N,C,O,H de la chaine peptidique et le numéro du résidu associé
    if line[17:20]=="PRO":
        continue
    elif line.startswith("ATOM") and line[21]==chain:
        if line[13:15]=="N ":
            coord.append([line[13:15],int(line[22:26]),float(line[30:38]),float(line[38:46]),float(line[46:54])])
            attribut.append([line[17:20],int(line[22:26]),"Coil"])
        elif line[13:15]=="C ":
            coord.append([line[13:15],int(line[22:26]),float(line[30:38]),float(line[38:46]),float(line[46:54])])
        elif line[13:15]=="O ":
            coord.append([line[13:15],int(line[22:26]),float(line[30:38]),float(line[38:46]),float(line[46:54])])
        elif line[13:15]=="H ":
            coord.append([line[13:15],int(line[22:26]),float(line[30:38]),float(line[38:46]),float(line[46:54])])
            verif.append(line[13:15])

if len(verif)==0:
    print("Votre fichier ne correspond pas aux critères d'utilisation.\n\nCe programme nécessite en premier argument le nom du fichier PDB à traiter, contenant les Hydrogènes, et en deuxième argument la lettre en majuscule de la chaine concernée.\n\nExemple : python3 FG_dssp.py 1btaFH.pdb A")
    sys.exit()
            
            
Hbond=[]                                #Liste des liaisons H

for i in range(3,len(coord)-3,4):       #Localisation et inclusion dans la liste de toutes les liaisons H au sein de la molécule
    for j in range(3,len(coord)-3,4):
        if i==j:
            continue
        else:
            nrj=QQ*((1/distance.euclidean(coord[i+2][2:5],coord[j][2:5]))
                    +(1/distance.euclidean(coord[i+1][2:5],coord[j+3][2:5]))
                    -(1/distance.euclidean(coord[i+2][2:5],coord[j+3][2:5]))
                    -(1/distance.euclidean(coord[i+1][2:5],coord[j][2:5])))*F
            if nrj < -0.5 and nrj > -3:
                Hbond.append([coord[i][1],coord[j][1],nrj])
           
t_Hbond=[]				# Sélection des liaisons les plus significatives
for i in range(0,len(Hbond)-1,1):
    if Hbond[i+1][0]==Hbond[i][0]:
        if abs(Hbond[i+1][0]-Hbond[i+1][1])>5:
            t_Hbond.append([Hbond[i+1][0],Hbond[i+1][1],Hbond[i+1][2]])
        elif abs(Hbond[i][0]-Hbond[i][1])>5:
            t_Hbond.append([Hbond[i][0],Hbond[i][1],Hbond[i][2]])
        elif Hbond[i+1][1]==Hbond[i+1][0]+4:
            t_Hbond.append([Hbond[i+1][0],Hbond[i+1][1],Hbond[i+1][2]])
        elif Hbond[i][1]==Hbond[i][0]+4:
            t_Hbond.append([Hbond[i][0],Hbond[i][1],Hbond[i][2]])
        elif Hbond[i+1][2]<Hbond[i][2]:
            t_Hbond.append([Hbond[i+1][0],Hbond[i+1][1],Hbond[i+1][2]])
        else:
            t_Hbond.append([Hbond[i][0],Hbond[i][1],Hbond[i][2]])
    elif Hbond[i-1][0]==Hbond[i][0]:
        continue
    else:
        t_Hbond.append([Hbond[i][0],Hbond[i][1],Hbond[i][2]])
        
T_Hbond=[]			# Suppression d'éventuels doublons identiques
for i in t_Hbond : 
    if i not in T_Hbond: 
        T_Hbond.append(i)     
                
Prevision=[]
for i in range(0,len(T_Hbond)-2,1):
    if T_Hbond[i][1]==T_Hbond[i][0]+4:                                            #Prévision des hélices, 4-Turn
        if T_Hbond[i+1][0]==T_Hbond[i][0]+1 and T_Hbond[i+1][1]==T_Hbond[i+1][0]+4:
            Prevision.append([T_Hbond[i][0],T_Hbond[i][1],T_Hbond[i][2],"Helix 4-T"])
        elif T_Hbond[i-1][0]==T_Hbond[i][0]-1 and T_Hbond[i-1][1]==T_Hbond[i-1][0]+4:
            Prevision.append([T_Hbond[i][0],T_Hbond[i][1],T_Hbond[i][2],"Helix 4-T"])
        else:
            Prevision.append([T_Hbond[i][0],T_Hbond[i][1],T_Hbond[i][2],"4-T"])
    elif T_Hbond[i][1]==T_Hbond[i][0]+3:                                              #Prévision des hélices, 3-Turn
        if T_Hbond[i+1][0]==T_Hbond[i][0]+1 and T_Hbond[i+1][1]==T_Hbond[i+1][0]+3:
            Prevision.append([T_Hbond[i][0],T_Hbond[i][1],T_Hbond[i][2],"Helix 3-T"])
        elif T_Hbond[i-1][0]==T_Hbond[i][0]-1 and T_Hbond[i-1][1]==T_Hbond[i-1][0]+3:
            Prevision.append([T_Hbond[i][0],T_Hbond[i][1],T_Hbond[i][2],"Helix 3-T"])
        else:
            Prevision.append([T_Hbond[i][0],T_Hbond[i][1],T_Hbond[i][2],"3-T"])
    elif T_Hbond[i][1]==T_Hbond[i][0]+5:                                            #Prévision des hélices, 5-Turn
        if T_Hbond[i+1][0]==T_Hbond[i][0]+1 and T_Hbond[i+1][1]==T_Hbond[i+1][0]+5:
            Prevision.append([T_Hbond[i][0],T_Hbond[i][1],T_Hbond[i][2],"Helix 5-T"])
        elif T_Hbond[i-1][0]==T_Hbond[i][0]-1 and T_Hbond[i-1][1]==T_Hbond[i-1][0]+5:
            Prevision.append([T_Hbond[i][0],T_Hbond[i][1],T_Hbond[i][2],"Helix 5-T"])
        else:
            Prevision.append([T_Hbond[i][0],T_Hbond[i][1],T_Hbond[i][2],"5-T"])
    else:
        if T_Hbond[i][0]+2==T_Hbond[i+2][0] and T_Hbond[i][1]+2==T_Hbond[i+2][1]:
            Prevision.append([T_Hbond[i][0],T_Hbond[i][1],T_Hbond[i][2],"B Para."])   #
        elif T_Hbond[i][0]-2==T_Hbond[i-2][0] and T_Hbond[i][1]-2==T_Hbond[i-2][1]:     # Prévision des bêtas parallèles
            Prevision.append([T_Hbond[i][0],T_Hbond[i][1],T_Hbond[i][2],"B Para."])   #
        elif T_Hbond[i][0]+2==T_Hbond[i+2][0] and T_Hbond[i][1]-2==T_Hbond[i+2][1]:
            Prevision.append([T_Hbond[i][0],T_Hbond[i][1],T_Hbond[i][2],"B AntiP."])  #
        elif T_Hbond[i][0]-2==T_Hbond[i-2][0] and T_Hbond[i][1]+2==T_Hbond[i-2][1]:     # Prévision des bêtas anti-parallèles
            Prevision.append([T_Hbond[i][0],T_Hbond[i][1],T_Hbond[i][2],"B AntiP."])  #
        else:
            Prevision.append([T_Hbond[i][0],T_Hbond[i][1],T_Hbond[i][2],"Turn"])      # Prévision d'un turn si aucune condition n'est rempli
            
for i in range(1,len(Prevision)-1,1):                                                 # Ajustement de prévision
    if Prevision[i-1][3]=="B Para." and Prevision[i+1][3]=="B Para.":
        Prevision[i][3]="B Para."
    elif Prevision[i-1][3]=="B AntiP." and Prevision[i+1][3]=="B AntiP.":
        Prevision[i][3]="B AntiP."

        
### Attribution générale des structures

j=0
for i in range(0,len(Prevision),1):        # Attribution
    while Prevision[i][0]!=attribut[j][1]:
        j=j+1
    if attribut[j][2]!="Coil":
        continue
    else:
        attribut[j][2]=Prevision[i][3] 

for i in range(1,len(attribut)-1,1):       # Ajustement de l'attribution 
    if attribut[i-1][2]=="B Para." and attribut[i+1][2]=="B Para.":
        attribut[i][2]="B Para."
    elif attribut[i-1][2]=="B AntiP." and attribut[i+1][2]=="B AntiP.":
        attribut[i][2]="B AntiP."
    elif attribut[i-1][2]=="Helix 4-T" and attribut[i+1][2]=="Helix 4-T":
        attribut[i][2]="Helix 4-T"

        
#### AFFICHAGE DES RÉSULTATS

sep="="*50
sep_="-"*50
real="Programme DSSP réalisé par :\nGASTRIN\t\tFrançois\nEtudiant en M2BI\nUniversité de Paris\n\nFichier PDB traité : {}\nChaine étudiée : {}\n\n".format(nom, chain) 
org="NOMENCLATURE :\n\nCoil : Random coil\n3-T/4-T/5-T : Turn-3, Turn-4, Turn-5\nHelix 3/4/5-T : Hélice composée de turn 3, 4 ou 5\nB Para. : Feuillet bêta Parallèle\nB AntiP. : Feuillet bêta Anti-Parallèle\n"
tabl="Residu\tNuméro\tStructure\n"
tabl_="_"*30
print(sep+"\n"+real+sep_+"\n\n"+org+sep+"\n\n"+tabl+tabl_+"\n")
for i in range(0,len(attribut),1):
    print("{}\t{}\t{}".format(attribut[i][0],attribut[i][1],attribut[i][2]))