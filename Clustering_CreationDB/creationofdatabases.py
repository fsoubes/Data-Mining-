# -*- coding: utf-8 -*-
# @Author: Soubès Franck, Debras Guillamaury, Frances Tristan
import json
import sys
import re
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import statistics
import numpy as np

def lecture(fichier):
    f=open(fichier, "r")
    data=json.load(f)
    return(data)
    f.close()

human =  "/home/franck1337/DEA/Homo_Sapiens/MERGED.json"
asp =   "/home/franck1337/DEA/Aspergillus/MERGED.json"

lecHum=lecture(human)

def save(dico, nom):
    fic=open(nom, "w")
    json.dump(dico,fic)
    fic.close()

def createtable(fichier):
    
    
    database = {}

    kd = { 'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2, 'X': 0, 'U': 0, 'B': 0, 'Z': 0}

    aa_mol = {"A":89.1,"R":174.2,"N":132.1,"D":133.1,"C":121.2,"E":147.1,"Q":146.2,"G":75.1,"H":155.2,"I":131.2,"L":131.2,"K":146.2,"M":149.2,"F":165.2,"P":115.1,"S":105.1,"T":119.1,"W":204.2,"Y":181.2,"V":117.1,"X":120.1,"U":120.1,"B":120.1,"Z":120.1}

    numbers = [aa_mol[key] for key in aa_mol]
    meanweight = statistics.mean(numbers)
    molecular_weight =[]
    longueur_sequence =[]
    phi = []
    aromaticity=[]
    hydrophobicity=[]
    Turn =[]
    Helix=[]
    Sheet=[]


    fl = open(fichier, 'r')
    print(fl)

    cptBeta = 0;
    line = fl.readline()
    
    for key in lecHum.keys():
        if "fasta" in lecHum[key]:
            for seq in lecHum[key]["fasta"]:
                if "Beta strand" in seq["type_seq"]:
                    cptBeta+=1

 
    parsed_json = json.loads(line)
    

    
    for k in parsed_json:
        sum_aa = 0.0
        sum_kd = 0.0
        cpt_cys = 0.0
        my_seq = parsed_json[k]["seq"]
        
        for aa in parsed_json[k]["seq"]:
            if aa in  aa_mol.keys():
                sum_aa += aa_mol[aa]
                sum_kd += kd[aa]
    
    
        analysed_seq = ProteinAnalysis(my_seq)
        try :
            molecular_weight.append(sum_aa)
        except :
            pass
    
        try :
            longueur_sequence.append(len(my_seq))
        except :
            pass                    
        try :
            phi.append(analysed_seq.isoelectric_point())
        except :
            pass
        try :
            aromaticity.append( analysed_seq.aromaticity())
        except :
            pass
        try :
            
            hydrophobicity.append(sum_kd/len(parsed_json[k]["seq"]))
        except :
            pass
        
        try :
            secondary = analysed_seq.secondary_structure_fraction()
            Helix.append(secondary[0])
            Turn.append(secondary[1])
            Sheet.append(secondary[2])
            
        except :
            pass

    meanw = np.mean(molecular_weight)
    print(meanw)
    meanl = np.mean(longueur_sequence)
    print(meanl)
    meanpi= np.mean(phi)
    print(meanpi)
    meanar = np.mean(aromaticity)
    print(meanar)
    meanhy = np.mean(hydrophobicity)
    print(meanhy)
    meanhe = np.mean(Helix)
    print(meanhe)
    meantu = np.mean(Turn)
    print(meantu)
    meansh = np.mean(Sheet)
    print(meansh)


    fl = open(fichier, 'r')
    print(fl)

    cptBeta = 0;
    line = fl.readline()
    
    for key in lecHum.keys():
        if "fasta" in lecHum[key]:
            for seq in lecHum[key]["fasta"]:
                if "Beta strand" in seq["type_seq"]:
                    cptBeta+=1

 
    parsed_json = json.loads(line)

    for k in parsed_json:
        um_aa = 0.0
        sum_kd = 0.0
        cpt_cys = 0.0
        database[k] = {}
        
        my_seq = parsed_json[k]["seq"]
        database[k]["seq"] = my_seq

    
        for aa in parsed_json[k]["seq"]:
            if aa in  aa_mol.keys():
                sum_aa += aa_mol[aa]
                sum_kd += kd[aa]

        for aa in parsed_json[k]["seq"]:
            if aa == 'C':
                cpt_cys +=1
           
        analysed_seq = ProteinAnalysis(my_seq)
        try :
            database[k]["molecularweight"] = (sum_aa)
        except :
            database[k]["molecularweight"] = meanw
            print "weigth replace by the mean"
        try :
            database[k]["longueur_sequence"] = len(my_seq)
        except :
            database[k]["longueur_sequence"] = meanl
            print "size replace by mean"
        try :
            database[k]["phi"] = analysed_seq.isoelectric_point()
        except :
            database[k]["phi"] = meanpi
            print "no phi"
        try :
            database[k]["aromaticity"] = analysed_seq.aromaticity()
        except :
            database[k]["aromaticity"] = meanar
            print "pas de aromaticity"
        try :
            database[k]["hydrophobicity"]= sum_kd/len(parsed_json[k]["seq"])
        except :
            database[k]["hydrophobicity"]= meanhy
            print "pas de hydro"

        try :
            secondary = analysed_seq.secondary_structure_fraction()
            database[k]["Helix"] = secondary[0]
            database[k]["Turn"] = secondary[1]
            database[k]["Sheet"] = secondary[2]

        except :
            database[k]["Helix"] = meanhe
            database[k]["Turn"] = meantu
            database[k]["Sheet"]=meansh
            print "no secondary structures"
        try:
            database[k]["cystein"] = cpt_cys
        except :
            database[k]["cystein"] = 0
            print "no cys"

    print("end")
    return database
    
    


def creation_csv(database):
    
    print("csv creation")
    
    csv =open("realtable.csv", "w")
    csv.write("uniprot_id;length;phi;weight;aromaticity;hydrophobicity;sheet;helix;turn\n")


    for key in database.keys():
        uniprot_id = key
        length = database[key]["longueur_sequence"]
        phi = database[key]["phi"]
        weight = database[key]["molecularweight"]
        aromaticity = database[key]["aromaticity"]
        hydrophobicity = database[key]["hydrophobicity"]
        sheet  = database[key]["Sheet"]
        helix = database[key]["Helix"]
        turn = database[key]["Turn"]
		
	cystein = database[key]["cystein"]
        row = str(uniprot_id) +   ";"+ str(length) + ";" + str(phi) + ";" + str(weight) + ";" + str(aromaticity) + ";" + str(hydrophobicity) + ";" + str(sheet) + ";" + str(helix) + ";" + str(turn) + ";" +str(cystein) +"\n"
        csv.write(row)
    print("end of csv")


structure = createtable(human)

structure = createtable(asp)
save(structure, "Human1.json")

creation_csv(structure)

        
"""

def nbProt(listeProt):
    return (len(listeProt))

def nbBeta(listeProt):
    cptBeta=0
    for value in listeProt:
        cptBeta+=value["Beta strand"]
    return (cptBeta)

def nbHelix(listeProt):
    cptHelix=0
    for value in listeProt:
        cptHelix+=value["Helix"]
    return (cptHelix)

def nbTurn(listeProt):
    cptTurn=0
    for value in listeProt:
        cptTurn+=value["Turn"]
    return (cptTurn)


def percentBeta(structure, nbBeta):
    return (nbBeta/structure*100)

def percentHelix(structure, nbHelix):
    return (nbHelix/structure*100)

def percentTurn(structure, nbTurn):
    return (nbTurn/structure*100)

def getMinMaxMass(protein_list):

    mass = 2 #index of mass in a protein list
    massMax = 0
    massMin = 1000 #arbitrary choosen
    for protein in protein_list:
        if not protein[mass] is None:
            if protein[mass] > massMax:
                massMax = protein[mass]
            if protein[mass] < massMin:
                massMin = protein[mass]
    return massMin, massMax



nbStructureAsp = nbStructure(lecAsp)
print(nbStructureAsp)
nbProtAsp = nbProt(listeAsp)

nbBetaAsp = nbBeta(listeAsp)
percentBetaAsp = percentBeta(nbStructureAsp, nbBetaAsp)

"""
