#-*- coding: utf-8 -*-
import json
import numpy as np
# @Author: Soubès Franck

def lecture(fichier):
    f=open(fichier, "r")
    database=json.load(f)
    return(database)
    f.close()




    
param =  ['molecularweight','longueur_sequence','Sheet','Turn','Helix','aromaticity','hydrophobicity','cystein','phi'] #our choices 
precluster = lecture("Human.json")

def distance(entry_file,output_file):
    cpt_clst = 0
    file = open(entry_file,'r')
    
    file2= open(output_file,"w")
    file2.write("clusterisation of Human proteins\n")
    for line in file:
        print(line)
        taille =[]
        poids =[]
        pi = []
        helix =[]
        turn=[]
        sheet=[]
        hydro=[]
        cys =[]
        arom =[]
        if line !=('\n'):
            cpt_clst +=1
            entry = 'cluster numero:' + str(cpt_clst)+('\n')
            file2.write(entry)
            #print line
            line = line.strip('\n')
            tab=line.split(', ')
            
            
            print tab
            for i in range(len(tab)):
                
                taille.append(precluster[tab[i]]['longueur_sequence'])
                
                poids.append(precluster[tab[i]]['molecularweight'])
                pi.append(precluster[tab[i]]['phi'])
                helix.append(precluster[tab[i]]['Helix'])
                turn.append(precluster[tab[i]]['Turn'])
                sheet.append(precluster[tab[i]]['Sheet'])
                hydro.append(precluster[tab[i]]['hydrophobicity'])
                cys.append(precluster[tab[i]]['cystein'])
                arom.append(precluster[tab[i]]['aromaticity'])
            """
            print(taille)
            if (taille != []):
                try:
                    mean_taille= np.average(taille)
                except:
                    pass
            if (poids != []):
                try:
                    mean_poids = np.average(poids)
                except:
                    pass
            if (pi != []):
                try:
                    mean_pi= np.average(pi)
                except:
                    pass
            if (helix != []):
                try:
                    mean_helix = np.average(helix)
                except:
                    pass
            if (turn != []):
                try:
                    mean_turn = np.average(turn)
                except:
                    pass
            if (sheet != []):
                try:
                    mean_sheet=np.average(sheet)
                except:
                    pass
            if (hydro != []):
                try:
                    mean_hydro=np.average(hydro)
                except:
                    pass
            if (cys != []):
                try:
                    mean_cys = np.average(cys)
                except:
                    pass
            if (arom != []):
                try:
                    mean_arom = np.average(arom)
                except:
                    pass
                
            if (taille != []):
                try:
                    std_taille= np.std(taille)
                except:
                    pass
            if (poids != []):
                try:
                    std_poids = np.std(poids)
                except:
                    pass
            if (pi != []):
                try:
                    std_pi= np.std(pi)
                except:
                    pass
            if (helix != []):
                try:
                    std_helix = np.std(helix)
                except:
                    pass
            if (turn != []):
                try:
                    std_turn = np.std(turn)
                except:
                    pass
            if (sheet != []):
                try:
                    std_sheet=np.std(sheet)
                except:
                    pass
            if (hydro != []):
                try:
                    std_hydro=np.std(hydro)
                except:
                    pass
            if (cys != []):
                try:
                    std_cys = np.std(cys)
                except:
                    pass
            if (arom != []):
                try:
                    std_arom = np.std(arom)
                except:
                    pass
                """
            mean_taille = np.average(taille)
            #print(mean_taille)
            cpt_i =0
            #print(taille)
            mean_poids = np.average(poids)
            #print(mean_poids)
            mean_pi = np.average(pi)
            #print(mean_pi)
            mean_helix = np.average(helix)
            #print(mean_helix)
            mean_turn = np.average(turn)
            #print(mean_turn)
            mean_sheet=np.average(sheet)
            #print (mean_sheet)
            mean_hydro=np.average(hydro)
            #print(mean_hydro)
            mean_cys = np.average(cys)
            #print(mean_cys)
            mean_arom = np.average(arom)
            #print(mean_arom)
            std_taille =np.std(taille)
            std_poids = np.std(poids)
            std_pi = np.std(pi)
            std_helix = np.std(helix)
            std_turn = np.std(turn)
            std_sheet=np.std(sheet)
            std_hydro=np.std(hydro)
            std_cys = np.std(cys)
            std_arom = np.std(arom)
            file2.write('moyenne taille :')
	    file2.write(str(mean_taille))
            file2.write(" ")
	    file2.write('ecart-type taille :')
	    file2.write(str(std_taille))
	    file2.write('\n')
            file2.write('moyenne poids :')
	    file2.write(str(mean_poids))
            file2.write(" ")
	    file2.write('ecart-type poids :')
	    file2.write(str(std_poids))
	    file2.write('\n')
	    file2.write('moyenne point isoéléctrique :')
	    file2.write(str(mean_pi))
            file2.write(" ")
	    file2.write('ecart-type point isoéléctrique :')
	    file2.write(str(std_pi))
	    file2.write('\n')
	    file2.write('moyenne helix :')
	    file2.write(str(mean_helix))
            file2.write(" ")
	    file2.write(' ecart-type helix :')
	    file2.write(str(std_helix))
	    file2.write('\n')
	    file2.write('moyenne turn : ')
	    file2.write(str(mean_turn))
            file2.write(" ")
	    file2.write(' ecart-type turn : ')
	    file2.write(str(std_turn))
	    file2.write('\n')
	    file2.write('moyenne feuillet : ')
	    file2.write(str(mean_sheet))
            file2.write(" ")
	    file2.write(' ecart-type feuillet : ')
	    file2.write(str(std_sheet))
	    file2.write('\n')
	    file2.write('moyenne hydrophobicité :')
	    file2.write(str(mean_hydro))
            file2.write(" ")
	    file2.write(' ecart-type  hydrophobicité : ')
	    file2.write(str(std_hydro))
	    file2.write('\n')
	    file2.write('moyenne cysteines :')
	    file2.write(str(mean_cys))
            file2.write(" ")
	    file2.write(' ecart-type cysteines : ')
	    file2.write(str(std_cys))
	    file2.write('\n')
            file2.write('moyenne aromaticité :')
	    file2.write(str(mean_arom))
            file2.write(" ")
	    file2.write(' ecart-type aromaticité : ')
	    file2.write(str(std_arom))
	    file2.write('\n')
	    #hydro_std.append(std_hydro)
	    #histo_hydro.append(moyenne_hydro)
    file.close()
    file2.close()
    
            

distance("clustered.txt","res_intra.txt")


        




