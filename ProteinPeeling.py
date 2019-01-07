

#------------------MASTER PROJECT : PROTEIN PEELING ---------------------- #
#-------------------- TIEO SONIA - M2BI -----------------------------------#

__author__ = "Sonia Tieo"

#*********************IMPORTATIONS *********************#

import os
import sys
import glob
import Bio.PDB
import numpy as np
import pylab
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import subprocess
import time

#*********************GLOBAL VARIABLES*********************#

D0 = 8.0 # CUTOFF dist< 8 A : no interaction
DELTA = 1.5 #logistic function parameter
CWD = os.getcwd()
#CWD = os.path.abspath(os.path.join(os.getcwd(), os.pardir))

#*********************FONCTIONS*********************#

def calc_residue_dist(residue_one, residue_two):
    """ C-alpha distance between two residues.
        residue_one, residue_two: <class 'Bio.PDB.Residue.Residue'>
        Returns: int C-alpha distance
    """
    diff_vector = residue_one["CA"].coord - residue_two["CA"].coord
    return np.sqrt(np.sum(diff_vector**2))


def ss_dict(model, pdb_filename) :
    """Secondary structure assigned by DSSP.
    model:<class 'Bio.PDB.Model.Model'>
    pdb_filename:str path and pdb file
    Returns:dict structure of each residue
    """
    dssp = Bio.PDB.DSSP(model, pdb_filename)
    dict_ss = {}
    for k in dssp.keys():
        if k[0] == CHAINE:
            dict_ss[(k[1][1])] = dssp[k][2]
    return dict_ss


def calc_dist_matrix(chain_one, chain_two, dict_ss) :
    """Matrix of C-alpha distances between two same chain
    chain_one,chain_two:<class 'Bio.PDB.Chain.Chain'>
    dict_ss:dict structure of each residue
    Returns:numpy.ndarray Distance Matrix
    """
    dist = np.zeros((len(chain_one), len(chain_two)), np.float)
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            if ('CA' in residue_one.child_dict and 'CA' in residue_two.child_dict):
                dist[row, col] = calc_residue_dist(residue_one, residue_two)
            else:
                dist[row, col] = 0
    # Remove (HETATOM or Res not assigned with DSSP)
    L=[]
    for residue_one in list(chain_one.get_residues()):
        L.append(residue_one.id[1])
    diff = list (set(L) - set(list(dict_ss)))
    new_matrix = np.delete(dist,diff,axis = 0)
    new_matrix2 = np.delete(new_matrix,diff,axis = 1)
    return new_matrix2


def dist_to_proba(d) :
    """Transform a distance to contact probability (logistic transformation)
    d:int
    Returns:int probability
    """
    return(1/( 1+np.exp((d-D0)/DELTA)))


def is_in_ss(a, dict_ss):
    """Return True if amino acid is in Secondary Structure"""
    # T,C and S are not secondary structure
    return (dict_ss[a] != 'T' and dict_ss[a] != '-' and
            dict_ss[a] != 'C' and dict_ss[a] != 'S')



def calc_PI_ab(contact_map, a, b, dict_ss):
    """Calcul PI for given residue number a and d
    contact_map: numpy.ndarray
    a,b : int residues number
    dict_ss: dict of secondary structure
    Returns: float PI value
    """
    #To not cut a protein in a secondary structure
    if (is_in_ss(a, dict_ss) or is_in_ss(b, dict_ss)):
        PI = 0
    else :
        start = 0
        end = contact_map.shape[1]
        cm_df = pd.DataFrame(data=contact_map)
        cm_df.index = np.arange(1, len(cm_df)+1)
        #Sum all probabilities in A, B, C blocks
        A = cm_df.iloc[a:b , a:b].values.sum()
        B = cm_df.iloc[b:end , b:end].values.sum()
        C = cm_df.iloc[b:end, a:b].values.sum()
        PI = (A*B - C*C)/((A+B)*(B+C))
    return PI


def all_PI_in_dict(contact_map ,dssp_dict, min_size_PU = 10 , max_size_PU = 50):
    """ Stock all start residus a with liste of PI values
    contact_map: numpy.ndarray
    dssp_dict: dict of secondary structure
    min_size_PU, max_size_PU : int
    Returns: Dict with start res a as key and a list of PI for each end res
    """
    dict_PI = {}
    for a in range(1,len(contact_map) -  max_size_PU  - 1):
        liste_PI =[]
        for m in range(a , a + max_size_PU):
            liste_PI.append(calc_PI_ab(contact_map, a, m, dssp_dict))
            dict_PI[a] = liste_PI
    #remove null list
    for i in range(1, len(dict_PI)):
        if sum(dict_PI[i]) == 0:
            dict_PI.pop(i, None)
    return dict_PI



def save_plot_allPI (dict_PI, pdb_code, max_size_PU):
    """Plot and save Pi values for each starting res
    dict_PI : dict all PUs values
    pdb_code : str
    max_size_PU : max size for a PU
    """
    for key,value in dict_PI.items():
        fig, ax = plt.subplots()
        ax.plot(list(range(key,key+max_size_PU)) , value )
        ax.set(xlabel='PI', ylabel='residues',
               title= key)
        ax.grid()
        dir = CWD + "/results/" + pdb_code + "/all_PI/"
        if not os.path.exists(dir):
            os.makedirs(dir)
        fig.savefig( dir  + "a-" + str(key) + "-PI.pdf")


def new_max_PI_in_dict(dict_PI):
    """Select unique max PI for each starting res
    Returns: dict with starting res as key and tuple (ending res,) as value
    """
    dict_PI_max={}
    for key in dict_PI.keys():
        dict_PI_max[key] = (key + dict_PI[key].index(max(dict_PI[key])), max(dict_PI[key]))
    return dict_PI_max


def new_combis_PU(dict_PI_max, min_size_PU= 10, max_size_PU = 50):
    """Make combination of PU with new_max_PI_in_dict
    dict_PI_max: dict with starting res as key and tuple (ending res,) as value
    Returns: list with all combination of splitting """
    all_combis_PU = []
    for k in dict_PI_max.keys():
        Ltmp=[k]
        while k in  dict_PI_max.keys() and k!= list(dict_PI_max.keys())[-1] :
            Ltmp.append(dict_PI_max[k][0])
            k = dict_PI_max[k][0]
        all_combis_PU.append(Ltmp)
    #Nettoyage
    new_all_combis_PU =[]
    for el in all_combis_PU:
        if ((len(el) != 1) and (el[0] >= min_size_PU or el[0] == 1) and
            (el[0] < max_size_PU) and (el[-1] <= (len(dssp_dict)-min_size_PU))):
            new_all_combis_PU.append(el)
    return new_all_combis_PU


#----------------------------------------------/!\  amelioration: not used here
def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def maxs_PI_in_dict(dict_PI):
    """ Select all max PI in range of 5 for each start"""
    dict_PI_max={}
    for key in dict_PI.keys():
        maxPI = max(l for l in list(chunks(dict_PI[key], 5)))
        maxPI = list(filter(lambda a: a != 0, maxPI))
        dict_PI_max[key]=[]
        for val,l in zip(dict_PI[key],range(len(dict_PI[key]))):
            if val in maxPI:
                dict_PI_max[key].append((l+key,val))
    return dict_PI_max


def test_next_PU(dict_PI_max,k):
    for i in range(len(dict_PI_max[k])):
        if dict_PI_max[k][i][0] in dict_PI_max.keys():
            return dict_PI_max[k][i][0]


def combis_PU(dict_PI_max):
    """Combi of PU with maxs_PI_in_dict"""
    #Choix combinaison PU ,attention
    #boucle while on prend juste en fait le premier PI potable
    # Critères à changer
    all_combis_PU = []
    for k in dict_PI_max.keys():
        Ltmp=[k]
        while test_next_PU(dict_PI_max,k) in  dict_PI_max.keys():
            Ltmp.append(test_next_PU(dict_PI_max,k))
            k = test_next_PU(dict_PI_max,k)
        all_combis_PU.append(Ltmp)
    #Nettoyage
    new_all_combis_PU =[]
    for el in all_combis_PU:
        if ((len(el) != 1) and (el[0] >= 10 or el[0] == 1) and
            (el[0] < 50) and (el[-1] <= (len(dssp_dict)-10) ) ):
            new_all_combis_PU.append(el)
    return new_all_combis_PU

#---------------------------------------------- end of  "amelioration not used here"


def calc_sigma_decoup(PUs):
    """ Calculate sigma for a given all PUs in given list
    PUs: list of res number that delimited a PU
    Returns: list sigma values
    """
    # sigma:independance
    #pLus le sigma est petit mieux c'est
    sigma_tmp=[]
    for i in range(0,len(PUs) - 2 ):
        alpha=0.43
        a = PUs[i]
        b = PUs[i+1]
        c = PUs[i+2]
        cm_df = pd.DataFrame(data=contact_map)
        PU_A = cm_df.iloc[a:b , a:b]
        PU_C = cm_df.iloc[b:c , a:b]
        PU_B = cm_df.iloc[b:c , b:c]
        numer = (PU_C.values.sum()*2)/((PU_A.shape[1]**alpha)*(PU_C.shape[1]**alpha))
        denom = (PU_A.values.sum()+PU_B.values.sum()+PU_C.values.sum())/(PU_A.shape[1]+PU_B.shape[1]+PU_C.shape[1])
        sigma = numer/denom
        sigma_tmp.append(sigma)
    return(sigma_tmp)

def calc_kappa_decoup(PUs):
    """ Calculate sigma for a given all PUs in given list
    PUs: list of res number that delimited a PU
    Returns: list kappa values
    """
    #kappa compacité
    #pLus kappa est grand mieux c'est
    kappa_tmp=[]
    for i in range(0,len(PUs) - 1):
        a =PUs[i]
        b = PUs[i+1]
        cm_df = pd.DataFrame(data=contact_map)
        PU = cm_df.iloc[a:b , a:b]
        kappa = PU.values.sum() / PU.shape[1]
        kappa_tmp.append(kappa)
    return kappa_tmp

def make_liste_PI(dict_PI_max, PUs):
    """List PI for all PUs in combination
    dict_PI_max : dict with starting res as key and tuple (ending res,PI) as value
    PUs : list of PUs
    Returns: list of PUs
    """
    liste_PI = []
    liste_PI_idx=[]
    for PU in PUs:
        if PU in dict_PI_max.keys():
            liste_PI.append(dict_PI_max[PU][1])
            liste_PI_idx.append(PU)
    return([liste_PI,liste_PI_idx])


def add_start_end_PUs(all_PUs, contact_map):
    """Add first and end of residus in list of all Pus
    all_PUs: list of all PUs combination
    contact_map: numpy.ndarray
    Returns new list of all PUs combination
    """
    all_PUs_bis=[]
    for PUs in all_PUs:
        if PUs[0] != 1:
            PUs = [1] + PUs
        if PUs[-1] != contact_map.shape[1]:
            PUs = PUs + [contact_map.shape[1]]
        all_PUs_bis.append(PUs)
    return all_PUs_bis


def best_split(all_PUs_sigma, all_PUs_kappa):
    """Give index of best protein splitting and write file with ranking splitting
    all_PUs_sigma : list of list of sigmas of all PUs combination
    all_PUs_kappa : list of list of kappas of all PUs combination
    Returns : int index of the best combination
    """
    #SUm all sigma,les kappas
    #SOrt simgas, kappa -> Give the ranking
    # Select the lowest ranking
    sigma_sum=[]
    kappa_sum=[]
    for i in range(0, len(all_PUs)):
        sigma_sum.append(sum(all_PUs_sigma[i]) / len(all_PUs_sigma[i]))
        kappa_sum.append(sum(all_PUs_kappa[i]) / len(all_PUs_kappa[i]))
    sigma_sum_sorted = [sorted(sigma_sum).index(v) for v in sigma_sum]
    kappa_sum_sorted = [sorted(kappa_sum, reverse=True).index(v) for v in kappa_sum]
    sigma_kappa_sum = [sum(x) for x in zip(sigma_sum_sorted, kappa_sum_sorted)]
    sk_rank = [sorted(sigma_kappa_sum).index(v) for v in sigma_kappa_sum]
    with open(CWD + "/results/" + pdb_code + "/ranking_spitting.txt", "w") as f:
        for sk in range(len(sigma_kappa_sum)):
            f.write( "PU:" + str(sk + 1) )
            f.write(" - Rank: " + str(sk_rank[sk] + 1))
            f.write(" - Split: " + str(all_PUs[sk]))
            f.write(" - Sigma Moy: " + str(round(sigma_sum[sk], 4)))
            f.write(" - Sigma Std: " + str(round( np.std(all_PUs_sigma[sk]), 4)))
            f.write(" - Kappa Moy: " + str((round(kappa_sum[sk], 4))))
            f.write(" - Kappa Std: " + str(round(np.std(all_PUs_kappa[sk]),4)) + "\n")
    return(sigma_kappa_sum.index(min(sigma_kappa_sum)))


def write_pdb_all_split(all_PUs):
    """Split PDB per PUs and write PDB for each PU
    all_PUs : all PUs combination
    """
    for PU in range(len(all_PUs)):
        choix = all_PUs[PU]
        for c in range(1,len(choix)):
            start_res= choix[c-1]
            end_res= choix[c]
            chain_id = CHAINE
            dir = CWD + "/results/" + pdb_code + "/PU" + str(PU + 1)
            file = dir + "/PU" + str(PU + 1) + "." + str(start_res) + "-" + str(end_res) +".pdb"
            if not os.path.exists(dir):
                os.makedirs(dir)
            io = Bio.PDB.PDBIO()
            io.set_structure(structure[0][CHAINE])
            io.save(file, ResSelect())


class ResSelect(Bio.PDB.Select):
    """Cut and Save PDB between start res and end res"""
    ##credits : https://stackoverflow.com/questions/22452748/remove-parts-from-a-pdb-file-using-python
    def accept_residue(self, res):
        if res.id[1] >= start_res and res.id[1] <= end_res:
            return True
        else:
            return False


def plot_contact_map(choix):
    """Plot contact map for the best split
    choix : list of PUs of best the best split
    """
    color=iter(plt.cm.rainbow(np.linspace(0,1,len(choix))))
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    plt.imshow(contact_map, cmap='gray_r', interpolation='None')
    for x in range(1, len(choix)):
        c=next(color)
        if  x != len(choix) - 1:
            plt.axvline(x=choix[x],color='lavender', linewidth=1,)
            plt.axhline(y=choix[x],color='lavender', linewidth=1,)
        ax1.add_patch(patches.Rectangle(
        (choix[x-1], choix[x-1]), # (x,y)
        choix[x] - choix[x-1] , # width
        choix[x] - choix[x-1] , # height
    # You can add rotation as well with 'angle'
    alpha=0.2, facecolor=c, edgecolor="black", linewidth=1, linestyle='solid'
    )
    )
        plt.text(choix[x-1] + (choix[x]-choix[x-1])/2.3,
                choix[x-1]+(choix[x]-choix[x-1])/2.5, x, fontweight = 'book')
    fig1.savefig("{0}/results/{1}/contactmap{2}.pdf".format(CWD, pdb_code,str(choix)))




def plot_criteria(PU):
    """PLot sigma, kappa, PI of PUs for each combination
    PU : list of a PUs combination in all PUs combinations
    """
    f, (PIs, ks) = plt.subplots(2, 1)
    PIs.plot(all_PUs_PI_index[PU], all_PUs_PI[PU], 'o-', label='PI' )
    PIs.legend(loc='best')
    ks.plot(all_PUs[PU][:-2], all_PUs_sigma[PU], 'o-',  label='sigma')
    #sigmas.legend(loc='best')
    ks.plot(all_PUs[PU][:-1], all_PUs_kappa[PU], 'o-',  label='kappa')
    ks.legend(loc='best')
    f.savefig( CWD + "/results/" + pdb_code + "/PU" + str(PU+ 1) +
                "/PU" +  str(PU+ 1) + "-criterions.pdf")
    plt.close()



#********************* MAIN *********************#

if __name__ == "__main__":
    start = time.time()
    file_path = sys.argv[1]
    #File
    pdb_code = file_path.split('/')[-1].split('.')[0]
    pdb_filename = file_path
    CHAINE = sys.argv[2]
    #Structure of PDB
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_code, pdb_filename)
    model = structure[0]
    # Create directory for results
    dir = CWD +"results/" + pdb_code
    if not os.path.exists(dir):
        os.makedirs(dir)
    #Assignation of Secondary Structure
    dssp_dict = ss_dict(model, pdb_filename)
    #contact probability map
    dist_matrix = calc_dist_matrix(model[CHAINE], model[CHAINE] , dssp_dict)
    contact_map = dist_to_proba(dist_matrix)
    #Dico de PI
    min_size_PU = int(sys.argv[3])
    max_size_PU = int(sys.argv[4])
    dict_PI = all_PI_in_dict(contact_map, dssp_dict, min_size_PU, max_size_PU )
    #save_plot_allPI(dict_PI, pdb_code, max_size_PU)
    dict_PI_max = new_max_PI_in_dict(dict_PI)
    #Define PU
    all_PUs = new_combis_PU(dict_PI_max, min_size_PU, max_size_PU)
    all_PUs = add_start_end_PUs(all_PUs, contact_map)
    # Criterions of PUs
    all_PUs_sigma=[]
    all_PUs_kappa=[]
    all_PUs_PI = []
    all_PUs_PI_index = []
    for PU_i in range(0,len(all_PUs)):
        all_PUs_sigma.append(calc_sigma_decoup(PUs = all_PUs[PU_i]))
        all_PUs_kappa.append(calc_kappa_decoup(PUs = all_PUs[PU_i]))
        all_PUs_PI_index.append(make_liste_PI(dict_PI_max, all_PUs[PU_i])[1])
        all_PUs_PI.append(make_liste_PI(dict_PI_max, all_PUs[PU_i])[0])
    #Write pdb of all splitting combinaison in folder PU1, PU2, PU3...
    #Inside each folder , file -> PU1.1-46.pdb (with start and begin of each PU)
    #write_pdb_all_split(all_PUs)
    for PU in range(len(all_PUs)):
        choix = all_PUs[PU]
        for c in range(1,len(choix)):
            start_res= choix[c-1]
            end_res= choix[c]
            chain_id = CHAINE
            dir = CWD + "/results/" + pdb_code + "/PU" + str(PU + 1)
            file = dir + "/PU" + str(PU + 1) + "." + str(start_res) +"-" + str(end_res) +".pdb"
            if not os.path.exists(dir):
                os.makedirs(dir)
            io = Bio.PDB.PDBIO()
            io.set_structure(structure[0][CHAINE])
            start_res= choix[c-1]
            io.save(file, ResSelect())
        plot_criteria(PU)
        #plot_contact_map(choix)
    #Choose best split
    best = best_split(all_PUs_sigma, all_PUs_kappa)
    bestt = all_PUs[best]
    print("Le meilleur découpage est: PU " + str(best + 1) )
    #Plot and save contact map divided per PUs
    plot_contact_map(bestt)
    #Call Pymol to plot all splitting combinations
    for  PU in range(len(all_PUs)):
        bashCommand = "pymol -qrc src/visu_all_pdbs.py -- {0}/results/{1}/{2}/".format(CWD ,pdb_code, "PU"+str(PU + 1))
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
    end = time.time()
    print("The Protein peeling is done after  " + str(end - start) +  " seconds")
