#  Functions for Superimposing

import numpy as np
from classes import CGR, ARG_class, ALA_class, ASN_class, ASP_class, Backbone_class, CYS_class, GLN_class, GLU_class, GLY_class, HIS_class, ILE_class,\
                    LEU_class, LYS_class, MET_class, PHE_class, PRO_class, SER_class, THR_class, TRP_class, TYR_class, VAL_class,\
                    FA, ARG_class_FA, ALA_class_FA, ASN_class_FA, ASP_class_FA, CYS_class_FA, GLN_class_FA, GLU_class_FA, GLY_class_FA, HIS_class_FA, ILE_class_FA,\
                    LEU_class_FA, LYS_class_FA, MET_class_FA, PHE_class_FA, PRO_class_FA, SER_class_FA, THR_class_FA, TRP_class_FA, TYR_class_FA, VAL_class_FA, Backbone_class_FA

def get_transform(X:'np.ndarray', Y:'np.ndarray'):
    X_avg = X.mean(axis=0)
    Y_avg = Y.mean(axis=0)
    
    X = X - X_avg
    Y = Y - Y_avg
    
    M = np.dot(np.transpose(Y), X)
    S, V, D = np.linalg.svd(M)
    
    A = np.transpose(np.dot(np.transpose(D), np.transpose(S)))
    if np.linalg.det(A) < 0:
        D[2] = -D[2]
        A = np.transpose(np.dot(np.transpose(D), np.transpose(S)))
    B = X_avg - np.dot(Y_avg, A)
    
    return A, B

def apply_transform(X:'np.ndarray', transform):
    A, B = transform
    return np.dot(X, A) + B

def Switch_suprimp(AA_type, N_type, class_AA_type, repr_atom): # Switch atoms for Superimposing
    if repr_atom == "CGR":
        if class_AA_type in ("P", "R"):
            str1 = AA_type + "_class.class_" + class_AA_type
        else:
            str1 = AA_type + "_class.class_" + N_type
        return eval(str1)
    elif repr_atom == "FA":
        if class_AA_type in ("P", "R"):
            str1 = AA_type + "_class_FA.class_" + class_AA_type
        else:
            str1 = AA_type + "_class_FA.class_" + N_type
        return eval(str1)

def Switch_type(AA_type, repr_atom, class_AA_type): # Switch coarse-grained representation of amino acid
    if repr_atom == "CGR":
        str1 = "CGR." + AA_type
        return eval(str1)
    elif repr_atom == "FA" and class_AA_type not in ("3'", "5'"):
        str1 = "FA." + AA_type
        return eval(str1)
    else:
        str1 = "CGR." + AA_type + "_stack"
        return eval(str1)
    
def Switch_AA_median_rmsd(AA_type, N_type, class_AA_type, repr_atom): # Switch class of interactions
    matrix_rmsd = "ERROR"
    matrix_pd   = "ERROR"
    median_rmsd = ""
    median_pd   = ""
    rmsd_rmsd   = ""
    rmsd_pd     = ""
    if class_AA_type == "3'" :
            class_AA_type = "3"
    elif class_AA_type == "5'" :
            class_AA_type = "5"
    if repr_atom == "CGR":
        if class_AA_type in ('P', 'R'):
            matrix_rmsd = AA_type + "_class.sample_" + class_AA_type
            median_rmsd = AA_type + "_class.median_" + class_AA_type
            rmsd_rmsd   = AA_type + "_class.rmsd_" + class_AA_type
            matrix_pd   = AA_type + "_class.sample_" + class_AA_type + "_pd"
            median_pd   = AA_type + "_class.median_" + class_AA_type + "_pd"
            rmsd_pd     = AA_type + "_class.rmsd_" + class_AA_type + "_pd"
        else:
            matrix_rmsd = AA_type + "_class.sample_" + N_type + class_AA_type
            median_rmsd = AA_type + "_class.median_" + N_type + class_AA_type
            rmsd_rmsd   = AA_type + "_class.rmsd_" + N_type + class_AA_type
            matrix_pd   = AA_type + "_class.sample_" + N_type + class_AA_type + "_pd"
            median_pd   = AA_type + "_class.median_" + N_type + class_AA_type + "_pd"
            rmsd_pd     = AA_type + "_class.rmsd_" + N_type + class_AA_type + "_pd"
    elif repr_atom == "FA":
        if class_AA_type in ('P', 'R'):
            matrix_rmsd = AA_type + "_class_FA.sample_" + class_AA_type
            median_rmsd = AA_type + "_class_FA.median_" + class_AA_type
            rmsd_rmsd   = AA_type + "_class_FA.rmsd_" + class_AA_type
            matrix_pd = AA_type + "_class_FA.sample_" + class_AA_type + "_pd"
            median_pd = AA_type + "_class_FA.median_" + class_AA_type + "_pd"
            rmsd_pd   = AA_type + "_class_FA.rmsd_" + class_AA_type + "_pd"
        else:
            matrix_rmsd = AA_type + "_class_FA.sample_" + N_type + class_AA_type
            median_rmsd = AA_type + "_class_FA.median_" + N_type + class_AA_type
            rmsd_rmsd   = AA_type + "_class_FA.rmsd_" + N_type + class_AA_type
            matrix_pd = AA_type + "_class_FA.sample_" + N_type + class_AA_type + "_pd"
            median_pd = AA_type + "_class_FA.median_" + N_type + class_AA_type + "_pd"
            rmsd_pd   = AA_type + "_class_FA.rmsd_" + N_type + class_AA_type + "_pd"

    if matrix_rmsd == "ERROR" or matrix_pd == "ERROR":
        return "ERROR", "ERROR", "ERROR", "ERROR", "ERROR", "ERROR"

    return eval(matrix_rmsd), eval(matrix_pd), eval(median_rmsd), eval(median_pd), eval(rmsd_rmsd), eval(rmsd_pd)

def Suprimp_F(NA_pair, AA_type, N_type, class_AA_type, repr_atom): # Return atoms of AA after superimposing

    AA = Switch_type(AA_type, repr_atom, class_AA_type) # Understand type of amino acid residue
    atoms_AA_sample = []
    temp_repr = repr_atom
    if temp_repr == "CGR":
        if len(NA_pair[2]) != 5 or len(NA_pair[5]) != len(AA): # Check number of atoms
            return atoms_AA_sample
        else:
            if class_AA_type in ('P', 'R'):        
                atoms_nucl_sample    = NA_pair[2][0:3].copy()
            else:
                atoms_nucl_sample    = NA_pair[2][2:5].copy()

            atoms_nucl_reference = Switch_suprimp(AA_type, N_type, class_AA_type, repr_atom) # Pick atoms for superimposing
            sample_atoms = NA_pair[2].copy()
            for i in range(0, len(AA)):
                sample_atoms.append(NA_pair[5][i])
            sample_atoms.append(NA_pair[6][0])
            sample_atoms.append(NA_pair[6][1])
            
            sample_atoms = np.array(sample_atoms)
            atoms_nucl_sample = np.array(atoms_nucl_sample)
            atoms_nucl_reference = np.array(atoms_nucl_reference)

            atoms_AA_sample      = apply_transform(sample_atoms, get_transform(atoms_nucl_reference, atoms_nucl_sample)) # Superimpose

            l = 5 + len(AA)
            atoms_AA_sample = np.array(atoms_AA_sample)

            if class_AA_type in ("3'", "5'") and AA_type in ("HIS", "PHE", "TYR", "TRP", "PRO"): # Additional check for stacking: base and sidechain should be parallel
                # Check planarity of AA atoms and base atoms
                atoms_AA_temp = atoms_AA_sample[(l-2):l]
                v_lin = atoms_AA_temp[0] - atoms_AA_temp[1]
                v_lin = np.array(v_lin)
                normalized_v_lin = v_lin / np.sqrt(np.sum(v_lin ** 2))
                v_plane1 = atoms_nucl_sample[1] - atoms_nucl_sample[0]
                v_plane1 = np.array(v_plane1)
                v_plane2 = atoms_nucl_sample[2] - atoms_nucl_sample[0]
                v_plane2 = np.array(v_plane2)

                v_norm = np.cross(v_plane1, v_plane2)
                v_norm = np.array(v_norm)
                normalized_v_norm = v_norm / np.sqrt(np.sum(v_norm ** 2))

                planarity = np.dot(normalized_v_lin, normalized_v_norm)

                if planarity > -0.1 and planarity < 0.1:
                    atoms_AA_sample = []
                    return atoms_AA_sample

            return atoms_AA_sample[5:l]
    elif temp_repr == "FA":
        if len(NA_pair[5]) != len(AA)+3 and AA_type != 'Backbone' and class_AA_type not in ("3'", "5'") and AA_type != 'GLY': # Check number of atoms
            return atoms_AA_sample
        atoms_nucl_sample = []

        if class_AA_type in ('P', 'R') and NA_pair[2][0][3] == 'P':
            for itemp in NA_pair[2]:
                if itemp[3] in ("P", "C4'", "C1'"):
                    atoms_nucl_sample.append(itemp[:-1])
        elif class_AA_type in ('P', 'R') and NA_pair[2][0][3] != 'P':
            for itemp in NA_pair[2]:
                if itemp[3] in ("O5'", "C4'", "C1'"):
                    atoms_nucl_sample.append(itemp[:-1])
        elif N_type in "AG":
            for itemp in NA_pair[2]:
                if itemp[3] in ("N9", "C6", "C2'"):
                    atoms_nucl_sample.append(itemp[:-1])
        else:
            for itemp in NA_pair[2]:
                if itemp[3] in ("N1", "C2", "C4"):
                    atoms_nucl_sample.append(itemp[:-1])

        if len(atoms_nucl_sample) < 3:
            return atoms_AA_sample

        atoms_nucl_reference = Switch_suprimp(AA_type, N_type, class_AA_type, repr_atom) # Pick atoms for superimposing
        sample_atoms = [ttt[:-1].copy() for ttt in NA_pair[2]]
        AA_coord = []
        prot_atom_type = []
        for row in NA_pair[5]:
            prot_atom_type.append(row[3])
            AA_coord.append(row[0:3])
        for i in range(0, len(AA_coord)):
            sample_atoms.append(AA_coord[i])
            
        sample_atoms = np.array(sample_atoms)
        atoms_nucl_sample = np.array(atoms_nucl_sample)
        atoms_nucl_reference = np.array(atoms_nucl_reference)

        try:
            atoms_AA_sample = apply_transform(sample_atoms, get_transform(atoms_nucl_reference, atoms_nucl_sample))
        except:
            AA_return = []
            return AA_return

        AA_return = []
        AA_temp = []
        nl = len(NA_pair[2])
        if AA_type == 'Backbone':
            AA_return = [atoms_AA_sample[nl+i] for i in [0, 3]]
            AA_temp = [AA_coord[i] for i in [0, 3]]
        elif AA_type == 'GLY' and class_AA_type not in ("3'", "5'"):
            AA_return = [atoms_AA_sample[nl+i] for i in [0, 1, 2, 3]]
            AA_temp = [AA_coord[i] for i in [0, 1, 2, 3]]
        elif AA_type == 'GLY' and class_AA_type in ("3'", "5'"):
            AA_return = [atoms_AA_sample[nl+i] for i in [1, 2]]
            AA_temp = [AA_coord[i] for i in [1, 2]]
        else:
            AA_return.append(atoms_AA_sample[1+nl])
            AA_temp.append(AA_coord[1])
            for i in range(1, len(AA)):
                AA_return.append(atoms_AA_sample[i+3+nl])
                AA_temp.append(atoms_AA_sample[i+3])

        AA_return = np.array(AA_return)

        if class_AA_type in ("3'", "5'") and AA_type in ("HIS", "PHE", "TYR", "TRP", "PRO"): # Additional check for stacking: base and sidechain should be parallel
                # Check planarity of AA atoms and base atoms
            atoms_AA_temp = []
            atoms_AA_temp.append(AA_return[0])
            atoms_AA_temp.append(AA_return[1])
            v_lin = atoms_AA_temp[0] - atoms_AA_temp[1]
            v_lin = np.array(v_lin)
            normalized_v_lin = v_lin / np.sqrt(np.sum(v_lin ** 2))
            v_plane1 = atoms_nucl_sample[2] - atoms_nucl_sample[0]
            v_plane1 = np.array(v_plane1)
            v_plane2 = atoms_nucl_sample[1] - atoms_nucl_sample[0]
            v_plane2 = np.array(v_plane2)

            v_norm = np.cross(v_plane1, v_plane2)
            v_norm = np.array(v_norm)
            normalized_v_norm = v_norm / np.sqrt(np.sum(v_norm ** 2))

            planarity = np.dot(normalized_v_lin, normalized_v_norm)

            if planarity > -0.1 and planarity < 0.1:
                AA_return = []
                return AA_return
            
        return AA_return

def Pdis(NA_pair, AA_type, N_type, class_AA_type, repr_atom): # Return atoms of AA with pairwaith distances
    AA = Switch_type(AA_type, repr_atom, class_AA_type)
    pd_matrix = []
    atoms_AA_sample = []

    if repr_atom == "CGR":
        atoms_nucl_sample = NA_pair[2].copy()
        if len(NA_pair[2]) != 5 or len(NA_pair[5]) != len(AA): # Check number of atoms
            return pd_matrix
        elif class_AA_type in ('P', 'R', 'W', 'H', 'S', "5'", "3'") and AA_type != 'Backbone':
            atoms_AA_sample   = NA_pair[5].copy()
        else:
            atoms_AA_sample.append(NA_pair[6][0])
            atoms_AA_sample.append(NA_pair[6][1])

        atoms_nucl_sample = np.array(atoms_nucl_sample)

        atoms_AA_sample = np.array(atoms_AA_sample)

        for i in range (0, atoms_AA_sample.shape[0]):
            tmp = []
            for j in range (0, atoms_nucl_sample.shape[0]):
                tmp.append(np.sqrt(np.sum((atoms_AA_sample[i] - atoms_nucl_sample[j]) ** 2))) # Calculate pd
            pd_matrix.append(tmp)
        
        return pd_matrix
    elif repr_atom == "FA":
        # if len(NA_pair[5]) != len(AA)+3 and AA_type != 'Backbone' and class_AA_type not in ("3'", "5'") and AA_type != 'GLY': # Check number of atoms
        #     return atoms_AA_sample
        # if class_AA_type in ('P', 'R'):
        #     atoms_nucl_sample = [ttt[:-1].copy() for ttt in NA_pair[2][:13]]
        # else:
        #     atoms_nucl_sample = [ttt[:-1].copy() for ttt in NA_pair[2][13:]]
        atoms_nucl_sample = []

        # print(N_type, class_AA_type)

        # for itemp in NA_pair[2]:
        #     print(itemp[3])

        if class_AA_type in ('P', 'R') and NA_pair[2][0][3] == 'P' and NA_pair[2][10][3] != "C1'":
            for itemp in NA_pair[2]:
                if N_type in ("A", "G"):
                    if itemp[3] in ('P',"OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N9'):
                        atoms_nucl_sample.append(itemp[:-1])
                elif N_type in ("C", "U"):
                    if itemp[3] in ('P',"OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N1'):
                        atoms_nucl_sample.append(itemp[:-1])
        elif class_AA_type in ('P', 'R') and NA_pair[2][0][3] == "O5'" and NA_pair[2][7][3] != "C1'":
            atoms_nucl_sample.append(NA_pair[2][0][:-1])
            atoms_nucl_sample.append(NA_pair[2][0][:-1])
            atoms_nucl_sample.append(NA_pair[2][0][:-1])
            for itemp in NA_pair[2]:
                if N_type in ("A", "G"):
                    if itemp[3] in ("O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N9'):
                        atoms_nucl_sample.append(itemp[:-1])
                elif N_type in ("C", "U"):
                    if itemp[3] in ("O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N1'):
                        atoms_nucl_sample.append(itemp[:-1])
        elif class_AA_type in ('P', 'R') and NA_pair[2][0][3] == "O5'" and NA_pair[2][7][3] == "C1'":
            atoms_nucl_sample.append(NA_pair[2][0][:-1])
            atoms_nucl_sample.append(NA_pair[2][0][:-1])
            atoms_nucl_sample.append(NA_pair[2][0][:-1])
            for itemp in NA_pair[2]:
                if N_type in ("A", "G"):
                    if itemp[3] in ("O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N9'):
                        atoms_nucl_sample.append(itemp[:-1])
                    if itemp[3] == "C2'":
                        atoms_nucl_sample.append(itemp[:-1])
                elif N_type in ("C", "U"):
                    if itemp[3] in ("O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N1'):
                        atoms_nucl_sample.append(itemp[:-1])
                    if itemp[3] == "C2'":
                        atoms_nucl_sample.append(itemp[:-1])
        elif class_AA_type in ('P', 'R') and NA_pair[2][0][3] == 'P' and NA_pair[2][10][3] == "C1'":
            for itemp in NA_pair[2]:
                if N_type in ("A", "G"):
                    if itemp[3] in ('P',"OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N9'):
                        atoms_nucl_sample.append(itemp[:-1])
                    if itemp[3] == "C2'":
                        atoms_nucl_sample.append(itemp[:-1])
                elif N_type in ("C", "U"):
                    if itemp[3] in ('P',"OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N1'):
                        atoms_nucl_sample.append(itemp[:-1])
                    if itemp[3] == "C2'":
                        atoms_nucl_sample.append(itemp[:-1])
        elif N_type == "A":
            for itemp in NA_pair[2]:
                if itemp[3] in ("C8","N7","C5","C6","N6","N1",'C2','N3',"C4"):
                    atoms_nucl_sample.append(itemp[:-1])
        elif N_type == "C":
            for itemp in NA_pair[2]:
                if itemp[3] in ("C2","O2","N3","C4","N4","C5",'C6'):
                    atoms_nucl_sample.append(itemp[:-1])
        elif N_type == "G":
            for itemp in NA_pair[2]:
                if itemp[3] in ("C8","N7","C5","C6","O6","N1",'C2','N2','N3',"C4"):
                    atoms_nucl_sample.append(itemp[:-1])
        elif N_type == "U":
            for itemp in NA_pair[2]:
                if itemp[3] in ("C2","O2","N3","C4","O4","C5",'C6'):
                    atoms_nucl_sample.append(itemp[:-1])

        # print(len(atoms_nucl_sample))

        if len(atoms_nucl_sample) != 13 and class_AA_type in ('P', 'R'):
            return atoms_AA_sample
        
        if len(atoms_nucl_sample) != len(eval("FA."+N_type))-13 and class_AA_type not in ('P', 'R'):
            return atoms_AA_sample

        AA_coord = []
        prot_atom_type = []
        for row in NA_pair[5]:
            prot_atom_type.append(row[3])
            AA_coord.append(row[0:3])
        if AA_type != 'Backbone' and AA_type != 'GLY':
            for i in range(0, len(prot_atom_type)): # Pick atoms for superimposing if we check stacking
                if prot_atom_type[i] in AA:
                    atoms_AA_sample.append(AA_coord[i])
        elif AA_type == 'Backbone':
            atoms_AA_sample.append(AA_coord[0])
            atoms_AA_sample.append(AA_coord[3])
        elif AA_type == 'GLY':
            for i in range(0, len(AA)):
                atoms_AA_sample.append(AA_coord[i])

        atoms_nucl_sample = np.array(atoms_nucl_sample)

        atoms_AA_sample = np.array(atoms_AA_sample)

        for i in range (0, atoms_AA_sample.shape[0]):
            tmp = []
            for j in range (0, atoms_nucl_sample.shape[0]):
                tmp.append(np.sqrt(np.sum((atoms_AA_sample[i] - atoms_nucl_sample[j]) ** 2))) # Calculate pd
            pd_matrix.append(tmp)
        
        return pd_matrix
            