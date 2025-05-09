# Parsing file with uslib2
from functions import Suprimp, SF, Contact_Map
from classes import CGR, FA
import ContExt_CGR
import ContExt_FA
from alive_progress import alive_bar
import numpy as np

def Parse_PDB(file, out_file, type_in, output, range_in, seq, mdl):
    res_ContExt = []
    seq_N_dot = []
    seq_P_dot = []
    contact_map = np.zeros((0, 0), dtype=int)

    if type_in == "CGR":
        res_ContExt, seq_N, seq_P = ContExt_CGR.ContactExtractor(inpfile1=file, Range=range_in, mask1 = [mdl,])  # Extraction of all doublets
    elif type_in == "FA":
        res_ContExt, seq_N, seq_P = ContExt_FA.ContactExtractor(inpfile1=file, Range=range_in, mask1 = [mdl,])  # Extraction of all doublets

    if seq_N and seq_P:
        seq_N_dot = ['.']*len(seq_N)
        seq_P_dot = ['.']*len(seq_P)
        contact_map = np.zeros((len(seq_N_dot), len(seq_P_dot)), dtype=int)

    res_file = open(out_file+str(output)+".csv", "w") # File with results of classifying
    # test_file = open("./clarnp_crossvalid/test_FA_base_SF.csv", 'a') # for testing only

    res_file.write('NA_id,P_id,score,type,N-AA,\n')

    if len(res_ContExt) == 0:
        res_file.write('No close contacts have been detected!')
        print('No close contacts have been detected!')

    def compute():
        num_pairs_classified = 0
        with alive_bar(len(res_ContExt)) as bar:
            for x in res_ContExt:
                type_score = type_in
                AA_arr = []
                AA_arr.append(x[4])
                N_type  = x[1]
                N_type1 = N_type
                if N_type in ["T", "DT", "TTP", "TYD", "TMP"]:
                    N_type = "U"
                    N_type1 = "T"
                elif N_type in ["UTP", "UDP", "UMP"]:
                    N_type = "U"
                    N_type1 = "U"
                elif N_type in ["DA", "ATP", "ADP", "AMP"]:
                    N_type = "A"
                    N_type1 = "A"
                elif N_type in ["DC", "CTP", "CDP", "CMP", "C5P"]:
                    N_type = "C"
                    N_type1 = "C"
                elif N_type in ["DG", "GTP", "GDP", "GMP"]:
                    N_type = "G"
                    N_type1 = "G"
                result_SF = {'P' : -1.0, 'R' : -1.0, "3'" : -1.0, "5'" : -1.0, 'H' : -1.0, 'S' : -1.0, 'W' : -1.0}
                map_class = {'P' : 1, 'R' : 2, "3'" : 6, "5'" : 6, 'H' : 5, 'S' : 4, 'W' : 3, 'M' : 7}
                
                AA_arr.append('Backbone')
                
                if len(x[2]) != 5 and type_score == "CGR" or len(x[5])-3 != FA.FA[x[4]] and type_score == "FA" and x[4] != 'GLY':
                    print('Problems with input! This pair skipped!')
                    print(x[0], x[1], len(x[2]), x[3], x[4], len(x[5]))
                else:
                    for AA_type in AA_arr:
                        for cl_type in CGR.classes_type:
                            class_AA_matrix_rmsd, class_AA_matrix_pd, class_AA_median_rmsd, class_AA_median_pd, class_AA_rmsd_rmsd, class_AA_rmsd_pd = Suprimp.Switch_AA_median_rmsd(AA_type, N_type, cl_type, type_score)
                            if class_AA_matrix_rmsd == 'ERROR' or class_AA_median_rmsd == 'ERROR' or class_AA_rmsd_rmsd == 'ERROR':
                                print(N_type, "-", AA_type,'No such amino acids for now in ClaRNP...')
                            else:
                                try:
                                    AA_sample_rmsd = Suprimp.Suprimp_F(x, AA_type, N_type, cl_type, type_score)
                                    AA_sample_pd   = Suprimp.Pdis(x, AA_type, N_type, cl_type, type_score)
                                    # print(AA_sample_rmsd, AA_sample_pd)
                                    if len(AA_sample_rmsd) == 0 or len(AA_sample_pd) == 0:
                                        # print(AA_sample_rmsd, AA_sample_pd)
                                        result_SF[cl_type] = -1.0
                                    else:
                                        t_rmsd = SF.Scoring_F(AA_sample_rmsd, class_AA_matrix_rmsd, class_AA_rmsd_rmsd, class_AA_median_rmsd)
                                        t_pd   = SF.Scoring_F(AA_sample_pd, class_AA_matrix_pd, class_AA_rmsd_pd, class_AA_median_pd)
                                        t_sf  = (t_rmsd + t_pd) / 2 # Scoring function
                                        result_SF[cl_type] = t_sf
                                except:
                                    print("Problem with this pair:", x[0], x[1], len(x[2]), x[3], x[4])
# ---------------------------------------------Only for testing------------------------------------------------------------------------------
                        # temp_cl = max(result_SF, key=result_SF.get)
                        # temp_sf = result_SF[temp_cl]
                        # temp_sf = float(f'{temp_sf:.2f}')
                        # if temp_sf > output: # and AA_type == 'Backbone': # for testing backbone interactions
                        #     test_file.write('\n')
                        #     test_file.write(file.split('/')[-2] + ',' + file.split('/')[-1] + ',' + AA_type + ',' + temp_cl + ',' + str(temp_sf))
# -------------------------------------------------------------------------------------------------------------------------------------------
                        result_SF_list = sorted(result_SF.items(), key=lambda x: x[1], reverse=True)
                        # print(result_SF_list)
                        # count_res = 0

                        for res_t in result_SF_list:
                            t = float(f'{res_t[1]:.2f}')
                            if res_t[0] in "PR":
                                k = t / FA.Tholds[AA_type][res_t[0]]
                                k = float(f'{k:.2f}')
                                if k > 1.0:
                                    result_SF[res_t[0]] = 1
                                else:
                                    result_SF[res_t[0]] = k
                            else:
                                k = t / FA.Tholds[AA_type][N_type+res_t[0][0]]
                                k = float(f'{k:.2f}')
                                if k > 1.0:
                                    result_SF[res_t[0]] = 1
                                else:
                                    result_SF[res_t[0]] = k

                        result_SF_list = sorted(result_SF.items(), key=lambda x: x[1], reverse=True)

                        for res_t in result_SF_list:
                            t = float(f'{res_t[1]:.2f}')
                            if t >= output:
                                res_file.write(x[0] + ',')
                                res_file.write(x[3] + ',')
                                Int_type = res_t[0]
                                if Int_type == "5'":
                                    Int_type = "+"
                                if Int_type == "3'":
                                    Int_type = "-"
                                res_file.write(str(t) + ',' + Int_type + ',')
                                if AA_type == 'Backbone':
                                    res_file.write(N_type1 + '-' + AA_arr[0] + '_' + AA_type + ',')
                                else:
                                    res_file.write(N_type1 + '-' + AA_type + ',')
                                if seq == 'yes':
                                    seq_N_dot[seq_N.index((x[0].split('-')[0], x[1], int(x[0].split('-')[1])))] = '*'
                                    seq_P_dot[seq_P.index((x[3].split('-')[0], FA.Protein[x[4]], int(x[3].split('-')[1])))] = '*'
                                    if contact_map[seq_N.index((x[0].split('-')[0], x[1], int(x[0].split('-')[1])))][seq_P.index((x[3].split('-')[0], FA.Protein[x[4]], int(x[3].split('-')[1])))] == 0:
                                        contact_map[seq_N.index((x[0].split('-')[0], x[1], int(x[0].split('-')[1])))][seq_P.index((x[3].split('-')[0], FA.Protein[x[4]], int(x[3].split('-')[1])))] = map_class[res_t[0]]
                                    else:
                                        contact_map[seq_N.index((x[0].split('-')[0], x[1], int(x[0].split('-')[1])))][seq_P.index((x[3].split('-')[0], FA.Protein[x[4]], int(x[3].split('-')[1])))] = map_class['M']
                                res_file.write('\n')

                        num_pairs_classified += 1
                bar()

        print("Pairs classified/extracted with ContExt :", int(num_pairs_classified/2), "/", len(res_ContExt))

        if seq == 'yes':
            res_file.write('\n' + ",,,,,")
            seq_N_dot1 = ''.join(seq_N_dot)
            pair_prev = seq_N[0]
            i_prev = 0
            i = 0
            for pair in seq_N:
                if pair[0] == pair_prev[0]:
                    print(pair[1], end="")
                    res_file.write(pair[1])
                    i+=1
                else:
                    print()
                    res_file.write('\n' + ",,,,,")
                    print(seq_N_dot1[i_prev:i])
                    res_file.write(seq_N_dot1[i_prev:i])
                    res_file.write('\n' + ",,,,,")
                    i_prev = i
                pair_prev = pair
            print()
            res_file.write('\n' + ",,,,,")
            print(seq_N_dot1[i_prev:i])
            res_file.write(seq_N_dot1[i_prev:i])

            res_file.write('\n' + ",,,,,")

            seq_P_dot1 = ''.join(seq_P_dot)
            pair_prev = seq_N[0]
            i_prev = 0
            i = 0
            for pair in seq_P:
                if pair[0] == pair_prev[0]:
                    print(pair[1], end="")
                    res_file.write(pair[1])
                    i+=1
                else:
                    print()
                    res_file.write('\n' + ",,,,,")
                    print(seq_P_dot1[i_prev:i])
                    res_file.write(seq_P_dot1[i_prev:i])
                    res_file.write('\n' + ",,,,,")
                    i_prev = i
                pair_prev = pair
            print()
            res_file.write('\n' + ",,,,,")
            print(seq_P_dot1[i_prev:i])
            res_file.write(seq_P_dot1[i_prev:i])

            print()

            Contact_Map.Save_Map(out_file, contact_map, seq_P, seq_N, output)

    compute()

    res_file.close()
    # test_file.close() # for testing only