#  Scoring function as combination of SF_RMSD and SF_PDis
import numpy as np

def Scoring_F(input_matrix, class_matrix, class_rmsd, class_median):

    min_X    = 0.0
    median_X = 0.0

    tmp = []
    tmp_class_rmsd = class_rmsd.copy()
    temp_class_rmsd = class_rmsd.copy()
    temp_class_median = class_median
    temp_class_matrix = class_matrix.copy()

    temp_class_rmsd = np.array(temp_class_rmsd)
    temp_class_matrix = np.array(temp_class_matrix)
    input_matrix = np.array(input_matrix)

    for sample in temp_class_matrix:
        temp = np.mean(np.sqrt(np.sum((input_matrix-sample) ** 2, axis = 1)))
        tmp_class_rmsd.append(temp)
        tmp.append(temp)

    tmp               = np.array(tmp)
    tmp_class_rmsd    = np.array(tmp_class_rmsd)

    tmp_class_rmsd    = tmp_class_rmsd[tmp_class_rmsd != 0.0]

    min_X             = np.min(tmp)
    # max_X             = np.max(tmp)
    # max_class         = np.max(temp_class_rmsd)
    median_X          = np.median(tmp_class_rmsd)

    res = temp_class_median / (min_X + median_X) #+ 0.25*max_class / (min_X + max_X)

    return res