import os
from alive_progress import alive_bar

def Compare_csv(in_dir, compp):
    comp_files = []
    comp_f = []
    for path, subdirs, files in os.walk(in_dir):
        for filename in files:
            if filename.endswith(".csv"):
                comp_files.append(filename)
                f = open(os.path.join(path, filename), "r")
                f_lines = [line.rstrip() for line in f]
                if len(comp_f) < 2:
                    comp_f = f_lines
                else:
                    if compp == "U":
                        comp_f = [item for item in f_lines if item in comp_f]
                    elif compp == "D":
                        comp_f = list(set(comp_f).symmetric_difference(set(f_lines)))
                    else:
                        exit()

    res_file = open(in_dir+"comparison_"+str(compp)+".txt", "w")
    if compp == "D":
        res_file.write('NA_id,P_id,score,type,N-AA,\n')
    for line in comp_f:
        res_file.write(line+"\n")

    print(comp_files)