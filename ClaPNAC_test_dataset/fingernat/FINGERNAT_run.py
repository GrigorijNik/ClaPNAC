import os

fingernat_PATH = "code/fingeRNAt.py"
out_dir = './my_test/'
files = [f for f in os.listdir(out_dir) if f.endswith(".pdb")]


for inp in files:
    query = "python3 {fingernat} -r {rec_file} -l {lig_file} -o {out_dir} -detail".format(fingernat=fingernat_PATH, rec_file=out_dir+inp, lig_file=out_dir+inp.split("_")[0]+"_prot_split.sdf", out_dir=out_dir)
    os.system(query)
    
print(query)
