import os

DSSR_PATH = "./x3dna-dssr"

files = [f for f in os.listdir('.') if f.endswith(".pdb")]
files_out = [f.split(".")[0]+"_out.json" for f in os.listdir('.') if f.endswith(".pdb")]

print(files)
print(files_out)

for inp in range(0, len(files)):
    query = "{dssr} -i={input_file} --json --idstr=long -o={output_file} --non-pair=nuc-protein".format(dssr=DSSR_PATH, input_file=files[inp], output_file=files_out[inp])
    os.system(query)
    
print(query)
