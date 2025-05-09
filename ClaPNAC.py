import os
from functions import Parsing_PDB, Compare_csv
import argparse

parser = argparse.ArgumentParser(description='Classification of AA-N contacts')
parser.add_argument('-c', '--comparison', help='comparison value: D for difference, U for union')
parser.add_argument('-i', '--in_dir', help='input directory')
parser.add_argument('-t', '--type_in', help='atom type representation: CGR or FA')
parser.add_argument('-m', '--model', help='enter model number: #1')
parser.add_argument('-o', '--output', help='output threshold: float number from 0 to 1')
parser.add_argument('-r', '--range_in', help='range for the ContExt: float number')
parser.add_argument('-s', '--seq', action='store_true', help='sequence output with * when have interaction and . otherwise')

args = parser.parse_args()

print('Start')

compp = args.comparison or "No"
in_dir  = args.in_dir + '/'
type_in = args.type_in or "FA"
output = float(args.output) if args.output else 0.5
range_in = float(args.range_in) if args.range_in else 3.5
seq = 'yes' if args.seq else 'no'
mdl = args.model if args.model else "#1"

count = 0

if compp == "No":
    for path, subdirs, files in os.walk(in_dir):
        for filename in files:
            if filename.endswith(".pdb") or filename.endswith(".cif"):
                count += 1
                print(filename + '\t' + str(count))
                        
                file = os.path.join(path, filename)
                        
                out_file = os.path.join(path, filename.split('.')[0] + '_' + filename.split('.')[1])

                out_file += "_"+str(range_in)+"_"

                if type_in == "CGR":
                    out_file += '_CGR_'
                elif type_in == "FA":
                    out_file += '_FA_'
                else:
                    print("ERROR input")
                    break
                        
                Parsing_PDB.Parse_PDB(file, out_file, type_in, output, range_in, seq, mdl)
else:
    Compare_csv.Compare_csv(in_dir, compp)
        