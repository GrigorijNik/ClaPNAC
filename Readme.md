 # ClaPNAC - program for annotation nucleotide (NA) - amino acid (Protein) contacts.

 Chosen representatives of each class could be found in ./classes_PDB

 As a result you will see the score from 0 to 1 for each class.

```
NA_id	P_id	score	type	N-AA

B-203   A-4     0.63    P       A-GLY

```

 ## Side libraries

 python 3.8

 ClaPNAC requiers **numpy**, **argparse**, **pandas**, **seaborn**, **matplotlib**, **alive_progress**.

 ## Usage

Takes all files in .pdb or .cif from input directory and annotate N-AA contacts.
The results are saving into csv files. ([inputfile]_[representation]_[threshold].csv)

-i input directory

-t full atom (FA) or coarse-grained (CGR) representation will be used

-o threshold for output scores, doublets esimated with the score lower than threshold will be not included in the report file

-r is a range number for ContExt, means the maximum distance for the contacts.

-s option is for sequence: it shows sequence and mark * residues making interactions with score > 0.5

Please see the examples in **test** directory

 ```
 python -h

 python -i [input directory] -t [FA or CGR] -o [float number] -r [float number] -s
 ```

https://doi.org/10.5281/zenodo.15423174
