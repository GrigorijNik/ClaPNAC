import sys
import os
import re
from scipy.spatial import KDTree
from CRTutils import DSSRformat
from classes import FA

N_set  = {'A' : ['P',"OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N9',"C8","N7","C5","C6","N6","N1",'C2','N3',"C4"],
          'ATP' : ['PA',"O1A","O2A","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N9',"C8","N7","C5","C6","N6","N1",'C2','N3',"C4"],
          'ADP' : ['PA',"O1A","O2A","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N9',"C8","N7","C5","C6","N6","N1",'C2','N3',"C4"],
          'AMP' : ['P',"O1P","O2P","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N9',"C8","N7","C5","C6","N6","N1",'C2','N3',"C4"],
          'G' : ['P',"OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N9',"C8","N7","C5","C6","O6","N1",'C2','N2','N3',"C4"],
          'GTP' : ['PA',"O1A","O2A","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N9',"C8","N7","C5","C6","O6","N1",'C2','N2','N3',"C4"],
          'GDP' : ['PA',"O1A","O2A","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N9',"C8","N7","C5","C6","O6","N1",'C2','N2','N3',"C4"],
          'GMP' : ['P',"OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N9',"C8","N7","C5","C6","O6","N1",'C2','N2','N3',"C4"],
          'C' : ['P',"OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N1',"C2","O2","N3","C4","N4","C5",'C6'],
          'CTP' : ['PA',"O1A","O2A","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N1',"C2","O2","N3","C4","N4","C5",'C6'],
          'C5P' : ['PA',"O1A","O2A","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N1',"C2","O2","N3","C4","N4","C5",'C6'],
          'CDP' : ['PA',"O1A","O2A","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N1',"C2","O2","N3","C4","N4","C5",'C6'],
          'CMP' : ['P',"OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N1',"C2","O2","N3","C4","N4","C5",'C6'],
          'U' : ['P',"OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N1',"C2","O2","N3","C4","O4","C5",'C6'],
          'UTP' : ['PA',"O1A","O2A","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N1',"C2","O2","N3","C4","O4","C5",'C6'],
          'UDP' : ['PA',"O1A","O2A","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N1',"C2","O2","N3","C4","O4","C5",'C6'],
          'UMP' : ['P',"OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N1',"C2","O2","N3","C4","O4","C5",'C6'],
          'DA': ['P',"OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N9',"C8","N7","C5","C6","N6","N1",'C2','N3',"C4"],
          'DC': ['P',"OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N1',"C2","O2","N3","C4","N4","C5",'C6'],
          'DT': ['P',"OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N1',"C2","O2","N3","C4","O4","C5",'C6'],
          'DG': ['P',"OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N9',"C8","N7","C5","C6","O6","N1",'C2','N2','N3',"C4"],
          'T' : ['P',"OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N1',"C2","O2","N3","C4","O4","C5",'C6'],
          'TTP' : ['PA',"O1A","O2A","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N1',"C2","O2","N3","C4","O4","C5",'C6'],
          'TYD' : ['PA',"O1A","O2A","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N1',"C2","O2","N3","C4","O4","C5",'C6'],
          'TMP' : ['P',"O1P","O2P","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",'N1',"C2","O2","N3","C4","O4","C5",'C6']}

AA_set = {'ARG' : ['N','CA','C','O','CB','CG','CD','NE','CZ','NH1','NH2'],
          'ALA' : ['N','CA','C','O','CB'],
          'ASN' : ['N','CA','C','O','CB','CG','OD1','ND2'],
          'ASP' : ['N','CA','C','O','CB','CG','OD1','OD2'],
          'CYS' : ['N','CA','C','O','CB','SG'],
          'GLN' : ['N','CA','C','O','CB','CG','CD','OE1','NE2'],
          'GLU' : ['N','CA','C','O','CB','CG','CD','OE1', 'OE2'],
          'GLY' : ['N','CA','C','O'],
          'HIS' : ['N','CA','C','O','CB','CG','ND1','CD2','CE1','NE2'],
          'ILE' : ['N','CA','C','O','CB','CG1','CG2','CD1'],
          'LEU' : ['N','CA','C','O','CB','CG','CD1','CD2'],
          'LYS' : ['N','CA','C','O','CB','CG','CD','CE','NZ'],
          'MET' : ['N','CA','C','O','CB','CG','SD','CE'],
          'PHE' : ['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ'],
          'PRO' : ['N','CA','C','O','CB','CG','CD'],
          'SER' : ['N','CA','C','O','CB','OG'],
          'THR' : ['N','CA','C','O','CB','OG1','CG2'],
          'TRP' : ['N','CA','C','O','CB','CG','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2'],
          'TYR' : ['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ','OH'],
          'VAL' : ['N','CA','C','O','CB','CG1','CG2']}

PDB_file = {}

# Pattern for atom mask parsing
ATOMSPATTERN = re.compile(r"^(#(([0-9]+)(_([0-9]+))?)?)?"+\
                          r"(\/([\w]*))?"+\
                          r"(:([A-Za-z0-9]+)?((_(-?[0-9]+))(_(-?[0-9]+))?)?)?"+\
                          r"(@([A-Za-z0-9]+'?)?((_(-?[0-9]+))(_(-?[0-9]+))?)?)?$")

    
def ParseAtomsFormat(kw, pattern = ATOMSPATTERN):
    """Returns a parsed atom mask object"""
    regsearch = re.findall(ATOMSPATTERN, kw)[0]
    
    modelmin   = regsearch[2]
    modelmax   = regsearch[4]
    chain      = regsearch[6]
    residue    = regsearch[8]
    resnummin  = regsearch[11]
    resnummax  = regsearch[13]
    atom       = regsearch[15]
    atomnummin = regsearch[18]
    atomnummax = regsearch[20]

    res =  {'MODELMIN':   int(modelmin)   if modelmin   else None,
            'MODELMAX':   int(modelmax)   if modelmax   else None,
            'CHAIN':      chain   if chain   else None,
            'RESIDUE':    residue if residue else None,
            'RESNUMMIN':  int(resnummin)  if resnummin  else None,
            'RESNUMMAX':  int(resnummax)  if resnummax  else None,
            'ATOM':       atom    if atom    else None,     
            'ATOMNUMMIN': int(atomnummin) if atomnummin else None,
            'ATOMNUMMAX': int(atomnummax) if atomnummax else None}

    if res['MODELMAX'] is None and not (res['MODELMIN'] is None):
        res['MODELMAX'] = res['MODELMIN']
    if res['RESNUMMAX'] is None and not (res['RESNUMMIN'] is None):
        res['RESNUMMAX'] = res['RESNUMMIN']
    if res['ATOMNUMMAX'] is None and not (res['ATOMNUMMIN'] is None):
        res['ATOMNUMMAX'] = res['ATOMNUMMIN']

    return res


def Allowed(atom, masks):
    """Is the atom allowed by any of the masks?"""
    for mask in masks:

        condition = (mask['ATOMNUMMIN'] is None or atom['id']                 >= mask['ATOMNUMMIN']) and\
                    (mask['ATOMNUMMAX'] is None or atom['id']                 <= mask['ATOMNUMMAX']) and\
                    (mask['ATOM']       is None or atom['auth_atom_id']       == mask['ATOM'])       and\
                    (mask['RESNUMMIN']  is None or atom['auth_seq_id']        >= mask['RESNUMMIN'])  and\
                    (mask['RESNUMMAX']  is None or atom['auth_seq_id']        <= mask['RESNUMMAX'])  and\
                    (mask['RESIDUE']    is None or atom['auth_comp_id']       == mask['RESIDUE'])    and\
                    (mask['CHAIN']      is None or atom['auth_asym_id']       == mask['CHAIN'])      and\
                    (mask['MODELMIN']   is None or atom['pdbx_PDB_model_num'] >= mask['MODELMIN'])   and\
                    (mask['MODELMAX']   is None or atom['pdbx_PDB_model_num'] <= mask['MODELMAX'])
        
        if condition:
            return True
    return False


def ParseAtomPDB(line, model):
    '''https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html'''
    
    atom = {}
    atom["id"]                 = int(line[6:11])
    atom["auth_atom_id"]       = line[12:16].strip()
    atom["label_alt_id"]       = line[16].strip()
    atom["auth_comp_id"]       = line[17:20].strip()
    atom["auth_asym_id"]       = line[20:22].strip()
    atom["auth_seq_id"]        = int(line[22:26])
    atom["pdbx_PDB_ins_code"]  = line[26].strip()
    atom["Cartn_x"]            = float(line[30:38])
    atom["Cartn_y"]            = float(line[38:46])
    atom["Cartn_z"]            = float(line[46:54])
    atom["pdbx_PDB_model_num"] = model

    if atom["pdbx_PDB_ins_code"] == '?':
        atom["pdbx_PDB_ins_code"] = ''

    return atom


def ParsePDB(inpfile, masks):
    """Returns a list of atoms"""
    atoms = []
    model = 1

    with open(inpfile) as file:
        for line in file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom = ParseAtomPDB(line, model)
                if Allowed(atom, masks):
                    atoms.append(atom)
            elif line.startswith('MODEL'):
                model = int(line.strip().split()[-1])
    return atoms

def ParsePDB_N(inpfile, masks):
    """Returns a list of atoms"""
    atoms = []
    seq_N = []
    model = 1

    with open(inpfile) as file:
        for line in file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom = ParseAtomPDB(line, model)
                if Allowed(atom, masks):
                    if atom['auth_comp_id'] in N_set.keys() and atom['auth_atom_id'] in N_set[atom['auth_comp_id']] and (atom['label_alt_id'] == "A" or atom['label_alt_id'] == ""):
                        atoms.append(atom)
                        if atom['auth_asym_id'] not in PDB_file.keys():
                            PDB_file[atom['auth_asym_id']] = {}
                        if atom['auth_comp_id'] not in PDB_file[atom['auth_asym_id']].keys():
                            PDB_file[atom['auth_asym_id']][atom['auth_comp_id']] = {}
                        if atom['auth_seq_id'] not in PDB_file[atom['auth_asym_id']][atom['auth_comp_id']].keys():
                            PDB_file[atom['auth_asym_id']][atom['auth_comp_id']][atom['auth_seq_id']] = []
                            seq_N.append((atom['auth_asym_id'], atom['auth_comp_id'], atom['auth_seq_id']))
                        t = []
                        t.append(atom['Cartn_x'])
                        t.append(atom['Cartn_y'])
                        t.append(atom['Cartn_z'])
                        t.append(atom['auth_atom_id'])
                        PDB_file[atom['auth_asym_id']][atom['auth_comp_id']][atom['auth_seq_id']].append(t)
            elif line.startswith('MODEL'):
                model = int(line.strip().split()[-1])
    return atoms, seq_N

def ParsePDB_P(inpfile, masks):
    """Returns a list of atoms"""
    atoms = []
    seq_P = []
    model = 1

    with open(inpfile) as file:
        for line in file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom = ParseAtomPDB(line, model)
                if Allowed(atom, masks):
                    if atom['auth_comp_id'] in AA_set.keys() and atom['auth_atom_id'] in AA_set[atom['auth_comp_id']] and (atom['label_alt_id'] == "A" or atom['label_alt_id'] == ""):
                        atoms.append(atom)
                        if atom['auth_asym_id'] not in PDB_file.keys():
                            PDB_file[atom['auth_asym_id']] = {}
                        if atom['auth_comp_id'] not in PDB_file[atom['auth_asym_id']].keys():
                            PDB_file[atom['auth_asym_id']][atom['auth_comp_id']] = {}
                        if atom['auth_seq_id'] not in PDB_file[atom['auth_asym_id']][atom['auth_comp_id']].keys():
                            PDB_file[atom['auth_asym_id']][atom['auth_comp_id']][atom['auth_seq_id']] = []
                            seq_P.append((atom['auth_asym_id'], FA.Protein[atom['auth_comp_id']], atom['auth_seq_id']))
                        t = []
                        t.append(atom['Cartn_x'])
                        t.append(atom['Cartn_y'])
                        t.append(atom['Cartn_z'])
                        t.append(atom['auth_atom_id'])
                        PDB_file[atom['auth_asym_id']][atom['auth_comp_id']][atom['auth_seq_id']].append(t)
            elif line.startswith('MODEL'):
                model = int(line.strip().split()[-1])
    return atoms, seq_P


def ParseAtomCIF(line,title):

    linesplit = line.strip().split()
    atom = {title[i]:linesplit[i] for i in range(len(title))}

    for frag in ("atom", "comp", "asym", "seq"):

        auth  =  "auth_{}_id".format(frag)
        label = "label_{}_id".format(frag)
        
        if auth not in atom and label in atom:
            atom[auth] = atom[label]

    for int_token in ("id", "auth_seq_id", "pdbx_PDB_model_num"):
        atom[int_token] = int(atom[int_token]) if int_token in atom else float("nan")

    for float_token in ("Cartn_x", "Cartn_y", "Cartn_z","occupancy","B_iso_or_equiv"):
        atom[float_token] = float(atom[float_token]) if float_token in atom else float("nan")
    
    if      "auth_atom_id" not in atom: atom["auth_atom_id"]      = ''
    if     "label_atom_id" not in atom: atom["label_atom_id"]     = ''
    if      "label_alt_id" not in atom: atom["label_alt_id"]      = ''
    if      "auth_comp_id" not in atom: atom["auth_comp_id"]      = ''
    if     "label_comp_id" not in atom: atom["label_comp_id"]     = ''
    if      "auth_asym_id" not in atom: atom["auth_asym_id"]      = ''
    if "pdbx_PDB_ins_code" not in atom: atom["pdbx_PDB_ins_code"] = ''

    if atom["pdbx_PDB_ins_code"] == '?':
        atom["pdbx_PDB_ins_code"] = ''

    atom["auth_atom_id"] = atom["auth_atom_id"].strip('"')
    atom["label_atom_id"] = atom["label_atom_id"].strip('"')

    return atom


def ParseCIF(inpfile, masks):
    """Returns a list of atoms"""
    atoms = []
    title = []

    with open(inpfile) as file:
        for line in file:

            if line.startswith("_atom_site."):
                title.append(line.strip().split('.')[-1])
            elif line.startswith('ATOM') or line.startswith('HETATM'):
                atom = ParseAtomCIF(line, title)
                if Allowed(atom, masks):
                    atoms.append(atom)
    return atoms

def ParseCIF_N(inpfile, masks):
    """Returns a list of atoms"""
    atoms = []
    seq_N = []
    title = []

    with open(inpfile) as file:
        for line in file:
            if line.startswith("_atom_site."):
                title.append(line.strip().split('.')[-1])
            elif line.startswith('ATOM') or line.startswith('HETATM'):
                atom = ParseAtomCIF(line, title)
                if Allowed(atom, masks):
                    if atom['auth_comp_id'] in N_set.keys() and atom['auth_atom_id'] in N_set[atom['auth_comp_id']]: # and (atom['label_alt_id'] == "A" or atom['label_alt_id'] == ''):
                        atoms.append(atom)
                        if atom['auth_asym_id'] not in PDB_file.keys():
                            PDB_file[atom['auth_asym_id']] = {}
                        if atom['auth_comp_id'] not in PDB_file[atom['auth_asym_id']].keys():
                            PDB_file[atom['auth_asym_id']][atom['auth_comp_id']] = {}
                        if atom['auth_seq_id'] not in PDB_file[atom['auth_asym_id']][atom['auth_comp_id']].keys():
                            PDB_file[atom['auth_asym_id']][atom['auth_comp_id']][atom['auth_seq_id']] = []
                            seq_N.append((atom['auth_asym_id'], atom['auth_comp_id'], atom['auth_seq_id']))
                        t = []
                        t.append(atom['Cartn_x'])
                        t.append(atom['Cartn_y'])
                        t.append(atom['Cartn_z'])
                        t.append(atom['auth_atom_id'])
                        PDB_file[atom['auth_asym_id']][atom['auth_comp_id']][atom['auth_seq_id']].append(t)
    return atoms, seq_N

def ParseCIF_P(inpfile, masks):
    """Returns a list of atoms"""
    atoms = []
    seq_P = []
    title = []

    with open(inpfile) as file:
        for line in file:
            if line.startswith("_atom_site."):
                title.append(line.strip().split('.')[-1])
            elif line.startswith('ATOM') or line.startswith('HETATM'):
                atom = ParseAtomCIF(line, title)
                if Allowed(atom, masks):
                    if atom['auth_comp_id'] in AA_set.keys() and atom['auth_atom_id'] in AA_set[atom['auth_comp_id']]: # and (atom['label_alt_id'] == "A" or atom['label_alt_id'] == ""):
                        atoms.append(atom)
                        if atom['auth_asym_id'] not in PDB_file.keys():
                            PDB_file[atom['auth_asym_id']] = {}
                        if atom['auth_comp_id'] not in PDB_file[atom['auth_asym_id']].keys():
                            PDB_file[atom['auth_asym_id']][atom['auth_comp_id']] = {}
                        if atom['auth_seq_id'] not in PDB_file[atom['auth_asym_id']][atom['auth_comp_id']].keys():
                            PDB_file[atom['auth_asym_id']][atom['auth_comp_id']][atom['auth_seq_id']] = []
                            seq_P.append((atom['auth_asym_id'], FA.Protein[atom['auth_comp_id']], atom['auth_seq_id']))
                        t = []
                        t.append(atom['Cartn_x'])
                        t.append(atom['Cartn_y'])
                        t.append(atom['Cartn_z'])
                        t.append(atom['auth_atom_id'])
                        PDB_file[atom['auth_asym_id']][atom['auth_comp_id']][atom['auth_seq_id']].append(t)
    return atoms, seq_P

def GuessFormat(inpfile):
    """Derives the format of the file; Returns 0 for PDB, 1 for CIF"""
    pdb = 0
    cif = 0

    with open(inpfile) as file:
        for line in file:
            if line.startswith('#'):
                cif += 1
            elif line.startswith('_atom.site'):
                cif += 1
            elif line.startswith('END'):
                pdb += 1
            elif line.startswith('MODEL'):
                pdb += 1

    return int(cif > pdb)

def ParseAtoms(inpfile, masks):
    """Takes masks and the input file, returns a list of the allowed atoms"""
    return (ParsePDB, ParseCIF)[GuessFormat(inpfile)](inpfile, masks)

def ParseAtoms_M(inpfile, type_atoms, masks):
    """Takes masks and the input file, returns a list of the allowed atoms"""
    if GuessFormat(inpfile):
        if type_atoms == 'N':
            return ParseCIF_N(inpfile, masks)
        if type_atoms == 'P':
            return ParseCIF_P(inpfile, masks)
    else:
        if type_atoms == 'N':
            return ParsePDB_N(inpfile, masks)
        if type_atoms == 'P':
            return ParsePDB_P(inpfile, masks)


def Atompairs(n1, n2, limit):
    """takes two sets of 3D coordinates and a threshold, returns all
    the pairs within the threshold"""
    tree1 = KDTree(n1)
    tree2 = KDTree(n2)

    dist = tree1.sparse_distance_matrix(tree2,
                                        limit,
                                        p=2,
                                        output_type='ndarray')
    dist.sort()
    return dist


def FormatAtom(atom):
    """atom -> (DSSR,ATOMID)"""
    return (DSSRformat(atom),atom["auth_atom_id"])


def ReturnContacts(contacts, atoms1, atoms2, skipequal=False):
    """Derives atom info for the contacts"""
    res = []

    for i,j,d in contacts:
        if skipequal and i == j:
            continue
        atom1 = FormatAtom(atoms1[i])
        atom2 = FormatAtom(atoms2[j])
        if atom1[0] != atom2[0]:
            res.append((d, FormatAtom(atoms1[i]), FormatAtom(atoms2[j])))

    return sorted(res)

def ReturnContacts_M(contacts, atoms1, atoms2, skipequal=False):
    """Derives atom info for the contacts"""
    res = []
    exist_res = []

    for i,j,d in contacts:
        ex_res = (str(atoms1[i]['auth_seq_id'])+atoms1[i]['auth_asym_id']+atoms1[i]['auth_comp_id'], str(atoms2[j]['auth_seq_id'])+atoms2[j]['auth_asym_id']+atoms2[j]['auth_comp_id'])
        if ex_res not in exist_res:
            exist_res.append(ex_res)
            temp1 = []
            temp2 = []
            temp2_b = []
            seq_id1 = atoms1[i]['auth_asym_id']
            res1 = atoms1[i]['auth_comp_id']
            res_id1 = atoms1[i]['auth_seq_id']
            seq_id2 = atoms2[j]['auth_asym_id']
            res2 = atoms2[j]['auth_comp_id']
            res_id2 = atoms2[j]['auth_seq_id']
            temp1 = PDB_file[seq_id1][res1][res_id1]
            temp2 = PDB_file[seq_id2][res2][res_id2]

            res.append((atoms1[i]["auth_asym_id"] + '-' + str(atoms1[i]['auth_seq_id']), atoms1[i]['auth_comp_id'], temp1,  atoms2[j]["auth_asym_id"] + '-' + str(atoms2[j]['auth_seq_id']), atoms2[j]['auth_comp_id'], temp2))
    return res


def ContactExtractor(inpfile1, inpfile2 = None, Range = 10.0, mask1 = "#", mask2 = None):
    """Takes input files, a threshold, and atom masks; returns pairs of allowed atoms within a threshold"""
    if inpfile2 is None:
        inpfile2 = inpfile1
    if mask2 is None:
        mask2 = mask1

    mask1 = [ParseAtomsFormat(kw) for kw in mask1]
    mask2 = [ParseAtomsFormat(kw) for kw in mask2]

    if not mask1:
        mask1 = [ParseAtomsFormat(''),]
    if not mask2:
        mask2 = [ParseAtomsFormat(''),]

    PDB_file.clear()
 
    atoms1, sequence_N = ParseAtoms_M(inpfile1, 'N', mask1)
    atoms2, sequence_P = ParseAtoms_M(inpfile1, 'P', mask1)

    onesetflag = False # if two sets are supposed to be identical - just duplicate the first one

    if not atoms1:
        print("No atoms found N")
        return [],[],[]

    if not atoms2:
        print("No atoms2 found P")
        return [],[],[]

    contacts = Atompairs([(a1['Cartn_x'],a1['Cartn_y'],a1['Cartn_z']) for a1 in atoms1],
                         [(a2['Cartn_x'],a2['Cartn_y'],a2['Cartn_z']) for a2 in atoms2],
                         Range)

    return ReturnContacts_M(contacts, atoms1, atoms2, onesetflag), sequence_N, sequence_P