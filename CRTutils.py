import sys

def DSSRformat(atom):
    '''1.A.G.22.B - DSSR-like format
       [model.chain.res.resnum.inscode]'''
    return '.'.join([str(x) for x in [atom["pdbx_PDB_model_num"],
                                      atom["auth_asym_id"],
                                      atom["auth_comp_id"],
                                      atom["auth_seq_id"],
                                      atom["pdbx_PDB_ins_code"]]])


def Strands(model, pairs):
    """Derive continuous RNA strands based on O3'-P atom pairs under 2 angstroms"""
    strands = []

    # Ignore residues that contain no O3' and no P atoms
    residues = {res['DSSR'] for res in model["RES"] if "P" in res['ATOMDICT'] or "O3'" in res['ATOMDICT']}

    # List of directed edges
    edges = [(p[1][0],p[2][0]) for p in pairs if p[0] < 2.0]

    # Adjective list
    adj = {x:y for x,y in edges}

    # Sources (strand starts) = all the residues with no incoming edges
    srcs = residues - set(adj.values())

    # If no sources found, but edges list is not empty - ignore the last edge
    # (this is for the case of a single circular RNA)
    # TODO: handle multiple circular RNAs
    if not srcs and edges:
        adj = {x:y for x,y in edges[:-1]}
        srcs = residues - set(adj.values())

    # Derive strands as continuous paths from a source
    while srcs:

        v = srcs.pop()
        curstrand = [v,]
        while v in adj:
            curstrand.append(adj[v])
            v = adj.pop(v)
        strands.append(tuple(curstrand))

    # Return strands sorted by (chain,resnum) of the residues
    return tuple(sorted(strands,key = lambda x: (x[0].split('.')[:2],
                                           int(x[0].split('.')[3]),
                                           x[0].split('.')[4])))


def CandPerms(pairs, strands, distlim = 8.0, minseqdist = 3, interonly = False, intraonly = False):
    """Returns a list of breaks"""
    dssrloc = {}

    perms = []

    # dssrloc = {dssr: (strand index, dssr index within the strand)}
    for i in range(len(strands)):
        for j in range(len(strands[i])):
            dssrloc[strands[i][j]] = (i, j)

    # Check all the O3'-P pairs
    for pair in sorted(pairs):

        dist, dssr1, dssr2 = pair[0], pair[1][0], pair[2][0]

        i1, j1 = dssrloc[dssr1]
        i2, j2 = dssrloc[dssr2]

        # If the distance is larger than 2 angstroms (i.e. the pair is non-neighboring)
        # and under the distlim threshold
        if 2.0 <= dist <= distlim:
            # check if the required conditions are met
            # (minseqdist, intraonly, intraonly)
            if (i1 != i2 and (interonly or not intraonly)) or\
               (i1 == i2 and abs(j2 - j1) >= minseqdist and not interonly):
                # consider the pair a break
                perms.append((dssr1, dssr2, dist))

    return perms


def ApplyPerm(perm, mutants):
    """Applies the break to each of the mutants"""

    # mutants is a dictionary
    # keys = strand lists
    # values = performed breaks
    for strands, path in mutants.items():
        
        dssrloc = {}
        for i in range(len(strands)):
            for j in range(len(strands[i])):
                dssrloc[strands[i][j]] = (i, j)

        i1, j1 = dssrloc[perm[0]]
        i2, j2 = dssrloc[perm[1]]

        # take all the non-affected strands as is
        mutstrands = [strands[x] for x in range(len(strands)) if x != i1 and x != i2]

        # inter-strand break
        if i1 != i2:

            strand1 = strands[i1][:j1+1]+strands[i2][j2:]
            strand2 = strands[i2][:j2]+strands[i1][j1+1:]

        # intra-strand break
        else:

            # here the difference is in where the j1 and j2
            # will be included
            if j1 < j2:

                strand1 = strands[i1][:j1+1]+strands[i1][j2:]
                strand2 = strands[i1][j1+1:j2]

            else:

                strand1 = strands[i1][j2:j1+1]
                strand2 = strands[i1][:j2] + strands[i1][j1+1:]

        # add new strands if not empty
        if strand1:
            mutstrands.append(strand1)
        if strand2:
            mutstrands.append(strand2)

        mutstrands = tuple(sorted(mutstrands,key = lambda x: (x[0].split('.')[:2],
                                                   int(x[0].split('.')[3]),
                                                   x[0].split('.')[4])))
        # update the list of performed breaks
        mutpath = path + [perm,]

        yield mutstrands, mutpath       


def ApplyPerms(perms, strands, limit = 3, force = False, specperm = []):
    """Applies the breaks to the list of strands"""

    # initiate the dictionary of mutants with the original strands
    mutants = {strands:[]}

    # in case of specific user-provided permutation 
    if specperm:
        for perm in specperm:
            mutants = dict(ApplyPerm(perm, mutants))
        return mutants        

    # lim = allowed maximum number of breaks within a permutation
    lim = limit if limit >= 0 else len(perms)

    for k in range(lim):

        # copy mutants to newmutants
        newmutants = {k:v for k,v in mutants.items()}

        for perm in perms:

            # Apply each of the breaks to each of the current mutants
            for mutant, breaks in ApplyPerm(perm, mutants):

                # if the mutant is already present in newmutants
                # then it was already derived by a non-longer
                # series of breaks 
                if mutant not in newmutants:
                    newmutants[mutant] = breaks

            # monitor the size of newmutants and warn
            # in case of size > 10,000 and the "force" param is disabled
            if len(newmutants) > 10**4 and not force:
                print("!!WARNING!! More than 10,000 permutations obtained"+\
                      " ({}), I won't go further, check the force parameter.".format(len(newmutants)),
                      file=sys.stderr)
                return newmutants

        # update the mutants dictionary
        mutants = newmutants

    return mutants


def Continuous(dssr1, dssr2):
    """two dssrs are considered continuous
    if their models and chains are the same
    and their resnum values are (i, i+1)"""
    m1, ch1, r1, n1, i1 = dssr1.split('.')
    m2, ch2, r2, n2, i2 = dssr2.split('.')
    return m1==m2 and ch1==ch2 and int(n2) - int(n1) == 1


def StrandPretty(strand):
    """Print the strand in a compressed form.
    In the compressed form all series of continuous dssrs
    are replaced with a pair of the first and the last dssrs
    of the series, separated by a colon"""
    strand = [x for x in strand if x!='&']
    
    if len(strand) < 2:
        return ','.join(strand)

    res = []

    within = False

    for i in range(len(strand)):

        if 0 < i < len(strand)-1 and\
           Continuous(strand[i-1],strand[i]) and Continuous(strand[i],strand[i+1]):
            pass
        elif i > 0 and Continuous(strand[i-1],strand[i]):
            res[-1].append(strand[i])
        else:
            res.append([strand[i],])
            
    return ','.join([':'.join(x) for x in res])


def StrandToSeq(strand):
    """Returns a single-letter formatted sequence from a list of dssrs"""
    seq = [x.split('.')[2].upper() if x!='&' else '&' for x in strand]
    seq = [x if x in {'A','G','C','U','T','&'} else 'n' for x in seq]
    return ''.join(seq)


def CombinePairsToDBN(newpairs, length):
    """Converts a list of (i,j) base pairs into a string in dot-bracket notation"""
    dbn = ['.']*length

    levels = ['()', '[]', '{}', '<>', 'Aa', 'Bb', 'Cc', 'Dd', 'Ee', 'Ff', 'Gg',
              'Hh', 'Ii', 'Jj', 'Kk', 'Ll', 'Mm', 'Nn', 'Oo', 'Pp', 'Qq',
              'Rr', 'Ss', 'Tt', 'Uu', 'Vv', 'Ww', 'Xx', 'Yy', 'Zz']

    groups = [[],]

    seen = set()
    pairs = set()

    for p in newpairs:

        if p[0] not in seen and p[1] not in seen:
            pairs.add(p)
            seen.add(p[0])
            seen.add(p[1])
    
    for pair in sorted(pairs):

        level = 0

        while any(v[0]<=pair[0]<=v[1]<=pair[1] or pair[0]<=v[0]<=pair[1]<=v[1] for v in groups[level]):
            level += 1
            if level == len(groups):
                groups.append([])

        groups[level].append(pair)

    for times in range(len(groups)-1):
    
        for i in range(len(groups)-1):

            test = [v for v in groups[i] if any(v[0]<=w[0]<=v[1]<=w[1] or w[0]<=v[0]<=w[1]<=v[1] for w in groups[i+1])]

            if len(test) < len(groups[i+1]):

                groups[i]   = [p for p in groups[i] if p not in test] + groups[i+1]
                groups[i+1] = test

    for i,g in enumerate(groups):

        for p in g:
            dbn[p[0]] = levels[i][0]
            dbn[p[1]] = levels[i][1]
            
    return ''.join(dbn)


def PairsToDBN(strand, bps):
    """Returns a dot-bracket string for a strand"""
    if not bps:
        return ''

    dssrloc = {}
    pairs = []

    for i, x in enumerate(strand):
        dssrloc[x] = i

    for v, w in bps:
        if v in dssrloc and w in dssrloc:
            pairs.append((min(dssrloc[v],dssrloc[w]),max(dssrloc[v],dssrloc[w])))

    dbn = CombinePairsToDBN(pairs,len(strand))

    return ''.join([dbn[i] if strand[i]!='&' else '&' for i in range(len(strand))])            
    

def PrintStrands(strands, bps=None, num=None):
    """Print a mutant's strands"""
    thestrand = []

    for strand in strands:
        for res in strand:
            thestrand.append(res)
        thestrand.append('&')
    thestrand = thestrand[:-1]  

    if num == None:
        for i in range(len(strands)):
            print('strand {}:\t'.format(i+1), StrandPretty(strands[i]))
            print('   seq {}:\t'.format(i+1), StrandToSeq(strands[i]))
            if bps:
                print('   dbn {}:\t'.format(i+1), PairsToDBN(strands[i], bps))
        print('strand All:\t', StrandPretty(thestrand))
        print('   seq All:\t', StrandToSeq(thestrand))
        if bps:
            print('   dbn All:\t', PairsToDBN(thestrand, bps) )
    else:
        for i in range(len(strands)):
            print('strand {}.{}:\t'.format(num+1,i+1), StrandPretty(strands[i]))
            print('   seq {}.{}:\t'.format(num+1,i+1), StrandToSeq(strands[i]))
            if bps:
                print('   dbn {}.{}:\t'.format(num+1,i+1), PairsToDBN(strands[i], bps))
        print('strand {}.All:\t'.format(num+1), StrandPretty(thestrand))
        print('   seq {}.All:\t'.format(num+1), StrandToSeq(thestrand))
        if bps:
            print('   dbn {}.All:\t'.format(num+1), PairsToDBN(thestrand, bps) )
    

def PrintMutant(ind, perms, strands, bps):
    """Print a mutant's content"""
    print('------------------------------')
    print("Permutation {} ({} strands)".format(ind+1, len(strands)))
    print()
    print("Breaks({}):\n{}".format(len(perms),'\n'.join(['-'.join(x[:2]) for x in perms])))
    print()      
    PrintStrands(strands, bps, ind)


def ParseBasePairs(bpsfile, strands):
    """Parses a list of base pairs"""
    dssrs = set()
    bps = set()

    for strand in strands:
        for res in strand:
            dssrs.add(res)

    with open(bpsfile) as file:
        for line in file:
            # for DSSR-output files
            if "cWW  cW-W" in line and (" WC " in line or " WB " in line):
                dssr1, dssr2 = line.strip().split()[1:3]
                dssr1 = dssr1.replace('..','.')
                dssr2 = dssr2.replace('..','.')
                if dssr1 in dssrs and dssr2 in dssrs:
                    bps.add((dssr1, dssr2))
            # for an arbitrary list of base pairs
            else:
                linesplit = line.strip().split()
                if len(linesplit) > 1:
                    dssr1, dssr2 = linesplit[:2]
                    if dssr1 in dssrs and dssr2 in dssrs:
                        bps.add((dssr1, dssr2))
    return bps
                    
                
