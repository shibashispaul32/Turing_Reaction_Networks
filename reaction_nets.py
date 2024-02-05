import numpy as np
from scipy.special import comb
import re
from itertools import product, combinations, chain, permutations
from collections import Counter
import networkx as nx
import pandas as pd

def checkeq(l1, l2, if_map=False):
    '''
    Check if two iterables are equivalent
    if_map: if True, return the mapping
    '''
    if len(l1) != len(l2):
        return False
    d1 = dict(sorted(Counter(l1).items(), key=lambda item: item[1]))
    d2 = dict(sorted(Counter(l2).items(), key=lambda item: item[1]))
    mdict = {}
    for k1, k2 in zip(d1, d2):
        a = d1[k1]
        b = d2[k2]
        if a != b:
            return False
        else:
            mdict[k1] = k2
    if if_map == True:
        return True, mdict
    else:
        return True
def checkpermute(l1, l2):
    '''
    Check if two iterables are permutations of each other
    '''
    return Counter(l1) == Counter(l2)
def check1in2(l1, l2, if_permute=False):
    '''
    Check if one iterable is a subset of another. If if_permute is True, consider permutations
    '''
    subsets = []
    for subset in list(combinations(l2, len(l1))):
        if if_permute:
            if checkpermute(subset, l1):
                subsets.append(subset)
        else:
            if checkeq(subset, l1):
                subsets.append(subset)
    return set(subsets)
class Node():
    def __init__(self, x, y, ix, iy, i):
        '''
        A Node object holds information about a characteristic complex
        '''
        self.x = x # x: x-coordinate of the node
        self.y = y # y: y-coordinate of the node
        self.ix = ix # ix: index of x-coordinate of the node
        self.iy = iy # iy: index of y-coordinate of the node
        self.id = i # id: index of the node
        self.paths = [] # a list of reaction networks that yield the complex
        self.pathnets = [] # a networkx graph of reaction networks. Complexes with equivalent configurations are merged
        self.modnets = [] # a networkx graph of reaction networks. Complexes with equivalent configurations are separated
        self.num_mmi = 0 # Instances of mRNA-microRNA systems
        self.num_protein = 0 # Instances of protein systems
        self.num_bifs = [] # number of bifurcations
        self.models = [] # ODE model objects
def node_from_pattern(nodes, pattern):
    '''
    Return a Node object from a list of Node objects based on a pattern
    '''
    for node in nodes:
        if checkeq(node.pattern, pattern):
            return node
def enumerate_FFL(G):
    '''
    Return FFL motifs from a directed graph
    '''
    # Create an empty directed graph to hold the motif instances
    motifs = nx.DiGraph()
    # Loop over all nodes in the graph
    for node in G.nodes():
        # Get the neighbors of the current node
        neighbors = list(G.neighbors(node))
        # If the current node has at least two neighbors...
        if len(neighbors) >= 2:
            # Loop over all pairs of neighbors
            for i in range(len(neighbors)):
                for j in range(len(neighbors)):
                    # If there is an edge from the first neighbor to the second neighbor...
                    if G.has_edge(neighbors[i], neighbors[j]):
                        # Add the A->B, A->C->B motif to the motifs graph
                        motifs.add_edges_from([(node, neighbors[j]), (node, neighbors[i]), (neighbors[i], neighbors[j])])
    if len(motifs) > 0:
        if not nx.is_weakly_connected(motifs):
            # TODO: split components using nx.weakly_connected_components
            print("Multiple FFL modules.")
    return motifs
def add_deg_to_df(df, reactant, products, subunit_name=""):
    '''
    Add degradation reactions to a dataframe
    '''
    reaction = reactant+"->"+"+".join(products)
    parname = 'k_' + subunit_name + '*' + subunit_name.lower()+'_'+reactant
    niso = reactant.count(subunit_name)
    nselfiso = df[df.Species==reactant].SelfIso.unique()[0]
    df.loc[len(df.index)] = [reactant, 'Reactant', 1, niso*nselfiso, nselfiso, parname, reaction, 'F']
    for product in products:
        nselfiso = df[df.Species==product].SelfIso.unique()[0]
        df.loc[len(df.index)] = [product, 'Product', 1, niso*nselfiso, nselfiso, parname, reaction, 'F']
def degrade1mol(names_in_path, name):
    all_frags = [] # Each element is [Subunit to remove, [Recycled complexes/subunits]]
    for i in range(len(name)):
        rm_mol = name[i]
        newname = name[:i] + name[i+1:]
        frags = []
        for n in reversed(range(2, len(newname)+1)):
            for k in combinations(range(len(newname)), r=n):
                frag = ''.join([newname[l] for l in k])
                if_found_match = False
                for name_in_path in names_in_path:
                    b = checkpermute(name_in_path, frag)
                    if b:
                        if_found_match = True
                        rest_frag = ''.join([newname[l] for l in range(len(newname)) if not l in k])
                        frags.append(frag)
                        break
                if if_found_match:
                    newname = rest_frag
                    break
        all_frags.append([rm_mol, [l for l in newname]+frags])
    return all_frags
def make_nodes():
    '''
    Make a list of Node objects, each of which represents a characteristic complex
    Find reaction networks that yield each complex
    Store the reaction networks in each Node
    Return a list of Node objects
    '''
    patterns = ['A', 'AA', 'AB', 'AAA', 'AAB', 'ABC', 'AAAA', 'AABB', 'AAAB', 'AABC', 'ABCD']
    nodes = []
    n_nodes = [1, 2, 3, 5]
    for i in range(len(n_nodes)):
        n_node = n_nodes[i]
        xs = np.linspace(-(n_node-1)/2, (n_node-1)/2, n_node)
        for ix, x in enumerate(xs):
            node = Node(x, -i, ix, i, len(nodes))
            node.pattern = patterns[len(nodes)]
            nodes.append(node)
    def add_binding_to_df(df, reactants, product, reaction=None, if_rev=True):
        df.loc[len(df.index)] = [reactants[0], 'Reactant', 1, 1, 1, 'K*koff', reaction, 'F']
        df.loc[len(df.index)] = [reactants[1], 'Reactant', 1, 1, 1, 'K*koff', reaction, 'F']
        df.loc[len(df.index)] = [product, 'Product', 1, 1, 1, 'K*koff', reaction, 'F']
        if if_rev:
            df.loc[len(df.index)] = [product, 'Reactant', 1, 1, 1, 'koff', reaction, 'B']
            df.loc[len(df.index)] = [reactants[0], 'Product', 1, 1, 1, 'koff', reaction, 'B']
            df.loc[len(df.index)] = [reactants[1], 'Product', 1, 1, 1, 'koff', reaction, 'B']
        return df
    genpathid = 0
    for node in nodes:
        node.paths = []
        node.pathnets = []
        node.modnets = []
    for node in nodes[:]:
        if node.iy == 0:
            node.paths = []
            node.pathnets = [nx.DiGraph()]
            node.modnets = [nx.DiGraph()]
            node.paths = []
        elif node.iy == 1:
            df = pd.DataFrame({'Species':[], 'Type':[], 'Stoic':[], 'Iso':[], 'SelfIso':[], 'Par':[], 'Reaction':[], 'Dirc':[]})
            reaction_str = node.pattern[0]+'+'+node.pattern[1]+'<->'+node.pattern
            df = add_binding_to_df(df, [node.pattern[0],node.pattern[1]], node.pattern, reaction=reaction_str)
            node.paths = [df]
            g = nx.DiGraph()
            g.add_nodes_from([0, node.id])
            g.add_edge(0,node.id)
            node.pathnets.append(g)
            mg = nx.DiGraph()
            mg.add_edges_from([(node.pattern[0], node.pattern), (node.pattern[1], node.pattern)])
            node.modnets.append(mg)
        else:
            for i in range(1, node.iy):
                nodes_i = [n for n in nodes if n.iy == i]
                for node_i in nodes_i:
                    s = check1in2(node_i.pattern, node.pattern, if_permute=True)
                    patternl = list(node.pattern)
                    if len(s) > 0:
                        s1 = list(s)[0]
                        rm_ids = []
                        for a in s1:
                            rm_ids.append(patternl.index(a))
                        for index in sorted(rm_ids, reverse=True):
                            del patternl[index]
                        pnode1 = node_from_pattern(nodes, s1)
                        pnode2 = node_from_pattern(nodes, patternl)
                        b1, m1 = checkeq(pnode1.pattern, s1, if_map=True)
                        b2, m2 = checkeq(pnode2.pattern, patternl, if_map=True)
                        trantab_m1 = str.maketrans(''.join(m1.keys()), ''.join(m1.values()))
                        trantab_m2 = str.maketrans(''.join(m2.keys()), ''.join(m2.values()))
                        mod_pattern1 = ''.join(s1)
                        mod_pattern2 = ''.join(patternl)
                        c1paths, c2paths = [], []
                        if len(pnode1.paths) > 0:
                            for path in pnode1.paths:
                                modpath = path.copy()
                                modpath['Species'] = modpath['Species'].str.translate(trantab_m1)
                                modpath['Reaction'] = modpath['Reaction'].str.translate(trantab_m1)
                                c1paths.append(modpath)
                        else:
                            c1paths =  [pd.DataFrame().reindex_like(pnode2.paths[0]).dropna()]
                        if len(pnode2.paths) > 0:
                            for path in pnode2.paths:
                                modpath = path.copy()
                                modpath['Species'] = modpath['Species'].str.translate(trantab_m2)
                                modpath['Reaction'] = modpath['Reaction'].str.translate(trantab_m2)
                                c2paths.append(modpath)
                        else:
                            c2paths =  [pd.DataFrame().reindex_like(pnode1.paths[0]).dropna()]
                        j = 0
                        curr_lim = len(node.pathnets)
                        for (i1, i2) in product(range(len(c1paths)), range(len(c2paths))):
                            c1p, c2p = c1paths[i1], c2paths[i2]
                            if_same_pnodes = False
                            for g in node.pathnets[:curr_lim]:
                                if g.has_edge(pnode1.id, node.id) and g.has_edge(pnode2.id, node.id):
                                    if_same_pnodes = True
                                    continue
                            if if_same_pnodes == True:
                                continue
                            genpathid += 1
                            dfcomb = pd.concat([c1p, c2p])
                            reaction_str = mod_pattern1+'+'+mod_pattern2+'<->'+node.pattern
                            dfcomb = add_binding_to_df(dfcomb, [mod_pattern1, mod_pattern2], node.pattern, reaction=reaction_str)
                            node.paths.append(dfcomb)
                            g1, g2 = pnode1.pathnets[i1], pnode2.pathnets[i2]
                            gcomb = nx.compose(g1, g2)
                            gcomb.add_edges_from([(pnode1.id, node.id), (pnode2.id, node.id)])
                            node.pathnets.append(gcomb)
                            mg = nx.DiGraph()
                            for r in set(dfcomb.Reaction):
                                m = re.match(r'(\w+)\+(\w+)\<\-\>(\w+)', r)
                                mg.add_edges_from([(m.group(1), m.group(3)), (m.group(2), m.group(3))])
                            node.modnets.append(mg)
                            j += 1
                            # Check symmetrical binding via scafolding
                            a = dfcomb.Species.unique()
                            rev = list(permutations(a, 3))
                            sca_lv = 0
                            for c in list(rev):
                                if len(c) != 3:
                                    continue
                                c1, c2, c3 = c
                                if len(c1)==3 and len(c2)==2 and len(c3)==1 and not c2[0]==c2[1]:
                                    if c2 in c1 and not c3 in c2 and c3 in c1:
                                        sca_lv = 2
                                        break
                            if sca_lv == 2:
                                for d in list(rev):
                                    if len(d) != 3:
                                        continue
                                    d1, d2, d3 = d
                                    if len(d1)==4 and d2==c1 and len(d3)==1:
                                        if d2 in d1 and not d3 in d2 and d3 in d1:
                                            sca_lv = 3
                                            break
                            if sca_lv in [2, 3]:
                                if not c2[1]+c2[1] in a:
                                    sca = c2[1]
                                    su1, su2, trimer = c2[0], c3, c1
                                else:
                                    print('Error in finding scafold.')
                                dfcomb_new = add_binding_to_df(dfcomb.copy(), [sca, su2], sca+su2, reaction=sca+'+'+su2+'<->'+sca+su2)
                                dfcomb_new = add_binding_to_df(dfcomb_new, [sca+su2, su1], trimer, reaction=sca+su2+'+'+su1+'<->'+trimer)
                            #if sca_lv == 3:
                            if 0:
                                su3, tetramer = d3, d1
                                dimer1, dimer2, dimer3 = c2, sca+su2, sca+su3 # BC and BD
                                trimer2 = dimer2+su3 # BCD
                                trimer3 = su1+sca+su3 # ABD
                                tetramer = d1 # ABCD
                                dfcomb_new = add_binding_to_df(dfcomb_new, [sca, su3], dimer3, reaction=sca+'+'+su3+'<->'+dimer3)
                                dfcomb_new = add_binding_to_df(dfcomb_new, [dimer2, su3], trimer2, reaction=dimer2+'+'+su3+'<->'+trimer2)
                                dfcomb_new = add_binding_to_df(dfcomb_new, [dimer3, su2], trimer2, reaction=dimer3+'+'+su2+'<->'+trimer2)
                                dfcomb_new = add_binding_to_df(dfcomb_new, [dimer1, su3], trimer3, reaction=dimer1+'+'+su3+'<->'+trimer3)
                                dfcomb_new = add_binding_to_df(dfcomb_new, [dimer3, su1], trimer3, reaction=dimer3+'+'+su1+'<->'+trimer3)
                                dfcomb_new = add_binding_to_df(dfcomb_new, [trimer2, su1], tetramer, reaction=trimer2+'+'+su1+'<->'+tetramer)
                                dfcomb_new = add_binding_to_df(dfcomb_new, [trimer3, su2], tetramer, reaction=trimer3+'+'+su2+'<->'+tetramer)
                            if sca_lv > 0:
                                node.paths[-1] = dfcomb_new
                                mg = nx.DiGraph()
                                for r in set(dfcomb_new.Reaction):
                                    m = re.match(r'(\w+)\+(\w+)\<\-\>(\w+)', r)
                                    mg.add_edges_from([(m.group(1), m.group(3)), (m.group(2), m.group(3))])
                                node.modnets[-1] = mg
                            #print(node.pattern, len(node.paths), len(node.paths[-1].Species.unique()))
                    else:
                        continue
    for n in nodes:
        for i, net in enumerate(n.modnets):
            fflnodes = [x[1] for x in sorted((len(a), a) for a in enumerate_FFL(net).nodes)]
            if 0:
                if i == 3 and n.pattern=='AABC':
                    print('FFL nodes', fflnodes)
                elif i == 3 and n.pattern=='AAAB':
                    print('FFL nodes', fflnodes)
                elif i == 0 and n.pattern=='ABC':
                    print('FFL nodes', fflnodes)
                else:
                    continue
                print(n.pattern, fflnodes)
            nbs = max(len(fflnodes)-1, 1)
            if nbs > 1:
                ib = 0 # 
                nodelendiffs = set([len(b)-len(a) for a, b in zip(fflnodes[:-1], fflnodes[1:])])
                if 0 in nodelendiffs:
                    continue
                for a, b in zip(fflnodes[:-1], fflnodes[1:]):
                    print(n.pattern, i,  a, b, len(a), len(b))
                    df = n.paths[i]
                    nforms_r = comb(nbs-ib, 1)*comb(nbs, nbs-ib)
                    df.loc[df.Reaction.str.contains(fr'\b{a}\b.*\-\>\b{b}\b'), 'Iso'] = nforms_r
                    nforms_s = comb(nbs, ib)
                    df.loc[df.Species==a, 'SelfIso'] = nforms_s
                    ib += 1
    for node in nodes[:]:
        if len(node.paths) == 0:
            continue
        if len(node.paths[0]) == 0:
            continue
        for path in node.paths:
            names_in_path = [n for n in set(path.Species) if len(n)>1]
            for name in names_in_path:
                all_frags = degrade1mol(names_in_path, name)
                for frags in all_frags:
                    add_deg_to_df(path, name, frags[1], subunit_name=frags[0])
    for n in nodes: # Combine redundant rows and update stoic
        if len(node.paths) == 0:
            continue
        if len(node.paths[0]) == 0:
            continue
        for i, df in enumerate(n.paths):
            mdf = df.groupby(['Species', 'Type', 'Par', 'Reaction', 'Dirc', 'Iso', 'SelfIso'], sort=False, as_index=False)["Stoic"].count().reset_index().drop('index', axis=1)
            mdf['TrueStoic'] = [len(re.findall(r'\b'+y+r'\b', x)) for x, y in zip(mdf['Reaction'], mdf['Species'])]
            n.paths[i] = mdf
    # Diversify Ks (dissociation constants)
    for n in nodes:
        for i, df in enumerate(n.paths):
            m = df.Reaction.str.extract(r'\<\-\>(\w+)')
            newpar = np.where(m[0].isna(), df.Par, [re.sub('K', 'K_'+str(m[0].iloc[x]), df.Par.iloc[x]) for x in df.Par.index])
            df['Par'] = newpar
    return nodes

def make_odes(df):
    '''
    Make ODEs based on reactions in a dataframe
    '''
    vars, pars = {}, []
    for v in df.Species.unique():
        vars[v] = ''
    s = set(df.Species.sum())
    if len(s) > 0:
        for svar in s:
            vars[svar] = 's'+svar
            vars[svar] += '-k_'+svar+'*'+svar
            pars.extend(['s'+svar, 'k_'+svar])
    dfrd = df.groupby(['Reaction', 'Dirc'], as_index=False).agg('count')[['Reaction', 'Dirc']]
    for ird in range(len(dfrd)):
        reaction, dirc = dfrd.iloc[ird, 0], dfrd.iloc[ird, 1]
        rrows = df[(df.Reaction==reaction)&(df.Dirc==dirc)]
        crows = rrows[rrows.Type=='Reactant']
        prows = rrows[rrows.Type=='Product']
        if crows.shape[0] == 2:
            ra, rb = crows.Species
            tsa, tsb = crows.TrueStoic
            base_rate = "*".join([ra]*tsa) +"*"+ "*".join([rb]*tsb)
        else:
            ra = crows.Species.iloc[0]
            tsa = crows.TrueStoic.iloc[0]
            base_rate = "*".join([ra]*tsa)
        for r in crows.Species.unique():
            iso, siso, p = crows[crows.Species==r].loc[:, ['Iso', 'SelfIso', 'Par']].iloc[0,:]
            tiso = str(int(iso)//int(siso))
            if crows.shape[0] == 1:
                tiso = tiso + '*' + str(crows.TrueStoic.iloc[0])
            vars[r] += '-' + p + re.sub(r"\*1", "", '*' + tiso) + '*' + base_rate
        for r in prows.Species.unique():
            siso, p = prows[prows.Species==r].loc[:, ['SelfIso', 'Par']].iloc[0,:]
            tiso = str(int(iso)//int(siso)) + '*' + str(prows[prows.Species==r].TrueStoic.iloc[0])
            vars[r] += '+' + p + re.sub(r"\*1", "", '*' + tiso) + '*' + base_rate
    pars += sorted(list({m.group(0) for m in re.finditer(r"\b(\w+)\b", "*".join(df.Par.unique()))}))
    return vars, pars

if __name__ == "__main__":
    nodes = make_nodes()
    for n in nodes:
        for ip, p in enumerate(n.paths[::]):
            vs, pars = make_odes(p)
            print(vs)
            print(pars)
