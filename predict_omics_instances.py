######### Script for omics-wide prediction and functional enrichment analysis
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
from collections import Counter
import pandas as pd
import networkx as nx
import re
import textwrap
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import matplotlib as mpl
import reaction_nets as rn

def create_axes(x, y, ax, w, x_off=0, sc=1, for_icon=True):
    '''
    Create a new axes inside the given axes.
    x, y: the 'center' position of the new axes
    ax: the parent axes
    w: the width of the new axes
    w_off: the offset of the new axes from the center position
    sc: the scale of the new axes
    for_icon: whether the new axes is for an icon
    '''
    width = w
    # transform the data coordinates to display coordinates
    x_display, y_display = ax.transData.transform((x, y))
    # transform the display coordinates to figure-relative coordinates
    x_fig, y_fig = ax.figure.transFigure.inverted().transform((x_display, y_display))
    # create the new inset axes
    new_ax = ax.figure.add_axes([x_fig - (x_off+0.5)*width, y_fig - 0.5*width, width*sc, width*sc])
    new_ax.spines['left'].set_position('zero')
    new_ax.spines['bottom'].set_position('zero')
    new_ax.spines['right'].set_color('none')
    new_ax.spines['top'].set_color('none')
    new_ax.spines['left'].set_visible(False)
    new_ax.spines['bottom'].set_visible(False)
    new_ax.spines['right'].set_visible(False)
    new_ax.spines['top'].set_visible(False)
    new_ax.xaxis.set_ticks_position('bottom')
    new_ax.yaxis.set_ticks_position('left')
    if for_icon:
        new_ax.set_xticks([])
        new_ax.set_yticks([])
        new_ax.set_aspect('equal')
    return new_ax
def plot_hatched_circle(x, y, radius, ax, hp, c='k'):
    '''
    Plot a hatched circle
    x, y: the center position of the circle
    radius: the radius of the circle
    ax: the axes
    hp: the hatch pattern
    c: the color of the edge
    '''
    circle = Circle((x, y), radius, facecolor='none',
                     edgecolor=c, lw=1, zorder=-20)
    ax.add_patch(circle)
    if hp in ['w', 'k']:
        circle.set_facecolor(hp.replace('k', c).replace('w', 'none'))
    else:
        circle.set_hatch(hp)
def plot_polygon_vertices(n, ax, hp, scr=1, s=0.4, c='k'):
    '''
    Plot the vertices of a polygon as hatched circles
    n: the number of vertices
    ax: the axes
    hp: the hatch pattern
    scr: the scaling factor
    s: the size of the hatched circles (radius)
    c: the color of the edge
    '''
    if n == 1:
        plot_hatched_circle(0, 0, s, ax, hp=hp, c=c)
    else:
        # Calculate the angles between vertices
        if n == 4:
            angles = np.linspace(-np.pi/4, 7*np.pi/4, n+1)[:-1]
            r = 1/np.sqrt(2)
        elif n == 3:
            angles = np.linspace(-np.pi/6, 11*np.pi/6, n+1)[:-1]
            r = 1/np.sqrt(3)
        elif n == 2:
            angles = np.linspace(-np.pi/2*0, 4*np.pi/2, n+1)[:-1]
            r = 1/2
        # Calculate the x and y coordinates of the vertices
        x = np.cos(angles)*r*scr
        y = np.sin(angles)*r*scr
        for xn, yn, hp in zip(x, y, hp):
            #print(xn, yn)
            plot_hatched_circle(xn, yn, s, ax, hp=hp, c=c)
    # Set the axis limits
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    # Hide the axis labels and ticks
    ax.axis('off')
    # Show the plot
patterns = ['A', 'AA', 'AB', 'AAA', 'AAB', 'ABC', 'AAAA', 'AABB', 'AAAB', 'AABC', 'ABCD']
def id2xy(i, nodes):
    '''
    Convert the node id to the x and y coordinates
    '''
    for node in nodes:
        if node.id == i:
            return (node.x, node.y)
# Hatch patterns for A, B, C, D subunits in various complexes
hps = {'0_0': ['w'], '1_0':['w', 'w'], '1_1':['w', 'k'], '2_0':['w', 'w', 'w'], 
       '2_1': ['w', 'k', 'k'], '2_2': ['w', 'k', '//////'],
       '3_0': ['w']*4, '3_1': ['w']*2+['k']*2, '3_2':['w']+['k']*3,
       '3_3': ['w', 'k', 'k', '//////'], '3_4':['w', 'k', '//////', 'ooo']}


################## Make nodes #########################################

nodes = rn.make_nodes()

#######################################################################
############## Plot icons #############################################

fig, ax = plt.subplots(figsize=(6, 6))
ax.set_xlim(-2.2, 2.2)
ax.set_ylim(-3.1, 0.1)
n_nodes = [1, 2, 3, 5]
for i in range(len(n_nodes)):
    n_node = n_nodes[i]
    xs = np.linspace(-(n_node-1)/2, (n_node-1)/2, n_node)
    ax.scatter(xs, [-1*i]*n_node, s=1600, c='w', edgecolor='w', zorder=0)
links = [[0,1], [0,2], [1,3], [2,4], [2,5], [3,6], [4,7], [4,8], [4,9], [5,9], [5,10]]
for node in nodes:
    x, y = id2xy(node.id, nodes)
    ax1 = create_axes(x, y, ax, 0.08, x_off=0.0, sc=0.99)
    plot_polygon_vertices(node.iy+1, ax1, hp=hps[str(node.iy)+'_'+str(node.ix)])
for link in links:
    xa, ya = id2xy(link[0], nodes)
    xb, yb = id2xy(link[1], nodes)
ax.set_axis_off()
ax.set_ylabel('Number of binding events')
ax.set_xlabel('Number of molecular types')
ax.set_xticks([])
ax.set_yticks([])
#fig.savefig('./figures/comp_atlas_icons.svg', format='svg', dpi=600)
fig.savefig('./temp.png', format='png', dpi=600)
#plt.show()



#######################################################################
################### Plot paths ########################################
fig, axs = plt.subplots(figsize=(13, 8), ncols=6, nrows=4)
fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)
def plot_icons(ax, nodes=[], c='k'):
    ax.set_xlim(-2.2, 2.2)
    ax.set_ylim(-3.1, 0.1)
    n_nodes = [1, 2, 3, 5]
    for i in range(len(n_nodes)):
        n_node = n_nodes[i]
        xs = np.linspace(-(n_node-1)/2, (n_node-1)/2, n_node)
        ax.scatter(xs, [-1*i]*n_node, s=1600, c='w', edgecolor='w', zorder=0)
    for node in nodes:
        x, y = id2xy(node.id, nodes)
        ax1 = create_axes(x, y, ax, 0.08, x_off=0.0, sc=1.0)
        plot_polygon_vertices(node.iy+1, ax1, scr=0.7, s=0.3, hp=hps[str(node.iy)+'_'+str(node.ix)], c=c)
    ax.set_axis_off()
    ax.set_xticks([])
    ax.set_yticks([])
pos = dict(zip([node.id for node in nodes], [(node.x, node.y) for node in nodes]))
i = 0
for node in nodes[:]:
    for pathnet in node.pathnets:
        ax = axs.flat[i]
        plot_icons(ax, nodes=nodes, c='lightgray')
        plot_icons(ax, nodes=[n for n in nodes if n.id in pathnet.nodes])
        nx.set_node_attributes(pathnet, pos, name="pos")
        nx.draw_networkx_edges(pathnet, pos, node_size=400, node_shape='s',
                                edge_color='r', #zorder=10,
                               connectionstyle="arc3, rad=0.4", ax=ax)
        i += 1
#ax.spines['right', 'top', 'bottom', 'left'].set_visible(False)
fig.savefig('./temp.png', format='png', dpi=600)
#plt.show()


################################################################################
################## Check patterns in HuMap2 complexes ##########################

#prots_mem = pd.read_csv('protein_class_Predicted_membrane.tsv', sep='\t').Gene.values
prots_mem = pd.read_csv('https://www.proteinatlas.org/search/protein_class%3APredicted+membrane+proteins?format=tsv&download=yes', sep='\t').Gene.values
#prots_sec = pd.read_csv('protein_class_Predicted_secreted.tsv', sep='\t').Gene.values
prots_sec = pd.read_csv('https://www.proteinatlas.org/search/protein_class%3APredicted+secreted+proteins?format=tsv&download=yes', sep='\t').Gene.values
prots_mem_sec = set(prots_mem) | set(prots_sec)
#df = pd.read_csv('humap2_complexes_20200809.txt')
df = pd.read_csv('http://humap2.proteincomplexes.org/static/downloads/humap2/humap2_complexes_20200809.txt')
tcomplex, tcomplex_ext, tcomplex_sec = set(), set(), set()
allps = set()
psets = dict(zip(patterns, [set()]*len(patterns)))
complexes_ext = []
for us, gs in zip(df.Uniprot_ACCs, df.genenames):
    gnames = gs.split(' ')
    distgns = sorted(Counter(gnames).values())
    allps = allps.union(gnames)
    for p in patterns:
        psets[p] = psets[p].union(*rn.check1in2(p, gnames))
    if len(gnames) > 2:
        if len(set(gnames)) < len(gnames):
            pass
        tcomplex = tcomplex.union(gnames)
        if all([p in prots_mem_sec for p in gnames]):
            tcomplex_ext = tcomplex_ext.union(gnames)
            complexes_ext.append(gnames)
        if all([p in prots_sec for p in gnames]):
            tcomplex_sec = tcomplex_sec.union(gnames)
#print(len(tcomplex), len(tcomplex_ext), len(tcomplex_sec))
#print(len(allps))
for k, v in psets.items():
    rn.node_from_pattern(nodes, k).num_protein = len(v)
    print(k, len(v))
rn.node_from_pattern(nodes, 'A').num_protein = 20000
##############   String interactions #############################
# Download STRING data from https://string-db.org/cgi/download?sessionId=bpfpPhMIEEjk&species_text=Homo+sapiens
#dfprot = pd.read_csv('9606.protein.info.v12.0.txt', sep='\t')
dfprot = pd.read_csv('https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz', sep='\t', compression='gzip')
#dfs = pd.read_csv('./9606.protein.physical.links.full.v12.0.txt', sep=' ')
dfs = pd.read_csv('https://stringdb-downloads.org/download/protein.physical.links.full.v12.0/9606.protein.physical.links.full.v12.0.txt.gz', sep=' ', compression='gzip')

# Gene lists (including background) for enrichment analysis
#pd.Series(list(tcomplex_ext)).to_csv('./output_gene_lists_tcomplex_ext', index=False, header=False)
#pd.Series(list(prots_mem_sec)).to_csv('./output_gene_lists_prots_mem_sec_background', index=False, header=False)
if 0: # Read Enrichment results from GO website
    dfenr = pd.read_csv('./Enrichment_Ext_BgdExt.txt', sep='\t', skiprows=6) 
    dfenr.columns
    dfenr.iloc[:,0].str.contains('development|formation|Wnt|morphogen|epithel').sum()
    dfedp = dfenr[dfenr.iloc[:,0].str.contains('development|formation|Wnt|morphogen|epithel')]
    dfedp.iloc[:,0]
    print(dfenr.iloc[:,0].values)

#############################################################################
################## Check patterns in TargetScan complexes ########################
# Download from https://www.targetscan.org/vert_80/vert_80_data_download/Summary_Counts.default_predictions.txt.zip
#df = pd.read_csv('Summary_Counts.default_predictions.txt', sep='\t')
dfts = pd.read_csv('https://www.targetscan.org/vert_80/vert_80_data_download/Summary_Counts.default_predictions.txt.zip', sep='\t', compression='zip')
dfh = dfts[dfts['Species ID'] == 9606]
dfhg = dfh.groupby('Gene Symbol').agg('sum')
dfhg[(dfhg['Total num conserved sites'])> 1].shape
targ_2sites_list = dfhg[(dfhg['Total num conserved sites'])> 1].index.tolist()
targ_3sites_list = dfhg[(dfhg['Total num conserved sites'])> 2].index.tolist()
rn.node_from_pattern(nodes, 'A').num_mmi = 20000
for p in ['AA', 'AAA', 'AAAA', 'AABB']:
    rn.node_from_pattern(nodes, p).num_mmi = 0
rn.node_from_pattern(nodes, 'AB').num_mmi = len(set(dfh['Gene Symbol']))
rn.node_from_pattern(nodes, 'ABC').num_mmi = len(set(targ_2sites_list))
rn.node_from_pattern(nodes, 'AAB').num_mmi = len(set(dfh[dfh['Total num conserved sites'] >1].index))
rn.node_from_pattern(nodes, 'AAAB').num_mmi = len(set(dfh[dfh['Total num conserved sites'] >2].index))
rn.node_from_pattern(nodes, 'ABCD').num_mmi = len(set(targ_2sites_list))
gene_mi_df = dfh.groupby(['Gene Symbol', 'Representative miRNA'], as_index=False).agg('count')
gene_mi_df['ref_count'] = 1
gene_mi_df = gene_mi_df.groupby(['Gene Symbol']).agg('count')
len(set(gene_mi_df[gene_mi_df['Representative miRNA']>1].index).intersection(set(targ_2sites_list)))
rn.node_from_pattern(nodes, 'AABC').num_mmi = len(set(gene_mi_df[gene_mi_df['Representative miRNA']>1].index).intersection(set(targ_2sites_list)))
# Pattern-enabling mRNAs based on mRNA-miRNA interactions
pe_complex_mRNAs = set(targ_2sites_list) | set(dfh[dfh['Total num conserved sites'] >1]['Gene Symbol']) | set(dfh[dfh['Total num conserved sites'] >2]['Gene Symbol']) | set(gene_mi_df[gene_mi_df['Representative miRNA']>1].index).intersection(set(targ_2sites_list))

# For enrichment analysis
#pd.Series(list(pe_complex_mRNAs)).to_csv('./output_gene_lists_pe_complex_mRNAs', index=False, header=False)
if 0:
    dfenrm = pd.read_csv('./Enrichment_mRNA_BgdAll.txt', sep='\t', skiprows=5) 
    dfenrm.columns
    dfenrm.iloc[:,0].str.contains('development|formation|Wnt|morphogen|epithel').sum()
    dfedm = dfenrm[dfenrm.iloc[:,0].str.contains('development|formation|Wnt|morphogen|epithel|vesicle')]
    dfedm.iloc[:,0]
    print(dfenrm.iloc[:,0].values)
    dfedp.iloc[:,[-4, -2]] = dfedp.iloc[:,[-4, -2]].astype(float)
    dfedm.iloc[:,[-4, -2]] = dfedm.iloc[:,[-4, -2]].astype(float)
    dfedp.loc[:,'name'] = dfedp.iloc[:,0].str.extract(r'(.*) \(')[0].values
    dfedm.loc[:, 'name'] = dfedm.iloc[:,0].str.extract(r'(.*) \(')[0].values
    fig, ax = plt.subplots(figsize=(3, 6))
    fig.subplots_adjust(left=.65, top=0.98, bottom=0.07, right=0.99)
    nls = []
    for l in list(dfedp.iloc[:,-1].values):
        nl = textwrap.fill(l, 30, max_lines=2)
        print(nl)
        nl = nl[0].upper() + nl[1:]
        nls.append(nl)
    # Change the colors of the bars based on pvs
    pvs = dfedp.iloc[:,-2].values
    my_cmap = cm.get_cmap('winter')
    my_norm = Normalize(vmin=pvs.min(), vmax=pvs.max())
    ax.barh(nls, dfedp.iloc[:,-4].values, color=my_cmap(my_norm(pvs)), height=0.05)
    sizes = np.log(dfedp.iloc[:,2]).values
    sizes = sizes - (sizes.min()-0.5)
    ax.scatter(dfedp.iloc[:,-4].values, nls, s=10*sizes, color=my_cmap(my_norm(pvs)))
    ax.tick_params(axis='both', which='major', labelsize=8)
    print(nls)
    ax.set_xlim(0, dfedp.iloc[:,-4].max()+0.5)
    ax.set_xlabel('Fold enrichment', size=9)
    ax.set_ylim(-0.5, len(nls)-0.5)
    cax = fig.add_axes([0.76, 0.79, 0.03, 0.18])
    cb1 = mpl.colorbar.ColorbarBase(cax, cmap=my_cmap,
                                    norm=my_norm,
                                    orientation='vertical')
    cb1.set_ticks([0.01, 0.02, 0.03, 0.04])
    cb1.set_label('FDR')
    #dfedm.iloc[:, [-1, -4, -2]].rename(columns={dfedm.columns[-1]: 'GO term', dfedm.columns[-4]: 'Fold enrichment', dfedm.columns[-2]: 'FDR'}).to_csv('./output_gene_lists/GO_enrichment_mRNA_BgdAll', index=False, header=True, sep='\t')


#################### Plot instances with complexes ################################

fig, ax = plt.subplots(figsize=(6, 6))
fig.subplots_adjust(bottom=0.125)
ax.set_xlim(-2.2, 2.2)
ax.set_ylim(-3.1, 0.1)
n_nodes = [1, 2, 3, 5]
for i in range(len(n_nodes)):
    n_node = n_nodes[i]
    xs = np.linspace(-(n_node-1)/2, (n_node-1)/2, n_node)
    ax.scatter(xs, [-1*i]*n_node, s=1600, c='w', edgecolor='w', zorder=0)
links = [[0,1], [0,2], [1,3], [2,4], [2,5], [3,6], [4,7], [4,8], [4,9], [5,9], [5,10]]
for node in nodes[1:]:
    x, y = id2xy(node.id, nodes)
    y = y*0.8
    ax1 = create_axes(x, y, ax, 0.08, x_off=0.5, sc=0.80)
    plot_polygon_vertices(node.iy+1, ax1, hp=hps[str(node.iy)+'_'+str(node.ix)])
    ax2 = create_axes(x, y, ax, 0.08, x_off=-0.3, sc=0.95, for_icon=False)
    ax2.bar([0,], [1,], color='none', edgecolor='indianred')
    ax2.bar([0,], [node.num_mmi/20000,], color="indianred")
    ax2.bar([1,], [1,], color='none', edgecolor='b')
    ax2.bar([1,], [node.num_protein/20000,], color="b")
    ax2.text(1.4, 1, r'2$\times$10$^4$', size=8, ha='left', va='center')
    ax2.set_xlim(-0.5,2)
    ax2.set_yticks([])
    ax2.set_ylim(-0.1,1.1)
    ax2.set_xticks([0, 1])
    ax2.set_xticklabels([str(node.num_mmi), str(node.num_protein)], rotation=45, ha='right')
    for xtick, color in zip(ax2.get_xticklabels(), ['indianred', 'b']):
        xtick.set_color(color)
for link in links:
    xa, ya = id2xy(link[0], nodes)
    xb, yb = id2xy(link[1], nodes)
ax.spines[['right', 'top', 'bottom', 'left']].set_visible(False)
ax.set_ylabel('Comlexity of multimer               ', labelpad=20)
ax.set_xlabel('Number of molecular types', labelpad=-20)
ax.set_xticks([])
ax.set_yticks([])
plt.show()


#################### Get representative complex for each gene ################################
rep_complexes = {}
for g in tcomplex | set(dfh['Gene Symbol']):
    rep_complexes[g] = {'protein':'', 'RNA':''}
for us, gs in zip(df.Uniprot_ACCs, df.genenames):
    gnames = gs.split(' ')
    if len(gnames)>2 and all([g in tcomplex for g in gnames]):
        for g in gnames:
            if rep_complexes[g]['protein'] == '':
                rep_complexes[g]['protein'] = gs.replace(' ', ',')
dfmireps = dfh.loc[:,['Representative miRNA', 'Gene Symbol']].groupby('Gene Symbol').agg(','.join)
for g in dfmireps[dfmireps['Representative miRNA'].str.contains(',')].index:
    if rep_complexes[g]['RNA'] == '':
        rep_complexes[g]['RNA'] = g+','+dfmireps.loc[g, 'Representative miRNA']
i = 0
for g in rep_complexes:
    if rep_complexes[g]['protein'] != '':
        i += 1
print('protein', i)
i = 0
for g in rep_complexes:
    if rep_complexes[g]['RNA'] != '':
        i += 1
print('RNA', i)
with open('./output_gene_lists_TableS5.txt', 'w') as f:
    f.write('Gene\tRepresentative protein complex\tRepresentative RNA complex\n')
    for g in rep_complexes:
        f.write(g)
        f.write('\t')
        f.write(rep_complexes[g]['protein'])
        f.write('\t')
        f.write(rep_complexes[g]['RNA'])
        f.write('\n')