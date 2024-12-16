import argparse
parser=argparse.ArgumentParser(prog='SAMap_plot', description='Visualize the SAMap result!')
parser.add_argument("--name", type=str, help="The name abbreviation of the species. eg:DR_SE_OM_CP_AS_MA_PA")
parser.add_argument('--alignThr', action='store', help='Mapping score threshold.')
parser.add_argument("--path1", type=str, help="The path of the pkl file.")
parser.add_argument("--path2", type=str, help="The path of the result.")
args=parser.parse_args()

import os
import datetime
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder, sankey_plot, chord_plot, CellTypeTriangles)
from samalg import SAM
from samap.utils import save_samap, load_samap
import holoviews as hv

def make_keys_plots(name, path):
  keys={}
  name_list=name.split("_")
  for i in name_list:
    keys[i]="celltype"
    sm.sams[i].adata.obsm['X_umap']=sm.sams[i].adata.obsm['X_umap_samap']
    plt.figure(figsize=(100,150))
    sc.pl.umap(sm.sams[i].adata, color='celltype')
    out_umap_file=os.path.join(path, f"{i.strip()}_{today}_umap.png")
    plt.savefig(out_umap_file, bbox_inches='tight', dpi=600)
  return keys

def q(x):
    return np.array(list(x))

today=datetime.datetime.now().strftime("%Y%m%d")
sm=load_samap(args.path1)
keys=make_keys_plots(args.name, args.path2)
species=args.name

D,MappingTable=get_mapping_scores(sm, keys, n_top = 0)
out_file_1=os.path.join(args.path2, f"{args.name.strip()}_{today}_MappingTable.csv")
MappingTable.to_csv(out_file_1, sep='\t', index=True)

M=MappingTable
species_order=species.split("_")
align_thr=float(args.alignThr)

ids = np.array(species_order)
d = M.values.copy()
d[d<align_thr]=0
x,y = d.nonzero()
x,y = np.unique(np.sort(np.vstack((x,y)).T,axis=1),axis=0).T
values = d[x,y]
nodes = q(M.index)
node_pairs = nodes[np.vstack((x,y)).T]

sn1 = q([xi.split('_')[0] for xi in node_pairs[:,0]])
sn2 = q([xi.split('_')[0] for xi in node_pairs[:,1]])

filt_bak=[0]*len(ids)
for i in range(len(ids)-1):
        filt_bak[i]=np.logical_or(np.logical_and(sn1==ids[i],sn2==ids[i+1]),np.logical_and(sn1==ids[i+1],sn2==ids[i]))
filt=(sum(elem for elem in filt_bak))>0
x,y,values=x[filt],y[filt],values[filt]

d=dict(zip(ids,list(np.arange(len(ids)))))
depth_map = dict(zip(nodes,[d[xi.split('_')[0]] for xi in nodes]))
data =  nodes[np.vstack((x,y))].T
for i in range(data.shape[0]):
    if d[data[i,0].split('_')[0]] > d[data[i,1].split('_')[0]]:
        data[i,:]=data[i,::-1]
R = pd.DataFrame(data = data,columns=['source','target'])
R['Value'] = values
R1=R

def f(plot,element):
    plot.handles['plot'].sizing_mode='scale_width'
    plot.handles['plot'].x_range.start = -600
    plot.handles['plot'].x_range.end = 1500

try:
        from holoviews import dim
        import holoviews as hv
        hv.extension('bokeh',logo=False)
        hv.output(size=100)
except:
        raise ImportError('Please install holoviews-samap with `!pip install holoviews-samap`.')

sankey1 = hv.Sankey(R1, kdims=["source", "target"]).options(show_values=False)
sankey1.opts(cmap='Colorblind',label_position='outer', edge_line_width=0, show_values=False,
             node_padding=4,depth_map=depth_map, node_alpha=1.0, node_width=20,
             node_sort=True,frame_height=1000,frame_width=1200,bgcolor='snow',
             apply_ranges=True,hooks=[f])
sankey_path=os.path.join(args.path2, f"{args.name}_{args.alignThr}.sankey.html")
hv.save(sankey1, filename=sankey_path, fmt='html')

gpf=GenePairFinder(sm,keys=keys)
gene_pairs=gpf.find_all(thr=0.2)
out_file_2=os.path.join(args.path2, f"{args.name.strip()}_{today}_genepairs.csv")
gene_pairs.to_csv(out_file_2, sep='\t', index=True)

plt.figure(figsize=(100,150))
sc.pl.umap(sm.samap.adata, color='species')
out_file_3=os.path.join(args.path2, f"{args.name.strip()}_{today}_umap.png")
plt.savefig(out_file_3, bbox_inches='tight', dpi=600)
