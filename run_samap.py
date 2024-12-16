#!/dellfsqd2/ST_OCEAN/USER/yangqi3/01_Software/mamba/mambaforge/envs/samap/bin/python
# -*- coding: utf-8 -*-
import argparse
parser = argparse.ArgumentParser(prog = 'SAMap_pipeline', description = 'Run the SAMap algorithms!')
parser.add_argument("--name", type=str, help="The name abbreviation of the species. eg:DR_SE_OM_CP_AS_MA_PA")
parser.add_argument("--path1", type=str, help="The path of the h5ad file. eg:/dellfsqd2/ST_OCEAN/USER/yangqi3/02_Project/02_nervous_system/02.cross-species/00.data/03.h5ad/")
parser.add_argument("--path2", type=str, help="The path of the blast result. eg:/dellfsqd2/ST_OCEAN/USER/yangqi3/02_Project/02_nervous_system/02.cross-species/02.samap/01.blast/03.final_maps/")
parser.add_argument("--path3", type=str, help="The path of the result. eg:/dellfsqd2/ST_OCEAN/USER/yangqi3/02_Project/02_nervous_system/02.cross-species/02.samap/03.result/01.pkl/")
args = parser.parse_args()

import os
import datetime
from samap.mapping import SAMAP
from samalg import SAM
import pandas as pd
from samap.utils import save_samap

def create_sams(name, path):
  sams={}
  name_list=name.split("_")
  for i in name_list:
    file_path=os.path.join(path, f"{i.strip()}.h5ad")
    sams[i]=file_path
  return sams

sams=create_sams(args.name, args.path1)

sm=SAMAP(
  sams,
  f_maps=args.path2
  )
samap=sm.run()

today=datetime.datetime.now().strftime("%Y%m%d")
out_file_path=os.path.join(args.path3, f"{args.name.strip()}_{today}.pkl")
save_samap(sm,out_file_path)
