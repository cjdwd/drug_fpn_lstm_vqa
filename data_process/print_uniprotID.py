import sys
sys.path.append('../')

import os
from utils.utils import *
from utils.data import *

dataset = 'kd'

if dataset == 'kiba':
    root_dir = '/GPUFS/nsccgz_ywang_ywj/caojindong/TDC/'
    contactmap_path = root_dir + 'my_data/kiba/kiba_contactmap/'
    split = get_kiba_split()
    split_names = ['train','valid','test']
    uniprot_ids = []
    for split_name in split_names:
        for i, row in split[split_name].iterrows():
            if row['Target_ID'] not in uniprot_ids:
                uniprot_ids.append(row['Target_ID'])

elif dataset == 'kd':
    root_dir = '/GPUFS/nsccgz_ywang_ywj/caojindong/TDC/'
    contactmap_path = root_dir + 'my_data/bindingdb/Kd/Kd_contactmap/'
    split = get_kd_split()
    split_names = ['train', 'valid', 'test']
    uniprot_ids = []
    for split_name in split_names:
        for i, row in split[split_name].iterrows():
            if row['Target_ID'] not in uniprot_ids:
                uniprot_ids.append(row['Target_ID'])

elif dataset == 'dude':
    root_dir = '/GPUFS/nsccgz_ywang_ywj/caojindong/TDC/'
    contactmap_path = root_dir + 'my_data/dude/dude_contactmap/'
    split = get_dude_split()
    split_names = ['train', 'test']
    for split_name in split_names:
        for i, row in split[split_name].iterrows():
            if row['Target_ID'] not in uniprot_ids:
                uniprot_ids.append(row['Target_ID'])