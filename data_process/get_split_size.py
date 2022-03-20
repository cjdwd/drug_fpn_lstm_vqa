import sys
sys.path.append('../')

import os
from utils.utils import *
from utils.data import *

dataset = 'dude'

if dataset == 'kiba':
    root_dir = '/GPUFS/nsccgz_ywang_ywj/caojindong/TDC/'
    contactmap_path = root_dir + 'my_data/kiba/kiba_contactmap/'
    split = get_kiba_split()
    split_names = ['train','valid','test']
    uniprot_ids = []
    for split_name in split_names:
        print(split_name)
        print(len(split[split_name]))

elif dataset == 'kd':
    root_dir = '/GPUFS/nsccgz_ywang_ywj/caojindong/TDC/'
    contactmap_path = root_dir + 'my_data/bindingdb/Kd/Kd_contactmap/'
    split = get_kd_split()
    split_names = ['train', 'valid', 'test']
    uniprot_ids = []
    for split_name in split_names:
        print(split_name)
        print(len(split[split_name]))

elif dataset == 'dude':
    root_dir = '/GPUFS/nsccgz_ywang_ywj/caojindong/TDC/'
    contactmap_path = root_dir + 'my_data/dude/dude_contactmap/'
    split = get_dude_split()
    split_names = ['train', 'test']
    for split_name in split_names:
        print(split_name)
        print(len(split[split_name]))