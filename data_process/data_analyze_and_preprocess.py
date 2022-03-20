import sys
sys.path.append('../')

import os
import pandas as pd
from utils.utils import *
from utils.data import *

dataset = 'kd'

root_dir = '../data/gene_seq/'+dataset+'/'
files = os.listdir(root_dir)
df = pd.DataFrame(columns=['Target_ID','gene_seq'])
uniprotID_gene_map = {}
avg_seq_len = 0
long_count = 0
for file in files:
    gene_seqs = open(root_dir+file).readlines()
    seq = ""
    for gene_seq in gene_seqs[1:]:
        seq = seq + gene_seq.strip()
    df = df.append({'Target_ID':file.strip().split('.')[0],'gene_seq':seq},ignore_index=True)
    uniprotID_gene_map[file.strip().split('.')[0]] = seq
    avg_seq_len += len(seq)
    if len(seq) > 204800:
        long_count += 1

print(df)
df.to_csv('../data/gene_seq/'+dataset+'.csv',index=False)

if dataset == 'kiba':
    root_dir = '/GPUFS/nsccgz_ywang_ywj/caojindong/TDC/'
    contactmap_path = root_dir + 'my_data/kiba/kiba_contactmap/'
    split = get_kiba_split()
    split_names = ['train','valid','test']

elif dataset == 'kd':
    root_dir = '/GPUFS/nsccgz_ywang_ywj/caojindong/TDC/'
    contactmap_path = root_dir + 'my_data/bindingdb/Kd/Kd_contactmap/'
    split = get_kd_split()
    split_names = ['train', 'valid', 'test']

elif dataset == 'dude':
    root_dir = '/GPUFS/nsccgz_ywang_ywj/caojindong/TDC/'
    contactmap_path = root_dir + 'my_data/dude/dude_contactmap/'
    split = get_dude_split()
    split_names = ['train', 'test']

print(split['train'])
smileLettersPath = '/GPUFS/nsccgz_ywang_ywj/caojindong/drugVQA/data/DUDE/voc/combinedVoc-wholeFour.voc'
smiles_letters = getLetters(smileLettersPath)
for split_name in split_names:
    split[split_name]['gene_seq'] = ''
    for index, row in split[split_name].iterrows():
        length = line2voc_arr(row['SMILES'], smiles_letters)[1]
        if row['Target_ID'] not in uniprotID_gene_map.keys() or length > 77:
            split[split_name].drop(index, axis=0, inplace=True)
        # else:
        #     split[split_name].loc[split[split_name]['Target_ID']==row['Target_ID']]\
        #         = uniprotID_gene_map[row['Target_ID']]

    split[split_name].reset_index(inplace=True)
    split[split_name].drop(['index'], axis=1, inplace=True)
    split[split_name].to_csv('../data/processed_csv/'+dataset+'_'+split_name+'.csv')

print(split['train']['gene_seq'])
print("avg seq len: ", avg_seq_len/len(df))
print("long count: ", long_count)