import sys
sys.path.append('../')

import torch
from sklearn import metrics
import warnings
from config import BIN_config_DBPE
from utils.stream import BIN_Data_Encoder, protein2emb_encoder, drug2emb_encoder, gene2emb_encoder
from torch.utils.data import Dataset, DataLoader
from tdc.multi_pred import DTI
import os
from utils.utils import *
from PIL import Image
from sklearn import metrics
import copy
from sklearn.metrics import roc_auc_score, average_precision_score, f1_score, roc_curve, confusion_matrix, \
    precision_score, recall_score
from sklearn.model_selection import KFold
from sklearn.metrics import auc as AUC
import pandas as pd

warnings.filterwarnings("ignore")
torch.cuda.set_device(1)
smileLettersPath = '/GPUFS/nsccgz_ywang_ywj/caojindong/drugVQA/data/DUDE/voc/combinedVoc-wholeFour.voc'
smiles_letters = getLetters(smileLettersPath)


def get_split(conf):
    if conf.read_processed_data:
        split = {}
        for split_name in conf.split_names:
            split[split_name] = pd.read_csv(conf.processed_data_path+conf.dataset+'_'+split_name+'.csv')

    else:
        if conf.dataset == 'kiba':
            split = get_kiba_split()

        elif conf.dataset == 'kd':
            split = get_kd_split()

        elif conf.dataset == 'dude':
            split = get_dude_split()

    return split


def get_kiba_split():
    data_kiba = DTI(name='KIBA')
    data_kiba.convert_to_log(form='binding')
    data_kiba.binarize(threshold=7.932, order='descending')
    split = data_kiba.get_split()

    root_dir = '/GPUFS/nsccgz_ywang_ywj/caojindong/TDC/'
    contactmap_path = root_dir + 'my_data/kiba/kiba_contactmap/'
    contactmap_files = os.listdir(contactmap_path)

    split_names = ['train', 'valid', 'test']
    uniprotID_sequence_map = {}

    for i in split_names:
        split[i] = split[i].dropna(axis=0, how='any')
        for index, row in split[i].iterrows():
            if row['Target_ID'] + '.npy' not in contactmap_files:
                split[i].drop(index, axis=0, inplace=True)
            if row['Target_ID'] not in uniprotID_sequence_map.keys():
                uniprotID_sequence_map[row['Target_ID']] = row['Target']
        split[i].reset_index(inplace=True)
        split[i].drop(['index'], axis=1, inplace=True)

    split['train'].rename(columns={'Drug': 'SMILES', 'Target': 'Target Sequence', 'Y': 'Label'}, inplace=True)
    split['valid'].rename(columns={'Drug': 'SMILES', 'Target': 'Target Sequence', 'Y': 'Label'}, inplace=True)
    split['test'].rename(columns={'Drug': 'SMILES', 'Target': 'Target Sequence', 'Y': 'Label'}, inplace=True)
    return split


def get_kd_split():
    data_Kd = DTI(name='BindingDB_Kd')
    data_Kd.convert_to_log(form='binding')
    data_Kd.binarize(threshold=5.0, order='descending')
    split = data_Kd.get_split()

    root_dir = '/GPUFS/nsccgz_ywang_ywj/caojindong/TDC/'
    contactmap_path = root_dir + 'my_data/bindingdb/Kd/Kd_contactmap/'
    contactmap_files = os.listdir(contactmap_path)

    split_names = ['train', 'valid', 'test']
    uniprotID_sequence_map = {}

    for i in split_names:
        split[i] = split[i].dropna(axis=0, how='any')
        for index, row in split[i].iterrows():
            if row['Target_ID'] + '.npy' not in contactmap_files:
                split[i].drop(index, axis=0, inplace=True)
            if row['Target_ID'] not in uniprotID_sequence_map.keys():
                uniprotID_sequence_map[row['Target_ID']] = row['Target']
        split[i].reset_index(inplace=True)
        split[i].drop(['index'], axis=1, inplace=True)

    split['train'].rename(columns={'Drug': 'SMILES', 'Target': 'Target Sequence', 'Y': 'Label'}, inplace=True)
    split['valid'].rename(columns={'Drug': 'SMILES', 'Target': 'Target Sequence', 'Y': 'Label'}, inplace=True)
    split['test'].rename(columns={'Drug': 'SMILES', 'Target': 'Target Sequence', 'Y': 'Label'}, inplace=True)
    return split


def get_dude_split():
    root_dir = '/GPUFS/nsccgz_ywang_ywj/caojindong/TDC/'
    contactmap_path = root_dir + 'my_data/dude/dude_contactmap/'
    contactmap_files = os.listdir(contactmap_path)

    split_names = ['train', 'test']
    df_train = pd.read_csv(root_dir + 'my_data/dude/dude_train.csv')
    df_test = pd.read_csv(root_dir + 'my_data/dude/dude_test.csv')
    split = {}
    split['train'] = df_train
    split['test'] = df_test
    for i in split_names:
        split[i] = split[i].dropna(axis=0, how='any')
        for index, row in split[i].iterrows():
            if row['Target ID'] + '.npy' not in contactmap_files:
                split[i].drop(index, axis=0, inplace=True)
        split[i].reset_index(inplace=True)
        split[i].drop(['index'], axis=1, inplace=True)

    split['train'].rename(columns={'Target ID': 'Target_ID'}, inplace=True)
    split['test'].rename(columns={'Target ID': 'Target_ID'}, inplace=True)
    return split


def getDataSet(data_csv):
    for index, row in data_csv.iterrows():
        length = line2voc_arr(row['SMILES'], smiles_letters)[1]
        if length > 77:
            data_csv.drop(index, axis=0, inplace=True)
    data_csv.reset_index(inplace=True)
    data_csv.drop(['index'], axis=1, inplace=True)
    return data_csv


class tdcDataset(Dataset):
    # Initialize your data, download, etc.
    def __init__(self, data_csv, conf, gene_csv):
        self.data_csv = data_csv
        self.len = len(data_csv)
        self.uniprotID_contactmap_dict = {}
        self.uniprotID_size_dict = {}
        self.contactmap_path = conf.contactmap_path
        self.uniprotID_gene_dict = {}
        self.conf = conf
        self.gene_csv = gene_csv
        self.Target_ID_gene_map = {}
        self.index_drug_map = {}
        for index,row in gene_csv.iterrows():
            self.Target_ID_gene_map[row['Target_ID']] = row['gene_seq']
        # print(gene_csv)

    def __getitem__(self, index):
        label = self.data_csv.iloc[index]['Label']
        smiles = self.data_csv.iloc[index]['SMILES']
        uniprotID = self.data_csv.iloc[index]['Target_ID']
        # p = self.data_csv.iloc[index]['Target Sequence']
        # d = self.data_csv.iloc[index]['SMILES']

        # print(self.gene_csv.loc[self.gene_csv['Target_ID'] == self.data_csv.iloc[index]['Target_ID']])
        # p_v, mask_p = protein2emb_encoder(p)
        # d_v, mask_d = drug2emb_encoder(d)
        gene_seq = self.Target_ID_gene_map[uniprotID]
        # print(uniprotID)
        # print(gene_seq)
        if uniprotID not in self.uniprotID_contactmap_dict.keys():
            contactmap = np.load(self.contactmap_path + uniprotID + '.npy')
            img = Image.fromarray(np.uint8(contactmap.squeeze() * 255))
            contactmap = np.asarray(img.resize((224, 224))) / 255
            size = contactmap.shape[1]
            contactmap = torch.from_numpy(np.expand_dims(contactmap, 0))
            contactmap = create_variable(contactmap)
            contactmap = torch.repeat_interleave(contactmap, repeats=3, dim=0)
            try:
                gene_emb = gene2emb_encoder(gene_seq,self.conf.gene_len)
            except:
                print(uniprotID)
            self.uniprotID_contactmap_dict[uniprotID] = contactmap
            self.uniprotID_size_dict[uniprotID] = torch.tensor(size)
            self.uniprotID_gene_dict[uniprotID] = torch.tensor(gene_emb).cuda().int()

        if index not in self.index_drug_map.keys():
            input_text, seq_lengths, y = make_variables_per_item(smiles, int(label), smiles_letters)
            self.index_drug_map[index] = input_text

        return self.uniprotID_contactmap_dict[uniprotID], self.index_drug_map[index], \
               torch.tensor(int(label)).cuda().float(), self.uniprotID_size_dict[uniprotID], self.uniprotID_gene_dict[uniprotID]

    def __len__(self):
        return self.len

def print_metrics(y_label,y_pred):
    #     print(y_label,y_pred)
    fpr, tpr, thresholds = roc_curve(y_label, y_pred)

    precision = tpr / (tpr + fpr)
    #     print(tpr,fpr,precision)
    f1 = 2 * precision * tpr / (tpr + precision + 0.00001)
    #     print(f1)
    thred_optim = thresholds[5:][np.argmax(f1[5:])]
    #     thred_optim = thresholds[np.argmax(f1)]
    print("optimal threshold: " + str(thred_optim))

    y_pred_s = [1 if i else 0 for i in (y_pred >= thred_optim)]

    auc_k = AUC(fpr, tpr)
    print("AUROC:" + str(auc_k))
    print("AUPRC: " + str(average_precision_score(y_label, y_pred)))

    cm1 = confusion_matrix(y_label, y_pred_s)
    print('Confusion Matrix : \n', cm1)
    print('Recall : ', recall_score(y_label, y_pred_s))
    print('Precision : ', precision_score(y_label, y_pred_s))

    total1 = sum(sum(cm1))
    #####from confusion matrix calculate accuracy
    accuracy1 = (cm1[0, 0] + cm1[1, 1]) / total1
    print('Accuracy : ', accuracy1)

    sensitivity1 = cm1[0, 0] / (cm1[0, 0] + cm1[0, 1])
    print('Sensitivity : ', sensitivity1)

    specificity1 = cm1[1, 1] / (cm1[1, 0] + cm1[1, 1])
    print('Specificity : ', specificity1)

    return auc_k

def init_hidden(batch_size, num_gpus):
    return (Variable(torch.zeros(int(4 * num_gpus), int(batch_size / num_gpus), 64).cuda()),
            Variable(torch.zeros(int(4 * num_gpus), int(batch_size / num_gpus), 64)).cuda())


def test(data_loader, model, conf):
    y_pred = []
    y_label = []
    model.eval()
    loss_accumulate = 0.0
    count = 0.0
    device = torch.device(0)
    criterion = conf.criterion

    for batch_idx, (contactmap, input_text, properties, size, gene_emb) in enumerate(data_loader):
        #         outputs = net(contactmap.cuda().float())
        #         print(outputs[0].shape,outputs[1].shape,outputs[2].shape,outputs[3].shape)
        #         final_conv = nn.Conv2d(256*4,256,7,bias=False).cuda()
        #         output = final_conv(torch.cat((outputs[0],outputs[1],outputs[2],outputs[3]),1))
        #         print(output.squeeze().shape)

        hidden_state = init_hidden(conf.batch_size, conf.num_gpus)
        if conf.model_type == 'drugVQA':
            contactmap = torch.unsqueeze(contactmap[:, 0, :, :], 1)
        output = model(input_text.cuda(), contactmap.cuda().float(), hidden_state, gene_emb)
        if conf.model_type == 'drugVQA':
            output = output[0]
        loss = criterion(output.type(torch.DoubleTensor).squeeze().to(device),
                         torch.tensor(properties).type(torch.DoubleTensor).squeeze().to(device))
        output = output.type(torch.DoubleTensor).squeeze()
        y_label = y_label + torch.tensor(properties).flatten().tolist()
        y_pred = y_pred + output.flatten().tolist()

        loss_accumulate += loss
        count += 1

        if (batch_idx % conf.print_step == 0):
            print('valid at iteration ' + str(batch_idx) + ' with loss ' + str(loss.cpu().detach().numpy()))

    return print_metrics(y_label,y_pred)


def test_model(data_loader, model_max, conf):
    with torch.set_grad_enabled(False):
        return test(data_loader, model_max, conf)


def get_data_loader(split, conf):
    ret = []
    gene_csv = pd.read_csv(conf.gene_csv_path)

    for i in conf.split_names:
        if conf.test_mode:
            sub_datacsv = split[i][0:1000]
        else:
            sub_datacsv = split[i]
        sub_dataset = tdcDataset(data_csv=sub_datacsv, conf=conf, gene_csv=gene_csv)
        sub_loader = DataLoader(dataset=sub_dataset, batch_size=conf.batch_size, shuffle=True, drop_last=True)
        ret.append(sub_loader)
    return ret

