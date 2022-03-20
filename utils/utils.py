import numpy as np
import re
import torch
import os
from torch.autograd import Variable
# from torch.utils.data import Dataset, DataLoader
smileLettersPath = '/GPUFS/nsccgz_ywang_ywj/caojindong/drugVQA/data/DUDE/voc/combinedVoc-wholeFour.voc'

def getLetters(path):
    with open(path, 'r') as f:
        chars = f.read().split()
    return chars
def create_variable(tensor):
    # Do cuda() before wrapping with variable
    if torch.cuda.is_available():
        return Variable(tensor)
    else:
        return Variable(tensor)
def replace_halogen(string):
    """Regex to replace Br and Cl with single letters"""
#     print(string)
    br = re.compile('Br')
    cl = re.compile('Cl')
    string = br.sub('R', string)
    string = cl.sub('L', string)
    return string
# Create necessary variables, lengths, and target
def make_variables(lines, properties,letters):
    sequence_and_length = [line2voc_arr(line,letters) for line in lines]
    vectorized_seqs = [sl[0] for sl in sequence_and_length]
    seq_lengths = torch.LongTensor([sl[1] for sl in sequence_and_length])
    return pad_sequences(vectorized_seqs, seq_lengths, properties)
def make_variables_per_item(line, p,letters):
    sequence_and_length = line2voc_arr(line,letters)
    vectorized_seq = sequence_and_length[0]
#     print(sequence_and_length)
    seq_length = torch.LongTensor([sequence_and_length[1]])
#     print(seq_length)
    return pad_sequence(vectorized_seq, seq_length, p)

def make_variables_seq(lines,letters):
    sequence_and_length = [line2voc_arr(line,letters) for line in lines]
    vectorized_seqs = [sl[0] for sl in sequence_and_length]
    seq_lengths = torch.LongTensor([sl[1] for sl in sequence_and_length])
    return pad_sequences_seq(vectorized_seqs, seq_lengths)
def line2voc_arr(line,letters):
    arr = []
    regex = '(\[[^\[\]]{1,10}\])'
    line = replace_halogen(line)
    char_list = re.split(regex, line)
    for li, char in enumerate(char_list):
        if char.startswith('['):
               arr.append(letterToIndex(char,letters)) 
        else:
            chars = [unit for unit in char]

            for i, unit in enumerate(chars):
                arr.append(letterToIndex(unit,letters))
    return arr, len(arr)
def letterToIndex(letter,smiles_letters):
    try:
        return smiles_letters.index(letter)
    except:
        print(letter)
        return 0
# pad sequences and sort the tensor

def pad_sequence(vectorized_seq, seq_len, p):
    seq_tensor = torch.zeros(77).long()
    # print(seq_len)    
    seq_len = seq_len[0]
    seq_tensor[0:seq_len] = torch.LongTensor(vectorized_seq)
    # Sort tensors by their length
    # seq_lengths, perm_idx = seq_lengths.sort(0, descending=True)
    # seq_tensor = seq_tensor[perm_idx]
    target = torch.LongTensor([p])
    # Also sort the target (countries) in the same order
    # Return variables
    # DataParallel requires everything to be a Variable
    return create_variable(seq_tensor),create_variable(seq_len),create_variable(target)

def pad_sequences(vectorized_seqs, seq_lengths, properties):
    seq_tensor = torch.zeros((len(vectorized_seqs), seq_lengths.max())).long()
    for idx, (seq, seq_len) in enumerate(zip(vectorized_seqs, seq_lengths)):
        seq_tensor[idx, :seq_len] = torch.LongTensor(seq)

    # Sort tensors by their length
    seq_lengths, perm_idx = seq_lengths.sort(0, descending=True)
    seq_tensor = seq_tensor[perm_idx]

    # Also sort the target (countries) in the same order
    target = properties.double()
    if len(properties):
        target = target[perm_idx]
    # Return variables
    # DataParallel requires everything to be a Variable
    return create_variable(seq_tensor),create_variable(seq_lengths),create_variable(target)
def pad_sequences_seq(vectorized_seqs, seq_lengths):
    seq_tensor = torch.zeros((len(vectorized_seqs), seq_lengths.max())).long()
    for idx, (seq, seq_len) in enumerate(zip(vectorized_seqs, seq_lengths)):
        seq_tensor[idx, :seq_len] = torch.LongTensor(seq)

    # Sort tensors by their length
    seq_lengths, perm_idx = seq_lengths.sort(0, descending=True)
#     print(seq_tensor)
    seq_tensor = seq_tensor[perm_idx]
    # Return variables
    # DataParallel requires everything to be a Variable
    return create_variable(seq_tensor), create_variable(seq_lengths)

def construct_vocabulary(smiles_list,fname):
    """Returns all the characters present in a SMILES file.
       Uses regex to find characters/tokens of the format '[x]'."""
    add_chars = set()
    for i, smiles in enumerate(smiles_list):
        regex = '(\[[^\[\]]{1,10}\])'
        smiles = ds.replace_halogen(smiles)
        char_list = re.split(regex, smiles)
        for char in char_list:
            if char.startswith('['):
                add_chars.add(char)
            else:
                chars = [unit for unit in char]
                [add_chars.add(unit) for unit in chars]

    print("Number of characters: {}".format(len(add_chars)))
    with open(fname, 'w') as f:
        f.write('<pad>' + "\n")
        for char in add_chars:
            f.write(char + "\n")
    return add_chars
def readLinesStrip(lines):
    for i in range(len(lines)):
        lines[i] = lines[i].rstrip('\n')
    return lines
def getProteinSeq(path,contactMapName):
    proteins = open(path+"/"+contactMapName).readlines()
    proteins = readLinesStrip(proteins)
    seq = proteins[1]
    return seq
def getProtein(path,contactMapName,contactMap = True):
    proteins = open(path+"/"+contactMapName.split('_')[0]).readlines()
    proteins = readLinesStrip(proteins)
    seq = proteins[1]
    if(contactMap):
        contactMap = []
        for i in range(2,len(proteins)):
            contactMap.append(proteins[i])
        return seq,contactMap
    else:
        return seq

smiles_letters_utils = getLetters(smileLettersPath)
def getTrainDataSet_origin(trainFoldPath):
    with open(trainFoldPath, 'r') as f:
        trainCpi_list = f.read().strip().split('\n')
    trainDataSet = [cpi.strip().split() for cpi in trainCpi_list]
    return trainDataSet#[[smiles, sequence, interaction],.....]

def getTrainDataSet(trainFoldPath):
    with open(trainFoldPath, 'r') as f:
        trainCpi_list = f.read().strip().split('\n')
    trainDataSet = []
    for i in trainCpi_list:
        length = line2voc_arr(i.strip().split()[0],smiles_letters_utils)[1]
        if length < 77:
            trainDataSet.append(i.strip().split())
#     trainDataSet = [cpi.strip().split() for cpi in trainCpi_list]
    return trainDataSet#[[smiles, sequence, interaction],.....]

def getTestProteinList(testFoldPath):
    testProteinList = readLinesStrip(open(testFoldPath).readlines())[0].split()
    return testProteinList#['kpcb_2i0eA_full','fabp4_2nnqA_full',....]
def getSeqContactDict_origin(contactPath,contactDictPath):# make a seq-contactMap dict 
    contactDict = open(contactDictPath).readlines()
    seqContactDict = {}
    for data in contactDict:
        _,contactMapName = data.strip().split(':')
        seq,contactMap = getProtein(contactPath,contactMapName)
        contactmap_np = [list(map(float, x.strip(' ').split(' '))) for x in contactMap]
        feature2D = np.expand_dims(contactmap_np, axis=0)
        feature2D = torch.FloatTensor(feature2D)    
        seqContactDict[seq] = feature2D
    return seqContactDict

def getSeqContactDict(contactPath,contactDictPath):# make a seq-contactMap dict
    contactDict = open(contactDictPath).readlines()
    seqContactDict = {}
    for data in contactDict:
        seq,contactMapName = data.strip().split(':')
        _,contactMap = getProtein(contactPath,contactMapName)
        contactmap_np = [list(map(float, x.strip(' ').split(' '))) for x in contactMap]
        feature2D = np.expand_dims(contactmap_np, axis=0)
        feature2D = torch.FloatTensor(feature2D)
        seqContactDict[seq] = feature2D
    return seqContactDict


def getDataDict(testProteinList,activePath,decoyPath,contactPath):
    dataDict = {}
    for x in testProteinList:#'xiap_2jk7A_full'
        xData = []
        protein = x.split('_')[0]
        print(protein)
        proteinActPath = activePath+"/"+protein+"_actives_final.ism"
        proteinDecPath = decoyPath+"/"+protein+"_decoys_final.ism"
        act = open(proteinActPath,'r').readlines()
        dec = open(proteinDecPath,'r').readlines()
        actives = [[x.split(' ')[0],1] for x in act] ######
        decoys = [[x.split(' ')[0],0] for x in dec]# test
        seq = getProtein(contactPath,x,contactMap = False)
        for i in range(len(actives)):
            xData.append([actives[i][0],seq,actives[i][1]])
        for i in range(len(decoys)):
            xData.append([decoys[i][0],seq,decoys[i][1]])
        print(len(xData))
        dataDict[x] = xData
    return dataDict
###BindingDB dataset utils
###
invalid_list = ["1kfxL","3oe6A","4zmaT","5w8lD","5is0E","4rewA","1f5fA","5ejvB"]
def getChemSmilesDict(chemPath,smilesPath):
    chem = open(chemPath).readlines()
    smiles = open(smilesPath).readlines()
    chemSmilesDict = {}
    for (x,y) in zip(chem,smiles):
        chemSmilesDict[x.strip()] = y.strip()
    return chemSmilesDict

def getProteinContactMapDict(distanceMapPath):
    proteinContactMapDict = {}
    distanceMapList = os.listdir(distanceMapPath)
#     print(distanceMapPath,distanceMapList)
    for mapName in distanceMapList:
        if mapName[0] == '.':
            continue
        _,distanceMap = getProtein(distanceMapPath,mapName)
#         for x in distanceMap:
#             print(x)
#             print(x.strip().split(' '))
#             time.sleep(1)
        distanceMap_np = [list(map(float, x.strip().split(' '))) for x in distanceMap]
        feature2D = np.expand_dims(distanceMap_np, axis=0)
        feature2D = torch.FloatTensor(feature2D)
        proteinContactMapDict[mapName] = feature2D
    return proteinContactMapDict

def getUniprotContactMapDict(uniprotPdbPath,distanceMapPath):
    data = open(uniprotPdbPath).readlines()
    uniprotPdbDict = {}
    proteinContactMapDict = getProteinContactMapDict(distanceMapPath)
    for line in data:
        uniprotId,pdbId = line.strip().split(",")
        if pdbId in proteinContactMapDict.keys():
            uniprotPdbDict[uniprotId] = proteinContactMapDict[pdbId]
    return uniprotPdbDict

smiles_letters_utils = getLetters(smileLettersPath)
def getTrainDataset_bindingDB(trainFoldPath,uniprotIdContactMapDict):
    chemPath = os.path.join(trainFoldPath,"chem")
    smilesPath = os.path.join(trainFoldPath,"chem.repr")
    chemSmilesDict = getChemSmilesDict(chemPath,smilesPath)
    pos = open(os.path.join(trainFoldPath,"edges.pos")).readlines()
    neg = open(os.path.join(trainFoldPath,"edges.neg")).readlines()
    dataset = []

    for line in pos:
        _,chem,_,uniprotId = line.strip().split(",")
        if uniprotId in uniprotIdContactMapDict.keys() and chem in chemSmilesDict.keys():
            length = line2voc_arr(chemSmilesDict[chem],smiles_letters_utils)[1]
            if length < 77:
                dataset.append((chemSmilesDict[chem],uniprotId,1))
    for line in neg:
        _,chem,_,uniprotId = line.strip().split(",")
        if uniprotId in uniprotIdContactMapDict.keys() and chem in chemSmilesDict.keys():
            length = line2voc_arr(chemSmilesDict[chem], smiles_letters_utils)[1]
            if length < 77:
                dataset.append((chemSmilesDict[chem], uniprotId, 0))
    return dataset

def getTrainDataset_bindingDB_dev(trainFoldPath,uniprotIdContactMapDict,devFoldPath):
    chemPath = os.path.join(devFoldPath,"chem")
    smilesPath = os.path.join(devFoldPath,"chem.repr")
    chemSmilesDict = getChemSmilesDict(chemPath,smilesPath)
    train_pos = open(os.path.join(trainFoldPath,"edges.pos")).readlines()
    train_neg = open(os.path.join(trainFoldPath,"edges.neg")).readlines()
    pos = open(os.path.join(devFoldPath,"edges.pos")).readlines()
    neg = open(os.path.join(devFoldPath,"edges.neg")).readlines()
    dataset = []
    
    train_uniprotId = []
    for line in train_pos:
        _,chem,_,uniprotId = line.strip().split(",")
        train_uniprotId.append(uniprotId)

    for line in train_neg:
        _,chem,_,uniprotId = line.strip().split(",")
        train_uniprotId.append(uniprotId)
                
    for line in pos:
        _,chem,_,uniprotId = line.strip().split(",")
        if uniprotId not in train_uniprotId:
            print("not in")
        if uniprotId in uniprotIdContactMapDict.keys() and chem in chemSmilesDict.keys() and uniprotId not in train_uniprotId:
            length = line2voc_arr(chemSmilesDict[chem],smiles_letters_utils)[1]
            if length < 77:
                dataset.append((chemSmilesDict[chem],uniprotId,1))
    for line in neg:
        _,chem,_,uniprotId = line.strip().split(",")
        if uniprotId in uniprotIdContactMapDict.keys() and chem in chemSmilesDict.keys() and uniprotId not in train_uniprotId:
            length = line2voc_arr(chemSmilesDict[chem], smiles_letters_utils)[1]
            if length < 77:
                dataset.append((chemSmilesDict[chem], uniprotId, 0))
    return dataset

###human
###
def getSeqDistanceMapDict_human(seqPdbPath,distanceMapPath):
    data = open(seqPdbPath).readlines()
    seqDistanceMapDict = {}
    proteinContactMapDict = getProteinContactMapDict(distanceMapPath)
    for line in data[1:]:
        seq,pdbId = line.strip().split(",")
        if pdbId in proteinContactMapDict.keys():
            seqDistanceMapDict[seq] = proteinContactMapDict[pdbId]
    return seqDistanceMapDict

def getTrainDataset_human(trainFoldPath,seqDistanceMapDict):
    data = open(trainFoldPath).readlines()
    dataset = []
    print(len(data))
    for line in data:
        chem = line.strip().split(" ")[0]
        length = line2voc_arr(chem,smiles_letters_utils)[1]
        if length < 77:
            seq = line.strip().split(" ")[1]
            p = line.strip().split(" ")[2]
            if seq in seqDistanceMapDict.keys():
                dataset.append((chem,seq,int(p)))
    return dataset

def getProteinSeqDict(contactDictPath):
    contactDict = open(contactDictPath).readlines()
    proteinSeqDict = {}
    for data in contactDict:
        seq, contactMapName = data.strip().split(':')
        protein = contactMapName.split("_")[0]
        proteinSeqDict[protein] = seq
    return proteinSeqDict