import torch

class kiba_conf():
    def __init__(self):
        self.dataset = 'kiba'
        self.batch_size = 192
        self.model_type = 'drugVQA_fpn'
        self.num_gpus = 4
        self.test_mode = False
        self.epochs = 60
        self.test_epo = 1
        self.print_step = 128
        self.split_names = ['train','valid','test']
        self.root_dir = '/GPUFS/nsccgz_ywang_ywj/caojindong/TDC/'
        self.contactmap_path = self.root_dir + 'my_data/kiba/kiba_contactmap/'
        self.save_prefix = './model_pkl/' + self.dataset + '/' + self.model_type + '_'

class kd_conf():
    def __init__(self):
        self.dataset = 'kd'

        ###train settings
        self.batch_size = 192
        self.model_type = 'drugVQA_fpn'
        self.num_gpus = 4
        self.device_ids = [1,0,2,3]
        self.test_mode = False
        self.epochs = 60
        self.test_epo = 1
        self.print_step = 1
        # self.criterion = torch.nn.functional.binary_cross_entropy_with_logits
        self.criterion = torch.nn.BCELoss()

        ###data settings
        self.gene_csv_path = './data/gene_seq/'+self.dataset+'.csv'
        self.gene_len = 204800
        self.read_processed_data = True
        self.split_names = ['train', 'valid', 'test']
        self.processed_data_path = './data/processed_csv/'
        self.root_dir = '/GPUFS/nsccgz_ywang_ywj/caojindong/TDC/'
        self.contactmap_path = self.root_dir + 'my_data/bindingdb/Kd/Kd_contactmap/'
        self.save_prefix = './model_pkl/' + self.dataset + '/' + self.model_type + '_'

class dude_conf():
    def __init__(self):
        self.dataset = 'dude'
        self.batch_size = 288
        self.model_type = 'drugVQA_fpn'
        self.num_gpus = 4
        self.test_mode = False
        self.epochs = 60
        self.test_epo = 10
        self.print_step = 128
        self.split_names = ['train', 'test', 'test']
        self.root_dir = '/GPUFS/nsccgz_ywang_ywj/caojindong/TDC/'
        self.contactmap_path = self.root_dir + 'my_data/dude/dude_contactmap/'
        self.save_prefix = './model_pkl/' + self.dataset + '/' + self.model_type + '_'


def BIN_config_DBPE():
    config = {}
    config['batch_size'] = 96
    config['input_dim_drug'] = 23532
    config['input_dim_target'] = 16693
    config['train_epoch'] = 13
    config['max_drug_seq'] = 50
    config['max_protein_seq'] = 545
    config['emb_size'] = 384
    config['dropout_rate'] = 0.1

    # DenseNet
    config['scale_down_ratio'] = 0.25
    config['growth_rate'] = 20
    config['transition_rate'] = 0.5
    config['num_dense_blocks'] = 4
    config['kernal_dense_size'] = 3

    # Encoder
    config['intermediate_size'] = 1536
    config['num_attention_heads'] = 12
    config['attention_probs_dropout_prob'] = 0.1
    config['hidden_dropout_prob'] = 0.1
    config['flat_dim'] = 78192
    return config

def drugVQA_conf():
    modelArgs = {}

    modelArgs['batch_size'] = 1
    modelArgs['lstm_hid_dim'] = 64
    modelArgs['d_a'] = 32
    modelArgs['r'] = 10
    modelArgs['n_chars_smi'] = 260
    modelArgs['n_chars_seq'] = 21
    modelArgs['dropout'] = 0.2
    modelArgs['in_channels'] = 8
    modelArgs['cnn_channels'] = 32
    modelArgs['cnn_layers'] = 4
    modelArgs['emb_dim'] = 30
    modelArgs['dense_hid'] = 64
    modelArgs['task_type'] = 0
    modelArgs['n_classes'] = 1
    modelArgs['gene_emb_dim'] = 4
    modelArgs['gene_hidden_dim'] = 2
    modelArgs['gene_kernel_size'] = 512
    modelArgs['gene_stride'] = 32
    modelArgs['gene_out_channels'] = 2
    modelArgs['gene_linear_in_dim'] = 6400*4
    modelArgs['gene_linear_out_dim'] = 128

    return modelArgs