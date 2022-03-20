import sys
import torch
from torch.cuda.amp import autocast as autocast
import warnings
from model_fpn_bilstm import *
from config import *
from utils.utils import *
from utils.data import *
import copy
import time
torch.manual_seed(1)
np.random.seed(1)
device = torch.device(1)
torch.cuda.set_device(device)
conf_map = {'kd':kd_conf(),'kiba':kiba_conf(),'dude':dude_conf()}

warnings.filterwarnings("ignore")
dataset = sys.argv[1]
conf = conf_map[dataset]
split = get_split(conf)
(train_loader, valid_loader, test_loader) = get_data_loader(split, conf)

# criterion = torch.nn.BCELoss()
criterion = conf.criterion
bin_conf = BIN_config_DBPE()
drugVQA_conf = drugVQA_conf()

if conf.model_type == 'drugVQA_fpn':
    model = DrugVQA_fpn(drugVQA_conf, ResidualBlock, bin_conf).cuda()
elif conf.model_type == 'drugVQA':
    model = DrugVQA(modelArgs, ResidualBlock).cuda()
model = nn.DataParallel(model,device_ids=conf.device_ids)

print("dataset len: ",len(train_loader)*conf.batch_size)
max_auc = 0
lr = 7e-4
model_max = copy.deepcopy(model)

for i in range(conf.epochs):
    if i >= 20 and i % 10 == 0:
        lr = lr / 2
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    total_loss = 0
    y_pred = []
    y_label = []
    model = model.train()
    start = time.time()
    for batch_idx, (contactmap, input_text, properties, size, gene_emb) in enumerate(train_loader):

        hidden_state = init_hidden(conf.batch_size,conf.num_gpus)
        if conf.model_type == 'drugVQA':
            contactmap = torch.unsqueeze(contactmap[:,0,:,:],1)
        # with autocast():
        output = model(input_text.cuda(), contactmap.cuda().float(), hidden_state, gene_emb)
        if conf.model_type == 'drugVQA':
            output = output[0]
        loss = criterion(output.type(torch.DoubleTensor).squeeze().to(device),
                         torch.tensor(properties).type(torch.DoubleTensor).squeeze().to(device))

        y_label = y_label + torch.tensor(properties).flatten().tolist()
        y_pred = y_pred + output.flatten().tolist()

        total_loss += loss.data
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        if (batch_idx+1) % conf.print_step == 0:
            print('one step size: ', conf.batch_size*conf.print_step)
            print('one print step time: ', time.time() - start)
            print('Training at Epoch ' + str(i + 1) + ' iteration ' + str(batch_idx) + ' with loss ' + str(
                loss.cpu().detach().numpy()))
            start = time.time()
    print("train")
    print_metrics(y_label,y_pred)

    if (i+1) % conf.test_epo == 0:
        print("valid")
        with torch.set_grad_enabled(False):
            auc = test(valid_loader, model, conf)
            if auc > max_auc:
                model_max = copy.deepcopy(model)
                max_auc = auc

print("test")
test_model(test_loader, model_max, conf)
torch.save(model_max.module.state_dict(), conf.save_prefix +  'model_max.pth')



