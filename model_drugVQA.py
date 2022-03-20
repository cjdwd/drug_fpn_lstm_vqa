from __future__ import print_function
import torch.nn as nn
import torch
from torch.autograd import Variable
import torch.nn.functional as F

import torch.utils.data as Data
import numpy as np

import collections
import math
import copy
torch.manual_seed(1)
np.random.seed(1)

def conv3x3(in_channels, out_channels, stride=1):
    return nn.Conv2d(in_channels, out_channels, kernel_size=3,
                     stride=stride, padding=1, bias=False)
def conv5x5(in_channels, out_channels, stride=1):
    return nn.Conv2d(in_channels, out_channels, kernel_size=5,
                     stride=stride, padding=2, bias=False)
def conv1x1(in_channels, out_channels, stride=1):
    return nn.Conv2d(in_channels, out_channels, kernel_size=1,
                     stride=stride, padding=0, bias=False)
# Residual block
class ResidualBlock(nn.Module):
    def __init__(self, in_channels, out_channels, stride=1, downsample=None):
        super(ResidualBlock, self).__init__()
        self.conv1 = conv5x5(in_channels, out_channels, stride)
        self.bn1 = nn.BatchNorm2d(out_channels)
        self.elu = nn.ELU(inplace=True)
        self.conv2 = conv3x3(out_channels, out_channels)
        self.bn2 = nn.BatchNorm2d(out_channels)
        self.downsample = downsample

    def forward(self, x):
        residual = x
        out = self.conv1(x)
        out = self.bn1(out)
        out = self.elu(out)
        out = self.conv2(out)
        out = self.bn2(out)
        if self.downsample:
            residual = self.downsample(x)
        out += residual
        out = self.elu(out)
        return out

class DrugVQA(torch.nn.Module):
    """
    The class is an implementation of the DrugVQA model including regularization and without pruning. 
    Slight modifications have been done for speedup
    
    """
    def __init__(self,args,block):
        """
        Initializes parameters suggested in paper
 
        args:
            batch_size  : {int} batch_size used for training
            lstm_hid_dim: {int} hidden dimension for lstm
            d_a         : {int} hidden dimension for the dense layer
            r           : {int} attention-hops or attention heads
            n_chars_smi : {int} voc size of smiles
            n_chars_seq : {int} voc size of protein sequence
            dropout     : {float}
            in_channels : {int} channels of CNN block input
            cnn_channels: {int} channels of CNN block
            cnn_layers  : {int} num of layers of each CNN block
            emb_dim     : {int} embeddings dimension
            dense_hid   : {int} hidden dim for the output dense
            task_type   : [0,1] 0-->binary_classification 1-->multiclass classification
            n_classes   : {int} number of classes
 
        Returns:
            self
        """
        super(DrugVQA,self).__init__()
        self.batch_size = args['batch_size']
        self.lstm_hid_dim = args['lstm_hid_dim']
        self.r = args['r']
        self.type = args['task_type']
        self.in_channels = args['in_channels']
        #rnn
        self.embeddings = nn.Embedding(args['n_chars_smi'], args['emb_dim'])
        self.seq_embed = nn.Embedding(args['n_chars_seq'],args['emb_dim'])
        self.lstm = torch.nn.LSTM(args['emb_dim'],self.lstm_hid_dim,2,batch_first=True,bidirectional=True,dropout=args['dropout']) 
        self.linear_first = torch.nn.Linear(2*self.lstm_hid_dim,args['d_a'])
        self.linear_second = torch.nn.Linear(args['d_a'],args['r'])
        self.linear_first_seq = torch.nn.Linear(args['cnn_channels'],args['d_a'])
        self.linear_second_seq = torch.nn.Linear(args['d_a'],self.r)

        #cnn
        self.conv = conv3x3(1, self.in_channels)
        self.bn = nn.BatchNorm2d(self.in_channels)
        self.elu = nn.ELU(inplace=False)
        self.layer1 = self.make_layer(block, args['cnn_channels'], args['cnn_layers'])
        self.layer2 = self.make_layer(block, args['cnn_channels'], args['cnn_layers'])

        self.linear_final_step = torch.nn.Linear(self.lstm_hid_dim*2+args['d_a'],args['dense_hid'])
        self.linear_final = torch.nn.Linear(args['dense_hid'],args['n_classes'])

#         self.hidden_state = self.init_hidden()
#         self.seq_hidden_state = self.init_hidden()
        
    def softmax(self,input, axis=1):
        """
        Softmax applied to axis=n
        Args:
           input: {Tensor,Variable} input on which softmax is to be applied
           axis : {int} axis on which softmax is to be applied

        Returns:
            softmaxed tensors
        """
        input_size = input.size()
        trans_input = input.transpose(axis, len(input_size)-1)
        trans_size = trans_input.size()
        input_2d = trans_input.contiguous().view(-1, trans_size[-1])
        soft_max_2d = F.softmax(input_2d)
        soft_max_nd = soft_max_2d.view(*trans_size)
        return soft_max_nd.transpose(axis, len(input_size)-1)

    def init_hidden(self,batch_size):
        return (Variable(torch.zeros(4,batch_size,self.lstm_hid_dim).cuda()),Variable(torch.zeros(4,batch_size,self.lstm_hid_dim)).cuda())
    
    def make_layer(self, block, out_channels, blocks, stride=1):
        downsample = None
        if (stride != 1) or (self.in_channels != out_channels):
            downsample = nn.Sequential(
                conv3x3(self.in_channels, out_channels, stride=stride),
                nn.BatchNorm2d(out_channels))
        layers = []
        layers.append(block(self.in_channels, out_channels, stride, downsample))
        self.in_channels = out_channels
        for i in range(1, blocks):
            layers.append(block(out_channels, out_channels))
        return nn.Sequential(*layers)
        
    # x1 = smiles , x2 = contactMap
#     def forward(self,x1,x2,size,hidden_state):
#         smile_embed = self.embeddings(x1)   
# #         print(x1.shape,smile_embed.shape,hidden_state[0].shape)
#         outputs, hidden_state = self.lstm(smile_embed,hidden_state)    
#         sentence_att = F.tanh(self.linear_first(outputs))       
#         sentence_att = self.linear_second(sentence_att)       
#         sentence_att = self.softmax(sentence_att,1)       
#         sentence_att = sentence_att.transpose(1,2)        
#         sentence_embed = sentence_att@outputs
#         avg_sentence_embed = torch.sum(sentence_embed,1)/self.r  #multi head
         
#         holder = torch.empty((0,32)).cuda()
#         for i in range(size.shape[0]):
#             x_temp = x2[i,0:size[i],0:size[i]].unsqueeze(0)
#             pic = self.conv(x_temp)
#             pic = self.bn(pic)
#             pic = self.elu(pic)
#             pic = self.layer1(pic)
#             pic = self.layer2(pic)
#             pic_emb = torch.mean(pic,2)
#             pic_emb = pic_emb.permute(0,2,1)
#             seq_att = F.tanh(self.linear_first_seq(pic_emb))       
#             seq_att = self.linear_second_seq(seq_att)       
#             seq_att = self.softmax(seq_att,1)       
#             seq_att = seq_att.transpose(1,2)
#             seq_embed = seq_att@pic_emb  
#             avg_seq_embed = torch.sum(seq_embed,1)/self.r
#             holder = torch.cat([holder.to(avg_seq_embed.dtype),avg_seq_embed],dim=0)
# #         pic = self.conv(x2)
# #         pic = self.bn(pic)
# #         pic = self.elu(pic)
# #         pic = self.layer1(pic)
# #         pic = self.layer2(pic)
# #         pic_emb = torch.mean(pic,2)
# #         pic_emb = pic_emb.permute(0,2,1)
# #         seq_att = F.tanh(self.linear_first_seq(pic_emb))       
# #         seq_att = self.linear_second_seq(seq_att)       
# #         seq_att = self.softmax(seq_att,1)       
# #         seq_att = seq_att.transpose(1,2)
# #         seq_embed = seq_att@pic_emb      
# #         avg_seq_embed = torch.sum(seq_embed,1)/self.r
        
# #         sscomplex = torch.cat([avg_sentence_embed,avg_seq_embed],dim=1) 
#         sscomplex = torch.cat([avg_sentence_embed,holder],dim=1) 
#         sscomplex = F.relu(self.linear_final_step(sscomplex))
        
#         if not bool(self.type):
#             output = F.sigmoid(self.linear_final(sscomplex))
#             return output,seq_att
#         else:
#             return F.log_softmax(self.linear_final(sscomplex)),seq_att
    #x1 smiles, x2 contactmap    
    def forward(self,x1,x2,hidden_state):
        #########smiles
#         print('smiles: ',x1.shape)
        smile_embed = self.embeddings(x1)    
#         print('smile_embed: ',smile_embed.shape)
        outputs, hidden_state = self.lstm(smile_embed,hidden_state)    
        sentence_att = F.tanh(self.linear_first(outputs))       
        sentence_att = self.linear_second(sentence_att)       
        sentence_att = self.softmax(sentence_att,1)       
        sentence_att = sentence_att.transpose(1,2)        
        sentence_embed = sentence_att@outputs
        avg_sentence_embed = torch.sum(sentence_embed,1)/self.r  #multi head
#         print('avg_sentence_embed: ',avg_sentence_embed.shape)
        #########contactmap
        pic = self.conv(x2)
        pic = self.bn(pic)
        pic = self.elu(pic)
        pic = self.layer1(pic)
        pic = self.layer2(pic)
        pic_emb = torch.mean(pic,2)
        pic_emb = pic_emb.permute(0,2,1)
        seq_att = F.tanh(self.linear_first_seq(pic_emb))       
        seq_att = self.linear_second_seq(seq_att)       
        seq_att = self.softmax(seq_att,1)       
        seq_att = seq_att.transpose(1,2)
        seq_embed = seq_att@pic_emb      
        avg_seq_embed = torch.sum(seq_embed,1)/self.r
        
        sscomplex = torch.cat([avg_sentence_embed,avg_seq_embed],dim=1) 
        sscomplex = F.relu(self.linear_final_step(sscomplex))
        
        if not bool(self.type):
            output = F.sigmoid(self.linear_final(sscomplex))
            return output,seq_att
        else:
            return F.log_softmax(self.linear_final(sscomplex)),seq_att

class DrugVQA_3modal(torch.nn.Module):
    """
    The class is an implementation of the DrugVQA model including regularization and without pruning. 
    Slight modifications have been done for speedup
    
    """
    def __init__(self,args,block,config):
        """
        Initializes parameters suggested in paper
 
        args:
            batch_size  : {int} batch_size used for training
            lstm_hid_dim: {int} hidden dimension for lstm
            d_a         : {int} hidden dimension for the dense layer
            r           : {int} attention-hops or attention heads
            n_chars_smi : {int} voc size of smiles
            n_chars_seq : {int} voc size of protein sequence
            dropout     : {float}
            in_channels : {int} channels of CNN block input
            cnn_channels: {int} channels of CNN block
            cnn_layers  : {int} num of layers of each CNN block
            emb_dim     : {int} embeddings dimension
            dense_hid   : {int} hidden dim for the output dense
            task_type   : [0,1] 0-->binary_classification 1-->multiclass classification
            n_classes   : {int} number of classes
 
        Returns:
            self
        """
        super(DrugVQA_3modal,self).__init__()
        ###MolTrans
        self.max_p = config['max_protein_seq']
        self.input_dim_target = config['input_dim_target']
        self.emb_size = config['emb_size']
        self.dropout_rate = config['dropout_rate']
        self.pemb = Embeddings(self.input_dim_target, self.emb_size, self.max_p, self.dropout_rate)
        self.n_layer = 2
        self.hidden_size = config['emb_size']
        self.intermediate_size = config['intermediate_size']
        self.num_attention_heads = config['num_attention_heads']
        self.attention_probs_dropout_prob = config['attention_probs_dropout_prob']
        self.hidden_dropout_prob = config['hidden_dropout_prob']
        self.p_encoder = Encoder_MultipleLayers(self.n_layer, self.hidden_size, self.intermediate_size, self.num_attention_heads, self.attention_probs_dropout_prob, self.hidden_dropout_prob)
        
        ###drugVQA
        self.batch_size = args['batch_size']
        self.lstm_hid_dim = args['lstm_hid_dim']
        self.r = args['r']
        self.type = args['task_type']
        self.in_channels = args['in_channels']
        #rnn
        self.embeddings = nn.Embedding(args['n_chars_smi'], args['emb_dim'])
        self.seq_embed = nn.Embedding(args['n_chars_seq'],args['emb_dim'])
        self.lstm = torch.nn.LSTM(args['emb_dim'],self.lstm_hid_dim,2,batch_first=True,bidirectional=True,dropout=args['dropout']) 
        self.linear_first = torch.nn.Linear(2*self.lstm_hid_dim,args['d_a'])
        self.linear_second = torch.nn.Linear(args['d_a'],args['r'])
        self.linear_first_seq = torch.nn.Linear(args['cnn_channels'],args['d_a'])
        self.linear_second_seq = torch.nn.Linear(args['d_a'],self.r)

        #cnn
        self.conv = conv3x3(1, self.in_channels)
        self.bn = nn.BatchNorm2d(self.in_channels)
        self.elu = nn.ELU(inplace=False)
        self.layer1 = self.make_layer(block, args['cnn_channels'], args['cnn_layers'])
        self.layer2 = self.make_layer(block, args['cnn_channels'], args['cnn_layers'])

        self.linear_final_step = torch.nn.Linear(self.lstm_hid_dim*2+256,args['dense_hid'])
        self.linear_final = torch.nn.Linear(args['dense_hid'],args['n_classes'])
 
        self.bn1 = nn.BatchNorm1d(545)
        self.bn2 = nn.BatchNorm2d(1)
        self.bn3 = nn.BatchNorm1d(256)

        self.ln1 = torch.nn.Linear(384,108)
        self.ln2 = torch.nn.Linear(108,32)
        self.fpn=FPN([3,4,6,3])
        self.sigmoid = nn.Sigmoid()
        self.final_conv = nn.Conv2d(256*4,256,7,bias=False).cuda()
#         self.hidden_state = self.init_hidden()
#         self.seq_hidden_state = self.init_hidden()
#         self.pic_attention_encoder = Encoder(64,64,2,0.1,0.1)
#         self.pic_linear_after_attention = nn.Linear(64, 1)
#         self.mutilayer_sentence_attention_encoder = Encoder_MultipleLayers(6,30,30,1,0.1,0.1)
#         self.sentence_linear_after_attention = nn.Linear(30, 128)
    def softmax(self,input, axis=1):
        """
        Softmax applied to axis=n
        Args:
           input: {Tensor,Variable} input on which softmax is to be applied
           axis : {int} axis on which softmax is to be applied

        Returns:
            softmaxed tensors
        """
        input_size = input.size()
        trans_input = input.transpose(axis, len(input_size)-1)
        trans_size = trans_input.size()
        input_2d = trans_input.contiguous().view(-1, trans_size[-1])
        soft_max_2d = F.softmax(input_2d)
        soft_max_nd = soft_max_2d.view(*trans_size)
        return soft_max_nd.transpose(axis, len(input_size)-1)

    def init_hidden(self,batch_size):
        return (Variable(torch.zeros(4,batch_size,self.lstm_hid_dim).cuda()),Variable(torch.zeros(4,batch_size,self.lstm_hid_dim)).cuda())
    
    def make_layer(self, block, out_channels, blocks, stride=1):
        downsample = None
        if (stride != 1) or (self.in_channels != out_channels):
            downsample = nn.Sequential(
                conv3x3(self.in_channels, out_channels, stride=stride),
                nn.BatchNorm2d(out_channels))
        layers = []
        layers.append(block(self.in_channels, out_channels, stride, downsample))
        self.in_channels = out_channels
        for i in range(1, blocks):
            layers.append(block(out_channels, out_channels))
        return nn.Sequential(*layers)
        
#     def forward(self,x1,x2,hidden_state,d,p,d_mask,p_mask):
#         smile_embed = self.embeddings(x1)         
#         outputs, hidden_state = self.lstm(smile_embed,hidden_state)    
#         sentence_att = F.tanh(self.linear_first(outputs))       
#         sentence_att = self.linear_second(sentence_att)       
#         sentence_att = self.softmax(sentence_att,1)       
#         sentence_att = sentence_att.transpose(1,2)        
#         sentence_embed = sentence_att@outputs
#         avg_sentence_embed = torch.sum(sentence_embed,1)/self.r  #multi head
        
        
#         outputs = self.fpn(x2)
#         #output shape(n,256)
#         output = self.final_conv(torch.cat((outputs[0],outputs[1],outputs[2],outputs[3]),1)).squeeze()
        
#         output = self.bn3(output)
#         output = self.sigmoid(output)
#         #print(avg_sentence_embed.shape,avg_seq_embed.shape)
#         sscomplex = torch.cat([avg_sentence_embed,output],dim=1) 
#         sscomplex = F.relu(self.linear_final_step(sscomplex))
        
#         if not bool(self.type):
#             output = F.sigmoid(self.linear_final(sscomplex))
#             return output
#         else:
#             return F.log_softmax(self.linear_final(sscomplex))
        
    def get_protein_feature(self,x1,x2,hidden_state,d,p,d_mask,p_mask):
        outputs = self.fpn(x2)
        return outputs
    
    def get_drug_feature(self,x1,x2,hidden_state,d,p,d_mask,p_mask):
        smile_embed = self.embeddings(x1)    
#         print('smile_embed: ',smile_embed.shape)
        outputs, hidden_state = self.lstm(smile_embed,hidden_state)    
        sentence_att = F.tanh(self.linear_first(outputs))       
        sentence_att = self.linear_second(sentence_att)       
        sentence_att = self.softmax(sentence_att,1)       
        sentence_att = sentence_att.transpose(1,2)        
#         sentence_embed = sentence_att@outputs
#         avg_sentence_embed = torch.sum(sentence_embed,1)/self.r
        return sentence_att
    
    def forward(self,x1,x2,hidden_state,d,p,d_mask,p_mask):
        pool2 = nn.MaxPool2d(kernel_size=8,stride=8)
        pool3 = nn.MaxPool2d(kernel_size=4,stride=4)
        pool4 = nn.MaxPool2d(kernel_size=2,stride=2)
        
        smile_embed = self.embeddings(x1)    
#         print('smile_embed: ',smile_embed.shape)
        outputs, hidden_state = self.lstm(smile_embed,hidden_state)    
        sentence_att = F.tanh(self.linear_first(outputs))       
        sentence_att = self.linear_second(sentence_att)       
        sentence_att = self.softmax(sentence_att,1)       
        sentence_att = sentence_att.transpose(1,2)        
        sentence_embed = sentence_att@outputs
        avg_sentence_embed = torch.sum(sentence_embed,1)/self.r
        
        
        outputs = self.fpn(x2)
        #output shape(n,256)
        output = self.final_conv(torch.cat((pool2(outputs[0]),pool3(outputs[1]),pool4(outputs[2]),outputs[3]),1)).squeeze()
#         scores = output.clone()
#         #scores shape(n,4,64)
#         scores = scores.view((output.shape[0],4,-1))
#         scores = self.pic_attention_encoder(scores)
#         #scores shape(n,4,1)
#         scores = self.pic_linear_after_attention(scores)
#         #output shape(n,4,64)
#         output = output.view((output.shape[0],4,-1))
#         #output shape(n,4,64)
#         output = output*scores
#         #output shape(n,256)
#         output = output.view((output.shape[0],-1))
        
        
        output = self.bn3(output)
        output = self.sigmoid(output)
        #print(avg_sentence_embed.shape,avg_seq_embed.shape)
        sscomplex = torch.cat([avg_sentence_embed,output],dim=1) 
        sscomplex = F.relu(self.linear_final_step(sscomplex))
        
        if not bool(self.type):
            output = F.sigmoid(self.linear_final(sscomplex))
            return output
        else:
            return F.log_softmax(self.linear_final(sscomplex))
        
class SelfOutput(nn.Module):
    def __init__(self, hidden_size, hidden_dropout_prob):
        super(SelfOutput, self).__init__()
        self.dense = nn.Linear(hidden_size, hidden_size)
        self.LayerNorm = LayerNorm(hidden_size)
        self.dropout = nn.Dropout(hidden_dropout_prob)

    def forward(self, hidden_states, input_tensor):
        hidden_states = self.dense(hidden_states)
        hidden_states = self.dropout(hidden_states)
        hidden_states = self.LayerNorm(hidden_states + input_tensor)
        return hidden_states    

class SelfAttention(nn.Module):
    def __init__(self, hidden_size, num_attention_heads, attention_probs_dropout_prob):
        super(SelfAttention, self).__init__()
        if hidden_size % num_attention_heads != 0:
            raise ValueError(
                "The hidden size (%d) is not a multiple of the number of attention "
                "heads (%d)" % (hidden_size, num_attention_heads))
        self.num_attention_heads = num_attention_heads
        self.attention_head_size = int(hidden_size / num_attention_heads)
        self.all_head_size = self.num_attention_heads * self.attention_head_size

        self.query = nn.Linear(hidden_size, self.all_head_size)
        self.key = nn.Linear(hidden_size, self.all_head_size)
        self.value = nn.Linear(hidden_size, self.all_head_size)

        self.dropout = nn.Dropout(attention_probs_dropout_prob)

    def transpose_for_scores(self, x):
        new_x_shape = x.size()[:-1] + (self.num_attention_heads, self.attention_head_size)
        x = x.view(*new_x_shape)
        return x.permute(0, 2, 1, 3)

    def forward(self, hidden_states):
        mixed_query_layer = self.query(hidden_states)
        mixed_key_layer = self.key(hidden_states)
        mixed_value_layer = self.value(hidden_states)

        query_layer = self.transpose_for_scores(mixed_query_layer)
        key_layer = self.transpose_for_scores(mixed_key_layer)
        value_layer = self.transpose_for_scores(mixed_value_layer)

        # Take the dot product between "query" and "key" to get the raw attention scores.
        attention_scores = torch.matmul(query_layer, key_layer.transpose(-1, -2))
        attention_scores = attention_scores / math.sqrt(self.attention_head_size)

        attention_scores = attention_scores

        # Normalize the attention scores to probabilities.
        attention_probs = nn.Softmax(dim=-1)(attention_scores)

        # This is actually dropping out entire tokens to attend to, which might
        # seem a bit unusual, but is taken from the original Transformer paper.
        attention_probs = self.dropout(attention_probs)

        context_layer = torch.matmul(attention_probs, value_layer)
        context_layer = context_layer.permute(0, 2, 1, 3).contiguous()
        new_context_layer_shape = context_layer.size()[:-2] + (self.all_head_size,)
        context_layer = context_layer.view(*new_context_layer_shape)
        return context_layer    

class LayerNorm(nn.Module):
    def __init__(self, hidden_size, variance_epsilon=1e-12):

        super(LayerNorm, self).__init__()
        self.gamma = nn.Parameter(torch.ones(hidden_size))
        self.beta = nn.Parameter(torch.zeros(hidden_size))
        self.variance_epsilon = variance_epsilon

    def forward(self, x):
        u = x.mean(-1, keepdim=True)
        s = (x - u).pow(2).mean(-1, keepdim=True)
        x = (x - u) / torch.sqrt(s + self.variance_epsilon)
        return self.gamma * x + self.beta
    
class Attention(nn.Module):
    def __init__(self, hidden_size, num_attention_heads, attention_probs_dropout_prob, hidden_dropout_prob):
        super(Attention, self).__init__()
        self.self = SelfAttention(hidden_size, num_attention_heads, attention_probs_dropout_prob)
        self.output = SelfOutput(hidden_size, hidden_dropout_prob)

    def forward(self, input_tensor):
        self_output = self.self(input_tensor)
        attention_output = self.output(self_output, input_tensor)
        return attention_output    
    
class Intermediate(nn.Module):
    def __init__(self, hidden_size, intermediate_size):
        super(Intermediate, self).__init__()
        self.dense = nn.Linear(hidden_size, intermediate_size)

    def forward(self, hidden_states):
        hidden_states = self.dense(hidden_states)
        hidden_states = F.relu(hidden_states)
        return hidden_states

class Output(nn.Module):
    def __init__(self, intermediate_size, hidden_size, hidden_dropout_prob):
        super(Output, self).__init__()
        self.dense = nn.Linear(intermediate_size, hidden_size)
        self.LayerNorm = LayerNorm(hidden_size)
        self.dropout = nn.Dropout(hidden_dropout_prob)

    def forward(self, hidden_states, input_tensor):
        hidden_states = self.dense(hidden_states)
        hidden_states = self.dropout(hidden_states)
        hidden_states = self.LayerNorm(hidden_states + input_tensor)
        return hidden_states

class Encoder(nn.Module):
    def __init__(self, hidden_size, intermediate_size, num_attention_heads, attention_probs_dropout_prob, hidden_dropout_prob):
        super(Encoder, self).__init__()
        self.attention = Attention(hidden_size, num_attention_heads, attention_probs_dropout_prob, hidden_dropout_prob)
        self.intermediate = Intermediate(hidden_size, intermediate_size)
        self.output = Output(intermediate_size, hidden_size, hidden_dropout_prob)

    def forward(self, hidden_states):
        attention_output = self.attention(hidden_states)
        intermediate_output = self.intermediate(attention_output)
        layer_output = self.output(intermediate_output, attention_output)
        return layer_output    

    
class Encoder_MultipleLayers(nn.Module):
    def __init__(self, n_layer, hidden_size, intermediate_size, num_attention_heads, attention_probs_dropout_prob, hidden_dropout_prob):
        super(Encoder_MultipleLayers, self).__init__()
        layer = Encoder(hidden_size, intermediate_size, num_attention_heads, attention_probs_dropout_prob, hidden_dropout_prob)
        self.layer = nn.ModuleList([copy.deepcopy(layer) for _ in range(n_layer)])    

    def forward(self, hidden_states, output_all_encoded_layers=True):
        all_encoder_layers = []
        for layer_module in self.layer:
            hidden_states = layer_module(hidden_states)
            #if output_all_encoded_layers:
            #    all_encoder_layers.append(hidden_states)
        #if not output_all_encoded_layers:
        #    all_encoder_layers.append(hidden_states)
        return hidden_states

class LayerNorm(nn.Module):
    def __init__(self, hidden_size, variance_epsilon=1e-12):

        super(LayerNorm, self).__init__()
        self.gamma = nn.Parameter(torch.ones(hidden_size))
        self.beta = nn.Parameter(torch.zeros(hidden_size))
        self.variance_epsilon = variance_epsilon

    def forward(self, x):
        u = x.mean(-1, keepdim=True)
        s = (x - u).pow(2).mean(-1, keepdim=True)
        x = (x - u) / torch.sqrt(s + self.variance_epsilon)
        return self.gamma * x + self.beta


class Embeddings(nn.Module):
    """Construct the embeddings from protein/target, position embeddings.
    """
    def __init__(self, vocab_size, hidden_size, max_position_size, dropout_rate):
        super(Embeddings, self).__init__()
        self.word_embeddings = nn.Embedding(vocab_size, hidden_size)
        self.position_embeddings = nn.Embedding(max_position_size, hidden_size)

        self.LayerNorm = LayerNorm(hidden_size)
        self.dropout = nn.Dropout(dropout_rate)

    def forward(self, input_ids):
        seq_length = input_ids.size(1)
        position_ids = torch.arange(seq_length, dtype=torch.long, device=input_ids.device)
        position_ids = position_ids.unsqueeze(0).expand_as(input_ids)
        
        words_embeddings = self.word_embeddings(input_ids)
        position_embeddings = self.position_embeddings(position_ids)

        embeddings = words_embeddings + position_embeddings
        embeddings = self.LayerNorm(embeddings)
        embeddings = self.dropout(embeddings)
        return embeddings
    
import torch.nn as nn
import torch.nn.functional as F
import math
class FeatureSelectionModule(nn.Module):
    def __init__(self, in_chan, out_chan, norm="GN"):
        super(FeatureSelectionModule, self).__init__()
        self.conv_atten = Conv2d(in_chan, in_chan, kernel_size=1, bias=False, norm=get_norm(norm, in_chan))
        self.sigmoid = nn.Sigmoid()
        self.conv = Conv2d(in_chan, out_chan, kernel_size=1, bias=False, norm=get_norm('', out_chan))
        weight_init.c2_xavier_fill(self.conv_atten)
        weight_init.c2_xavier_fill(self.conv)

    def forward(self, x):
        atten = self.sigmoid(self.conv_atten(F.avg_pool2d(x, x.size()[2:])))
        feat = torch.mul(x, atten)
        x = x + feat
        feat = self.conv(x)
        return feat


class FeatureAlign_V2(nn.Module):  # FaPN full version
    def __init__(self, in_nc=128, out_nc=128, norm=None):
        super(FeatureAlign_V2, self).__init__()
        self.lateral_conv = FeatureSelectionModule(in_nc, out_nc, norm="")
        self.offset = Conv2d(out_nc * 2, out_nc, kernel_size=1, stride=1, padding=0, bias=False, norm=norm)
        self.dcpack_L2 = dcn_v2(out_nc, out_nc, 3, stride=1, padding=1, dilation=1, deformable_groups=8,
                                extra_offset_mask=True)
        self.relu = nn.ReLU(inplace=True)
        weight_init.c2_xavier_fill(self.offset)

    def forward(self, feat_l, feat_s, main_path=None):
        HW = feat_l.size()[2:]
        if feat_l.size()[2:] != feat_s.size()[2:]:
            feat_up = F.interpolate(feat_s, HW, mode='bilinear', align_corners=False)
        else:
            feat_up = feat_s
        feat_arm = self.lateral_conv(feat_l)  # 0~1 * feats
        offset = self.offset(torch.cat([feat_arm, feat_up * 2], dim=1))  # concat for offset by compute the dif
        feat_align = self.relu(self.dcpack_L2([feat_up, offset], main_path))  # [feat, offset]
        return feat_align + feat_arm
 
#ResNet的基本Bottleneck类
class Bottleneck(nn.Module):
    expansion=4#通道倍增数
    def __init__(self,in_planes,planes,stride=1,downsample=None):
        super(Bottleneck,self).__init__()
#         print(downsample)
        self.bottleneck=nn.Sequential(
            nn.Conv2d(in_planes,planes,1,bias=False),
            nn.BatchNorm2d(planes),
            nn.ReLU(inplace=True),
            nn.Conv2d(planes,planes,3,stride,1,bias=False),
            nn.BatchNorm2d(planes),
            nn.ReLU(inplace=True),
            nn.Conv2d(planes,self.expansion*planes,1,bias=False),
            nn.BatchNorm2d(self.expansion*planes),
        )
        self.relu=nn.ReLU(inplace=True)
        self.downsample=downsample
    def forward(self,x):
        identity=x
        out=self.bottleneck(x)
        if self.downsample is not None:
#             print(self.downsample)
            identity=self.downsample(x)
        out+=identity
        out=self.relu(out)
        return out
 
#FNP的类，初始化需要一个list，代表RESNET的每一个阶段的Bottleneck的数量
import torch.nn as nn
import torch.nn.functional as F
import math
class FeatureSelectionModule(nn.Module):
    def __init__(self, in_chan, out_chan, norm="GN"):
        super(FeatureSelectionModule, self).__init__()
        self.conv_atten = Conv2d(in_chan, in_chan, kernel_size=1, bias=False, norm=get_norm(norm, in_chan))
        self.sigmoid = nn.Sigmoid()
        self.conv = Conv2d(in_chan, out_chan, kernel_size=1, bias=False, norm=get_norm('', out_chan))
        weight_init.c2_xavier_fill(self.conv_atten)
        weight_init.c2_xavier_fill(self.conv)

    def forward(self, x):
        atten = self.sigmoid(self.conv_atten(F.avg_pool2d(x, x.size()[2:])))
        feat = torch.mul(x, atten)
        x = x + feat
        feat = self.conv(x)
        return feat


class FeatureAlign_V2(nn.Module):  # FaPN full version
    def __init__(self, in_nc=128, out_nc=128, norm=None):
        super(FeatureAlign_V2, self).__init__()
        self.lateral_conv = FeatureSelectionModule(in_nc, out_nc, norm="")
        self.offset = Conv2d(out_nc * 2, out_nc, kernel_size=1, stride=1, padding=0, bias=False, norm=norm)
        self.dcpack_L2 = dcn_v2(out_nc, out_nc, 3, stride=1, padding=1, dilation=1, deformable_groups=8,
                                extra_offset_mask=True)
        self.relu = nn.ReLU(inplace=True)
        weight_init.c2_xavier_fill(self.offset)

    def forward(self, feat_l, feat_s, main_path=None):
        HW = feat_l.size()[2:]
        if feat_l.size()[2:] != feat_s.size()[2:]:
            feat_up = F.interpolate(feat_s, HW, mode='bilinear', align_corners=False)
        else:
            feat_up = feat_s
        feat_arm = self.lateral_conv(feat_l)  # 0~1 * feats
        offset = self.offset(torch.cat([feat_arm, feat_up * 2], dim=1))  # concat for offset by compute the dif
        feat_align = self.relu(self.dcpack_L2([feat_up, offset], main_path))  # [feat, offset]
        return feat_align + feat_arm
 
#ResNet的基本Bottleneck类
class Bottleneck(nn.Module):
    expansion=4#通道倍增数
    def __init__(self,in_planes,planes,stride=1,downsample=None):
        super(Bottleneck,self).__init__()
#         print(downsample)
        self.bottleneck=nn.Sequential(
            nn.Conv2d(in_planes,planes,1,bias=False),
            nn.BatchNorm2d(planes),
            nn.ReLU(inplace=True),
            nn.Conv2d(planes,planes,3,stride,1,bias=False),
            nn.BatchNorm2d(planes),
            nn.ReLU(inplace=True),
            nn.Conv2d(planes,self.expansion*planes,1,bias=False),
            nn.BatchNorm2d(self.expansion*planes),
        )
        self.relu=nn.ReLU(inplace=True)
        self.downsample=downsample
    def forward(self,x):
        identity=x
        out=self.bottleneck(x)
        if self.downsample is not None:
#             print(self.downsample)
            identity=self.downsample(x)
        out+=identity
        out=self.relu(out)
        return out
 
#FNP的类，初始化需要一个list，代表RESNET的每一个阶段的Bottleneck的数量
class FPN(nn.Module):
    def __init__(self,layers):
        super(FPN,self).__init__()
        self.inplanes=64
        #处理输入的C1模块（C1代表了RestNet的前几个卷积与池化层）
        self.conv1=nn.Conv2d(3,64,7,2,3,bias=False)
        self.bn1=nn.BatchNorm2d(64)
        self.relu=nn.ReLU(inplace=True)
        self.maxpool=nn.MaxPool2d(3,2,1)
        #搭建自下而上的C2，C3，C4，C5
        self.layer1=self._make_layer(64,layers[0])
        self.layer2=self._make_layer(128,layers[1],2)
        self.layer3=self._make_layer(256,layers[2],2)
        self.layer4=self._make_layer(512,layers[3],2)
        #对C5减少通道数，得到P5
        self.toplayer=nn.Conv2d(2048,256,1,1,0)
        #3x3卷积融合特征
        self.smooth1=nn.Conv2d(256,256,3,1,1)
        self.smooth2=nn.Conv2d(256,256,3,1,1)
        self.smooth3=nn.Conv2d(256,256,3,1,1)
        #横向连接，保证通道数相同
        self.latlayer1=nn.Conv2d(1024,256,1,1,0)
        self.latlayer2=nn.Conv2d(512,256,1,1,0)
        self.latlayer3=nn.Conv2d(256,256,1,1,0)
        
        self.pool2 = nn.MaxPool2d(kernel_size=8,stride=8)
        self.pool3 = nn.MaxPool2d(kernel_size=4,stride=4)
        self.pool4 = nn.MaxPool2d(kernel_size=2,stride=2)
        
#         self.lateral_conv = Conv2d(512, 256, kernel_size=1, bias=True,
#                       norm=get_norm(norm, out_channels))
        
    def _make_layer(self,planes,blocks,stride=1):
        downsample=None
        if stride!=1 or self.inplanes != Bottleneck.expansion*planes:
            downsample=nn.Sequential(
                nn.Conv2d(self.inplanes,Bottleneck.expansion*planes,1,stride,bias=False),
                nn.BatchNorm2d(Bottleneck.expansion*planes)
            )
        layers=[]
        layers.append(Bottleneck(self.inplanes,planes,stride,downsample))
        self.inplanes=planes*Bottleneck.expansion
        for i in range(1,blocks):
            layers.append(Bottleneck(self.inplanes,planes))
        return nn.Sequential(*layers)
    #自上而下的采样模块
    def _upsample_add(self,x,y):
        _,_,H,W=y.shape
        return F.upsample(x,size=(H,W),mode='bilinear')+y
#     def forward(self,x):
#         #自下而上
#         c1=self.maxpool(self.relu(self.bn1(self.conv1(x))))
#         c2=self.layer1(c1)
#         c3=self.layer2(c2)
#         c4=self.layer3(c3)
#         c5=self.layer4(c4)
#         #自上而下
#         p5=self.toplayer(c5)
#         p4=self._upsample_add(p5,self.latlayer1(c4))
#         p3=self._upsample_add(p4,self.latlayer2(c3))
#         p2=self._upsample_add(p3,self.latlayer3(c2))
#         #卷积的融合，平滑处理
#         p4=self.smooth1(p4)
#         p3=self.smooth2(p3)
#         p2=self.smooth3(p2)
        
#         return self.pool2(p2),self.pool3(p3),self.pool4(p4),p5
    def forward(self,x):
        #自下而上
        c1=self.maxpool(self.relu(self.bn1(self.conv1(x))))
        c2=self.layer1(c1)
        c3=self.layer2(c2)
        c4=self.layer3(c3)
        c5=self.layer4(c4)
        #自上而下
        p5=self.toplayer(c5)
        p4=self._upsample_add(p5,self.latlayer1(c4))
        p3=self._upsample_add(p4,self.latlayer2(c3))
        p2=self._upsample_add(p3,self.latlayer3(c2))
        #卷积的融合，平滑处理
        p4=self.smooth1(p4)
        p3=self.smooth2(p3)
        p2=self.smooth3(p2)
        
        return p2,p3,p4,p5