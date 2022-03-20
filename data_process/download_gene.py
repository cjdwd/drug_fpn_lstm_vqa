import os

dataset = 'kd'
text = open('./'+dataset+'_uniprot_gene_map.txt').readlines()
id_map = {}
for line in text:
    uniprot_id, gene_id = line.strip().split()
    id_map[uniprot_id] = gene_id

id_existed = []
files = os.listdir('../data/gene_seq/'+dataset+'/')
for file in files:
    id_existed.append(file.strip().split('.')[0])

for uniprot_id in id_map.keys():
    print(uniprot_id)
    if uniprot_id not in id_existed:
        os.system(
            'esearch -db gene -query "' + id_map[uniprot_id] +
            ' [ID]" | efetch -format docsum | xtract -pattern '
            'GenomicInfoType -element ChrAccVer ChrStart ChrStop '
            '| xargs -n 3 sh -c \'efetch -db nuccore -format fasta '
            '-id "$0" -chr_start "$1" -chr_stop "$2"\' > '
            '../data/gene_seq/' + dataset +'/' + uniprot_id + '.txt')
