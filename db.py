from pymongo import MongoClient # Database connector
import pandas as pd

client = MongoClient('localhost', 27017)    #Configure the connection to the database
db = client.vispr    #Select the database
sequence = db.sequence #Select the collection
annotation = db.annotation


## Insert gene sequence
with open('ten_gene.txt') as f:
    genes = f.read().splitlines()

gene_info = pd.read_table('gene_info.txt', header=0, sep='\t')
ids = list(gene_info['ID'][gene_info['Gene'].isin(genes)])

df_seq = pd.read_table('ten_gene_seq.txt', header=None, sep='\t', index_col=0)
df_seq.columns = ['Seq']

df_gene = pd.read_table('gene_info.txt', header=0, sep='\t', index_col=3)
for id in ids:
    sequence.insert_one({"id":id, "name":df_gene.loc[id,'Gene'], "chr":df_gene.loc[id,'Chr'],
                         "start":str(df_gene.loc[id,'Start']), "stop":str(df_gene.loc[id,'End']),
                         "seq": df_seq.loc[id,'Seq']})


## Insert gene sequence whole database
gene_info = pd.read_table('gene_info.txt', header=0, sep='\t', index_col=3)
with open('gene_seq.txt', 'r') as fin:
    line = fin.readline() #file header
    line = fin.readline()
    while line:
        line = line.split()
        sequence.insert_one({"id": line[0], "name": gene_info.loc[line[0], 'Gene'], "chr": gene_info.loc[line[0], 'Chr'],
                             "start": str(gene_info.loc[line[0], 'Start']), "stop": str(gene_info.loc[line[0], 'End']),
                             "seq": line[1]})
        line = fin.readline()


## Insert gene annotation
df_major_trans = pd.read_table('main_transcripts.txt', header=0, sep='\t')

with open('ten_gene_annotation.txt', 'r') as fin:
    line = fin.readline()
    while line[0] == '#':  # skip file header
        line = fin.readline()
    line = line.split()
    while line:
        if line[2] == 'gene':
            gene = line[line.index('gene_id')+1].split("\"")[1]
            gene_name = line[line.index('gene_name')+1].split("\"")[1]
            line = fin.readline()
            line = line.split()
            while 'gene_id' in line and line[line.index('gene_id')+1].split("\"")[1] == gene:
                if line[2] == 'transcript':
                    transcript = line[line.index('transcript_id')+1].split("\"")[1]
                else:
                    line = fin.readline()
                    line = line.split()
                    continue
                if transcript not in list(df_major_trans['Trans_ID']):
                    line = fin.readline()
                    line = line.split()
                    continue
                exon_id = []
                exon_start = []
                exon_stop = []
                line = fin.readline()
                line = line.split()
                while 'transcript_id' in line and line[line.index('transcript_id') + 1].split("\"")[1] == transcript:
                    if line[2]=='exon':
                        exon_id.append(line[line.index('exon_id') + 1].split("\"")[1])
                        exon_start.append(line[3])
                        exon_stop.append(line[4])
                    line = fin.readline()
                    line = line.split()
                exon = {"gene_id":gene, "gene_name":gene_name, "transcript_id":transcript, "exon_id":exon_id,
                        "exon_start":exon_start, "exon_stop":exon_stop}
                #print(exon)
                annotation.insert_one(exon)
