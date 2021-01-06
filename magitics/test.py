import numpy as np
import pandas as pd
import os
from Bio import SeqIO



# import pandas as pd
# df=pd.read_csv('/home/ylucas/Bureau/annotated_fastas/1313.5461.PATRIC.features.tab',sep='\t')
# from Bio import SeqIO
# it=SeqIO.parse('/home/ylucas/Bureau/annotated_fastas/1313.5461.fna','fasta')
# ls=[]
# for truc in it:
#     ls.append(truc)

# path='/home/ylucas/magitics_data/cefepime/'
# ls_paths=os.listdir(path)
#
# ls_dfs=[]
# for folder in ls_paths:
#     for file in os.listdir(os.path.join(path, folder)):
#         ls_dfs.append(pd.read_csv(os.path.join(path,folder,file), sep= '\t'))
# df=pd.concat(ls_dfs, ignore_index=True)

#a=df['plfam_id'].value_counts()
#a=df['plfam_id'].hist(bins=10)

# with open('/home/ylucas/magitics_data/ls_file.txt', 'r') as f:
#     ls_file=f.readlines()
#
# with open('/home/ylucas/magitics_data/ls_aligns.txt','r') as f:
#     ls_aligns=f.readlines()
#
# ls_missing=[]
# for i in ls_file:
#     if i not in ls_aligns:
#         ls_missing.append(i)



from BCBio import GFF
class parse_gff_fasta_to_extract_CDS(object):
    def __init__(self):
        return

    def parse_gff(self, pathgff):
        import gffutils
        db = gffutils.create_db(data=pathgff, dbfn='memory',
                                force=True)
        return db

    def parse_fasta(self, pathfasta):
        import pyfaidx
        fasta = pyfaidx.Fasta(pathfasta)


    def run(self, path):
        pathgff=os.path.join(path, 'gnagna.gff')
        pathfasta=os.path.join(path, 'gnagna.fna')

        db=self.parse_gff(pathgff)
        fasta=self.parse_fasta(pathfasta)

        ls_seq = []
        ls_id=[]
        for cds in db.features_of_type('CDS', order_by='start'):
            ls_seq.append(cds.sequence(fasta))
            ls_id.append(cds.id)



# df=pd.read_csv('/home/ylucas/Bureau/annotated_fastas/1313.5461.PATRIC.features.tab',sep='\t')
# from Bio import SeqIO
# it=SeqIO.parse('/home/ylucas/Bureau/annotated_fastas/1313.5461.fna','fasta')
# ls=[]
# for truc in it:
#     ls.append(truc)

#GFF_fasta=parse_gff_fasta_to_extract_CDS()
#GFF_fasta.parse_gff('/home/ylucas/magitics_data/Prokka_annotation/PROKKA_11092020.gff')


"""
from collections import defaultdict
import pandas as pd
from Bio import SearchIO

filename = '/home/ylucas/magitics_data/hmmscan_out.txt'

attribs = ['accession', 'bias', 'bitscore', 'description', 'cluster_num', 'domain_exp_num',  'domain_included_num', 'domain_obs_num', 'domain_reported_num', 'env_num', 'evalue', 'id', 'overlap_num', 'region_num']

hits = defaultdict(list)

with open(filename) as handle:
    for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
      #print(queryresult.id)
      #print(queryresult.accession)
      #print(queryresult.description)
      for hit in queryresult.hits:
        for attrib in attribs:
          hits[attrib].append(getattr(hit, attrib))

df=pd.DataFrame.from_dict(hits)

import pandas as pd
path='/home/ylucas/magitics_data/cef_hmm_pred.txt'

df= pd.read_csv(path, sep='\t')
file=open(path,'r')
lines=file.readlines()
print(lines[:10])
"""

"""
path='/home/ylucas/magitics_data/cef_hmm_pred.txt'

from collections import defaultdict
import pandas as pd
from Bio import SearchIO
attribs = ['query_id', 'bitscore','evalue', 'id']

hits = defaultdict(list)

with open(path) as handle:
    for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
      #print(queryresult.id)
      #print(queryresult.accession)
      #print(queryresult.description)
      for hit in queryresult.hits:
        for attrib in attribs:
          hits[attrib].append(getattr(hit, attrib))

df=pd.DataFrame.from_dict(hits)

grouped=df.groupby('query_id')
uniques=df['query_id'].unique()

ls_evalue=[]
ls_score=[]
ls_ratio=[]

for locus_id in uniques:
    group=grouped.get_group(locus_id)
    try:
        ls_evalue.append(np.log10(group.iloc[0,2])-np.log10(group.iloc[1,2]))
        ls_score.append(np.log10(group.iloc[0, 2]))
        ls_ratio.append(group.iloc[0, 1]-group.iloc[1, 1])
    except:
        a=1

import matplotlib.pyplot as plt
ls_a=[]
for thing in ls_evalue:
    if thing>-3000:
        ls_a.append(thing)


c=0
for thing in ls_a:
    if thing<2:
        c+=1
plt.hist(ls_a, bins=100)
"""
import pandas as pd
import numpy as np
import argparse
from collections import defaultdict
from Bio import SearchIO, SeqIO

class hmmscan_to_pandas(object):
    def __init__(self):
        return

    def parse_file(self, path):
        attribs = ['query_id', 'bitscore', 'evalue', 'id']

        hits = defaultdict(list)
        with open(path) as handle:
            for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
              #print(queryresult.id)
              #print(queryresult.accession)
              #print(queryresult.description)
              for hit in queryresult.hits:
                for attrib in attribs:
                  hits[attrib].append(getattr(hit, attrib))
        df=pd.DataFrame.from_dict(hits)
        return df

    def get_score_values(self, df):
        grouped = df.groupby('query_id')
        uniques = df['query_id'].unique()

        ls_evalue = []
        ls_score = []
        ls_ratio = []

        for locus_id in uniques:
            group = grouped.get_group(locus_id)
            try:
                ls_evalue.append(np.log10(group.iloc[0, 2]) - np.log10(group.iloc[1, 2]))
                ls_score.append(group.iloc[0, 1])
                ls_ratio.append(group.iloc[0, 1] - group.iloc[1, 1])
            except:
                a = 1


    def get_annotation(self, df):
        grouped = df.groupby('query_id')
        uniques = df['query_id'].unique()

        dic_annotation={}
        for locus_id in uniques:
            group = grouped.get_group(locus_id)
            dic_annotation[group.iloc[0,0]]=group.iloc[0,3] #todo verifier coordinates
        return dic_annotation

pathhmmscan='/home/ylucas/287999_hmmscan.txt'
hmmscan_parse=hmmscan_to_pandas()
parsed_hmmscan=hmmscan_parse.parse_file(pathhmmscan)
dic_annotation=hmmscan_parse.get_annotation(parsed_hmmscan)
ls_scan=[]
for key in dic_annotation:
    ls_scan.append(dic_annotation[key])


pathpatric='/home/ylucas/287.999.PATRIC.features.tab'
df=pd.read_csv(pathpatric, sep='\t')

# Get names of indexes for which column Age has value 30
indexNames = df[ df['feature_type'] != 'CDS' ].index
# Delete these row indexes from dataFrame
df.drop(indexNames , inplace=True)
ls_pat=df['plfam_id'].to_list()

import difflib
sm=difflib.SequenceMatcher(None, ls_pat, ls_scan)
sm.ratio()

ls_wrong=[]

for thing in ls_pat:
    if thing in ls_scan:
        ls_scan.remove(thing)

grouped=parsed_hmmscan.groupby('query_id')
dic_wrong={}
dic_wrong_noratio={}
for wrongplf in ls_scan:
    susceptible=parsed_hmmscan.loc[parsed_hmmscan['id']==wrongplf]
    for index, row in susceptible.iterrows():
        CDS=row['query_id']
        group=grouped.get_group(CDS)
        if group.iloc[0,3]==wrongplf:
            dic_wrong_noratio[wrongplf] = [-np.log10(group.iloc[0, 2])]

            try:
                dic_wrong[wrongplf]=[-np.log10(group.iloc[0,2]), -np.log10(group.iloc[0,2])+np.log10(group.iloc[1,2])]
            except:
                a=0

dfwrong=pd.DataFrame.from_dict(dic_wrong, orient='index')
dfwrongreplaced=dfwrong.replace([np.inf, -np.inf], np.nan)
dfwrongreplaced.dropna(inplace=True)

fasta=SeqIO.parse('/home/ylucas/287.999.fna','fasta')
ls_contigs=[]
for record in fasta:
    ls_contigs.append(record)
