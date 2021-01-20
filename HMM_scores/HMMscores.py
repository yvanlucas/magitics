#TODO:
# For each strain:
# 1) Gather PLFAM ids and corresponding nuc sequence
# 2) HMMscan each locus against the corresponding HMM flatfile
# 3) get bitscore (susceptible - resistant) or only the one existings
# 4) put bitscore into a dic for each strain: strain: {plfam1:score1, ...}
# 5) Create dataframe


# old data path '../toy_dataset/susceptible'
# new data path '/scratch/MAGITICS_data/Pseudomonas_aeruginosa/levofloxacin/data/test/susceptible'

datapath= '/scratch/MAGITICS_data/Pseudomonas_aeruginosa/levofloxacin/data/test/susceptible'
import pandas as pd
import os, sys
from Bio import SeqIO
import numpy as np

# 1) Gather PLFAM ids and corresponding nuc sequence
class parse_PATRIC_and_fastas(object):
    def __init__(self):
        return

    def parse_patric_gff(self, pathtofile):
        print(pathtofile)
        df=pd.read_csv(pathtofile, sep='\t')
        return df

    def dic_data(self,row, contig_number):
        dic={}
        dic['contig']=contig_number -1
        dic['start']=int(row['start'])
        dic['end']=int(row['end'])
        dic['patric_id']=str(row['patric_id'])
        dic['plfam_id']=str(row['plfam_id'])
        dic['genome']=str(row['genome_name'])
        dic['genome_id']=str(row['genome_id'])
        return dic

    def get_locus_plfams(self, df):
        """
        get locus of gene (start end), contig number, patric_id and plfam if feature_type=='CDS'
        """
        ls=[]
        contig_number=0
        prev_contig=''
        for index, row in df.iterrows():
            if not row['accession'] == prev_contig:
                contig_number+=1
            prev_contig=row['accession']
            if row['feature_type']=='CDS':
                ls.append(self.dic_data(row, contig_number))
        return ls

    def write_plfam_fastas(self, pathtofasta, pathtogff, pathstraindir, label):
        df=self.parse_patric_gff(pathtogff)
        ls_dics=self.get_locus_plfams(df)

        fasta=SeqIO.parse(pathtofasta, 'fasta')
        ls_contigs=[]
        for contig in fasta:
            ls_contigs.append(contig)
        for dic in ls_dics:
            try:
                plfam=dic['plfam_id']
                prot_seq=ls_contigs[dic['contig']][dic['start']-1:dic['end']-1] #.translate()
                write_seq=''.join([aa for aa in prot_seq])
                with open(os.path.join(pathstraindir ,plfam), 'a') as f:
                    f.write('>'+dic['genome']+'|'+dic['patric_id']+'\n')
                    f.write(write_seq)
                    f.write('\n')
                    self.count+=1
            except Exception as e:
                print(e)
                print(dic['plfam_id'])
                print(dic['contig'])
                print(dic['start'])
                print(dic['end'])
                print(len(ls_contigs[dic['contig']]))

    def unique_listdir_PATRIC(self, filelist):
        ls_unique=[]
        for thing in filelist:
            ls_unique.append('.'.join(thing.split('.')[:2]))
        return list(set(ls_unique))

    def get_label(self, pathtofolder):
        label=pathtofolder.split('/')[-1]
        return label

# old path '../toy_dataset/susceptible'
#newpath '/scratch/MAGITICS_data/Pseudomonas_aeruginosa/levofloxacin/data/test/susceptible'
    def run(self, pathtofolder='../toy_dataset/susceptible'): #TODO: change path
        self.count=0
        strainlist=self.unique_listdir_PATRIC(os.listdir(pathtofolder))
        for strainID in strainlist:
            pathstraindir=os.path.join(pathtofolder, strainID)
            os.system("mkdir %s" % (pathstraindir))
            pathtofasta=os.path.join(pathtofolder, strainID+'.fna')
            pathtogff=os.path.join(pathtofolder, strainID+'.PATRIC.features.tab')
            label=self.get_label(pathtofolder)
            self.write_plfam_fastas(pathtofasta, pathtogff, pathstraindir, label)

gene_fam=parse_PATRIC_and_fastas()
gene_fam.run()

#######################################################################################
# 2) HMMscan each locus against the corresponding HMM flatfile
#pour compréhension: hmmscan --tblout PATRIC_scan_out.txt --cpu 10 hmm_database_flatfile 287999_new/PROKKA_12082020.ffn

import argparse
def get_strainname():
    parser=argparse.ArgumentParser()
    parser.add_argument('--strainname', type=str, default='31')
    arg=parser.parse_args()
    return arg.strainname


def hmm_scan_each_strain(strainname):

    pathtostrain=os.path.join('../toy_dataset/susceptible', strainname) #TODO: change path
    ls_plfam=os.listdir(pathtostrain)
    pathtoflatfiles=os.path.join('../toy_dataset/susceptible') #TODO: change path

    for plfam in ls_plfam:
        pathtoscan=os.path.join(pathtostrain,'hmmscan', plfam)
        cmd="hmmscan --tblout %s.txt --cpu 1 %s %s" %( pathtoscan, os.path.join(pathtoflatfiles, plfam), pathtoplfam+'.ffn')
        #TODO ^ change path of nuc sequence to scan
        os.system(cmd)

#FOR EACH STRAINNAME
strainname=get_strainname()

hmm_scan_each_strain(strainname)


# 3) Get Scores and label from HMMscan output
from Bio import SearchIO
from collections import defaultdict
import pickle
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

        dic={}

        for locus_id in uniques:
            group = grouped.get_group(locus_id)
            ls_id = []
            ls_score = []
            try:
                ls_score.append(group.iloc[0, 1])
                ls_score.append(group.iloc[1,1])
                ls_id.append(group.iloc[0,3])
                ls_id.append(group.iloc[1,3])
                score, plf=self.calculate_score(ls_score, ls_id)
                dic[plf]=score
            except:
                a = 1
        return dic

    def get_label(self, strainID):
        label='susceptible' #TODO: change by hand or automate it
    def calculate_score(self, ls_score, ls_id):
        if ls_id[0][:4]=='susc':
            score=ls_score[0]-ls_score[1]
            plf=ls_id[1][10:]
        elif ls_id[0][:4]=='resi':
            score=np.float64(ls_score[1])-np.float64(ls_score[0])
            plf=ls_id[0][10:]
        return score, plf


    def run(self, strainID):
        df_annotations=self.parse_file(os.path.join('../toy_dataset/susceptible',strainID,'hmmscan')) #TODO: change path
        dic=self.get_score_values(df_annotations)
        dic['label']=self.get_label(strainID)
        with open(os.path.join('../toy_dataset/susceptible', strainID, 'dic_scores.pkl'), "w") as f: #TODO: change path
            pickle.dump(dic, f)


#FOR EACH STRAINNAME
hmmscanparse=hmmscan_to_pandas()

hmmscanparse.run(strainname)
#4 Load all dics and append them in a pandas dataframe


class create_df_from_hmmscandict(object):
    def __init__(self):
        return

    def run(self, datadir):
        ls_straindir=os.listdir(datadir)

        ls_dicts=[]
        for straindir in ls_straindir:
            with open(os.path.join(datadir, straindir, 'dic_scores.pkl'), "r") as f:  # TODO: change path
                dic=pickle.load(f)
            ls_dicts.append(dic)

        self.df=pd.dataframe(ls_dicts)


create_df=create_df_from_hmmscandict()
create_df.run('../toy_dataset_susceptible')









"""

#!/bin/bash

#SBATCH --job-name=hmm_profile
#SBATCH --nodes=1
#SBATCH --mem=60gb
#SBATCH --ntasks-per-node=12
#SBATCH --array=1

module load hmmer

hmmscan --tblout PATRIC_scan_out.txt --cpu 10 hmm_database_flatfile 287999_new/PROKKA_12082020.ffn

"""