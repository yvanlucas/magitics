#TODO:
# For each strain:
# 1) Gather PLFAM ids and corresponding nuc sequence
# 2) HMMscan each locus against the corresponding HMM flatfile
# 3) get bitscore (susceptible - resistant) or only the one existings
# 4) put bitscore into a dic for each strain: strain: {plfam1:score1, ...}
# 5) Create dataframe

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

    def run(self, pathtofolder='../toy_dataset/susceptible'):
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

# 2) HMMscan each locus against the corresponding HMM flatfile
import argparse
def get_lenkmers():
    parser=argparse.ArgumentParser()
    parser.add_argument('--strainname', type=int, default=31)
    arg=parser.parse_args()
    return arg.strainname


def hmm_scan_each_strain():
    strainname=get_lenkmers()

    pathtostrain=os.path.join('../toy_dataset/susceptible', strainname)
    ls_plfam=os.listdir(pathtostrain)
    pathtoflatfiles=os.path.join('../toy_dataset/susceptible')

    for plfam in ls_plfam:
        pathtoplfam=os.path.join(pathtostrain, plfam)
        cmd="hmmscan --tblout %s.txt --cpu 1 %s %s" %( pathtoplfam, os.path.join(pathtoflatfiles, plfam), pathtoplfam+'.ffn')











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