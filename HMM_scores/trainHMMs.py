
#TODO 1 PARSE GFF+FASTA POUR FAIRE LISTE DE SUPERVISED PLFAM SEQ
import pandas as pd
import os, sys
from Bio import SeqIO
import numpy as np

class gene_family_file(object):
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

    def write_plfam_fastas(self, pathtofasta, pathtogff, label):
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
                with open(os.path.join('../toy_dataset/plfams/susceptible',plfam), 'a') as f:
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
            pathtofasta=os.path.join(pathtofolder, strainID+'.fna')
            pathtogff=os.path.join(pathtofolder, strainID+'.PATRIC.features.tab')
            label=self.get_label(pathtofolder)
            self.write_plfam_fastas(pathtofasta, pathtogff, label)

#gene_fam=gene_family_file()
#gene_fam.run()

#TODO 2: MSA + HMM PROFILE PAR PLFAM
"""
Execute the following bash commands:
* cd plfams/resistant/
* ls > ../ls_plfams_resistant.txt 
* cd ../..
* mkdir aligns
* cd aligns
* mkdir resistant
* mkdir susceptible
* cd ..
* mkdir hmm_profiles
* cd hmm_profiles
* mkdir resistant
* mkdir susceptible
* cd ..
write the following bash script using nano and save it under the name run.sh

run the script using:
* sbatch run.sh
------------------
#!/bin/bash
#SBATCH --job-name=musclehmmer
#SBATCH --nodes=1
#SBATCH --mem=8gb

#SBATCH --ntasks-per-node=1
#SBATCH --array=1-10000%100  #The range of the array will need to be updated 
                             # in order to cover the totality of the plfams collection

module load hmmer
module load muscle

FILENAME=$(head -n $(($SLURM_ARRAY_TASK_ID)) ls_plfams_resistant.txt |tail -1)

muscle -in plfams/$FILENAME -out aligns/$FILENAME
hmmbuild hmm_profiles/$FILENAME aligns/$FILENAME
------------------
"""

#TODO 3: HMM FLATFILE PAR COUPLE DE PROFILE DE PLFAM
#1. Get list of existing hmm_profiles
def write_plfam_lists(pathresistant='../toy_dataset/plfams/resistant', pathsusceptible='../toy_dataset/plfams/susceptible', pathtowrite='../toy_dataset/plfams'):
    ls_resistant_file=os.listdir(pathresistant)
    ls_resistant=[]
    for thing in ls_resistant_file:
        ls_resistant.append(thing.split('_')[-1])
    ls_susceptible=os.listdir(pathsusceptible)

    susceptible_only=open(os.path.join(pathtowrite, 'susceptible_only.txt'),'w')
    resistant_only=open(os.path.join(pathtowrite, 'resistant_only.txt'),'w')
    both_susc_and_res=open(os.path.join(pathtowrite, 'both_susc_and_res.txt'),'w')
    for thing in ls_susceptible:
        if thing in ls_resistant:
            both_susc_and_res.write(thing+'\n')
        else:
            susceptible_only.write(thing+'\n')

    for thing in ls_resistant:
        if thing not in ls_susceptible:
            resistant_only.write(thing+'\n')


    susceptible_only.close()
    resistant_only.close()
    both_susc_and_res.close()

#1.5 move and rename hmm_profiles to the same folder: susceptiblePLFXXX and resistantPLFXXX
"""
run the script using:
* sbatch run.sh
------------------
#!/bin/bash
#SBATCH --job-name=rename_profile
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-10000%100  #The range of the array will need to be updated 
                             # in order to cover the totality of the plfams collection


FILENAME=$(head -n $(($SLURM_ARRAY_TASK_ID)) ls_resistant.txt |tail -1)


mv resistant/$FILENAME all_profiles/resistant_$FILENAME
------------------
"""

#2 concatenate hmm_profiles of the same plfam
def concatenate_profiles(path)

"""
Execute the following bash commands:

* mkdir flatfiles
------------------
#!/bin/bash
#SBATCH --job-name=musclehmmer
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-10000%100  #The range of the array will need to be updated 
                             # in order to cover the totality of the plfams collection

module load hmmer
module load muscle

FILENAME=$(head -n $(($SLURM_ARRAY_TASK_ID)) susceptible_only.txt |tail -1)

muscle -in plfams/$FILENAME -out aligns/$FILENAME
hmmbuild hmm_profiles/$FILENAME aligns/$FILENAME
------------------

"""

