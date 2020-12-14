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
                with open(os.path.join('/scratch/MAGITICS_data/Pseudomonas/plfams_nuc',plfam), 'a') as f:
                    f.write('>'+dic['genome']+'|'+dic['patric_id']+'|Amikacin: susceptible \n')
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

    def run(self, pathtofolder='/scratch/MAGITICS_data/Pseudomonas/amikacin/Susceptible'):
        self.count=0
        strainlist=self.unique_listdir_PATRIC(os.listdir(pathtofolder))

        for strainID in strainlist:
            pathtofasta=os.path.join(pathtofolder, strainID+'.fna')
            pathtogff=os.path.join(pathtofolder, strainID+'.PATRIC.features.tab')
            label=self.get_label(pathtofolder)
            self.write_plfam_fastas(pathtofasta, pathtogff, label)

gene_fam=gene_family_file()
gene_fam.run()

class gff_to_pandas(object):
    def __init__(self):
        return

    def gffs_to_dic(self, pathtogff):
        with open(pathtogff, 'r') as f:
            lines=f.readlines()

        dic={'scaffold': [], 'source': [], 'type': [], 'start': [], 'end': [],
         'score': [], 'strand': [], 'frame': [], 'attributes': [], 'ID': [], 'eC_number':[], 'Name': [], 'gene': [],
         'inference': [], 'locus_tag': [],
         'product': [], 'prot_fam': [],
         'note': []}
        tabsplit=['scaffold','source','type','start','end','score','strand','frame','attributes']
        for line in lines[1:]:
            #print(line)
            if len(line.split('\t'))<3:
                continue

            line=line.split('\t')
            for thing, head in zip(line,tabsplit):

                dic[head].append(thing)
            attributes=dic['attributes'][-1].split(';')
            for thing in attributes:
                splitted=thing.split('=')
                dic[splitted[0]].append(splitted[1])
            inference=dic['inference'][-1].split(':')
            if len(inference)==5:
                dic['prot_fam'].append(inference[-1])
            else:
                dic['prot_fam'].append('protein')
        todel=['eC_number', 'Name','gene','inference','note']
        for thing in todel:
            del(dic[thing])
        return dic

    def lsdic_to_pandas(self, ls_dic):
        dfs=[]
        for key in ls_dic[0]:
            print(key)
            print(len(ls_dic[0][key]))
        for dic in ls_dic:
            dfs.append(pd.DataFrame.from_dict(dic))
        dfs=pd.concat(dfs, ignore_index=True)
        return dfs

    def describe(self, dfs):
        a=dfs.describe(include='all')
        count3=0
        for plfam in np.unique(dfs['prot_fam']):
            a=dfs['size'][dfs['prot_fam']==plfam]
            b,c = np.unique(a, return_counts=True)
            print(plfam)
            print(b)
            print(c)
            print(np.sum(c))
            print('---')
            if len(b)==1:
                count3+=1
        print(count3)
        print(len(np.unique(dfs['prot_fam'])))
        return a

    def run(self):
        ls_path=['/home/ylucas/Bureau/annotated_fastas/1313.5461.gff','/home/ylucas/Bureau/annotated_fastas/1313.5465.gff','/home/ylucas/Bureau/annotated_fastas/1313.5466.gff']
        ls_dic=[]

        for path in ls_path:
            ls_dic.append(self.gffs_to_dic(path))
        dfs=self.lsdic_to_pandas(ls_dic)
        types = {'start': np.int64, 'end': np.int64}
        dfs = dfs.astype(types)
        print(dfs.dtypes)
        dfs['size']=dfs['end']-dfs['start']
        self.dfs = dfs
        self.dic=ls_dic
        self.desc=self.describe(dfs)


#parse=gff_to_pandas()
#parse.run()