import pandas as pd
import numpy as np
import argparse
from collections import defaultdict
from Bio import SearchIO

####################################################################
#TODO 1. Parse existing gff using gff_to_pandas() class

class gff_to_pandas(object):
    def __init__(self):
        return

    def gffs_to_dic(self, pathtogff):
        with open(pathtogff, 'r') as f:
            lines=f.readlines()

        dic={'scaffold': [], 'source': [], 'type': [], 'start': [], 'end': [],
         'score': [], 'strand': [], 'frame': [], 'attributes': [], 'ID':[], 'eC_number':[],'Name':[],
             'db_xref':[], 'gene':[], 'inference':[], 'locus_tag':[], 'product':[], 'note':[], 'rpt_family':[],
             'rpt_type':[]}

        tabsplit=['scaffold','source','type','start','end','score','strand','frame','attributes']
        for line in lines[1:]:
            #print(line)
            if len(line.split('\t'))<3:
                continue

            line=line.split('\t')
            for thing, head in zip(line,tabsplit):

                dic[head].append(thing)
            attributes=dic['attributes'][-1].split(';')
            IDTHERE=0
            for thing in attributes:
                splitted=thing.split('=')
                dic[splitted[0]].append(splitted[1])
                if splitted[0]=='ID':
                    IDTHERE=1
            if IDTHERE==0:
                print(attributes)
                keyslist=list(dic.keys())
                keyslist.remove('ID')
                for key in keyslist:
                    del dic[key][-1]
        #     inference=dic['inference'][-1].split(':')
        #     if len(inference)==5:
        #         dic['prot_fam'].append(inference[-1])
        #     else:
        #         dic['prot_fam'].append('protein')
        todel=['eC_number', 'Name','gene','inference','note', 'db_xref','locus_tag','product', 'rpt_family', 'rpt_type']
        for thing in todel:
            del(dic[thing])
        return dic

    def lsdic_to_pandas(self, dic):
        for key in dic.keys():
            print(key)
            print(len(dic[key]))
        #dfs=[]
        # for key in ls_dic[0]:
        #     print(key)
        #     print(len(ls_dic[0][key]))
        #for dic in ls_dic:
        #    dfs.append(pd.DataFrame.from_dict(dic))

        #dfs=pd.concat(dfs, ignore_index=True)
        df=pd.DataFrame.from_dict(dic)
        return df

    def run(self, path):
        ls_dic=self.gffs_to_dic(path)
        dfs=self.lsdic_to_pandas(ls_dic)
        types = {'start': np.int64, 'end': np.int64}
        dfs = dfs.astype(types)
        dfs['size']=dfs['end']-dfs['start']
        return dfs


##################################################################################
#TODO 2. Parse hmmscan output to get locus <-> annotation correspondance

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


##############################################################################################

#todo 3. Add annotation as a column to the gff parsed as pandas and save dataframe as csv
class annotate_and_save(object):
    def __init__(self):
        return

    def run(self, df, dic_annotation):

        df['PLFAM']=0
        for key in dic_annotation:
            df.loc[df['ID']==key, 'PLFAM']=dic_annotation[key]

        return df




def get_lenkmers():
    parser=argparse.ArgumentParser()
    parser.add_argument('--pathgff', type=str, default='/home/ylucas/magitics_data/NMDC60014126_prokka/NMDC60014126.gff')
    parser.add_argument('--pathhmmscan', type=str, default='/home/ylucas/magitics_data/NMDC60014126_prokka/cef_hmm_pred.txt')
    parser.add_argument('--outpath', type=str, default='/home/ylucas/magitics_data/NMDC60014126_prokka/PLFAM_annotation.csv')
    arg=parser.parse_args()
    return arg

if __name__=='__main__':
    args= get_lenkmers()

    pathgff=args.pathgff
    gff_parse=gff_to_pandas()
    parsed_gff=gff_parse.run(pathgff)

    pathhmmscan=args.pathhmmscan
    hmmscan_parse=hmmscan_to_pandas()
    parsed_hmmscan=hmmscan_parse.parse_file(pathhmmscan)
    dic_annotation=hmmscan_parse.get_annotation(parsed_hmmscan)

    annotation=annotate_and_save()
    parsed_gff=annotation.run(parsed_gff, dic_annotation)

    parsed_gff.to_csv(args.outpath, index=False)




