"""
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import log_loss

X_all = np.random.randn(5000, 2)
y_all = ((X_all[:, 0]/X_all[:, 1]) > 1.5)*2 - 1
X_train, X_test, y_train, y_test = train_test_split(X_all, y_all, test_size=0.5, random_state=42)


clf = GradientBoostingClassifier(n_estimators=10, learning_rate=0.01, max_depth=3, random_state=0)
clf.fit(X_train, y_train)
y_pred = clf.predict_proba(X_test)[:, 1]
print('Accuracy for a GBM: {}'.format(clf.score(X_test, y_test)))
print("Test logloss: {}".format(log_loss(y_test, y_pred)))

def compute_loss(y_true, scores_pred):
    return log_loss(y_true, sigmoid(scores_pred))

def sigmoid(x):
    return 1 / (1 + np.exp(-x))

cum_preds = np.array([x for x in clf.staged_decision_function(X_test)])[:, :, 0]

def prune_prediction(cum_pred, ls_index):
    preds_out=cum_preds[-1,:]
    for i in ls_index: #i can't be 0 but who would prune first tree of boosting
        preds_out=preds_out - (cum_preds[i-1,:]-cum_preds[i,:])
    return preds_out


import pandas as pd

df=pd.read_csv('/home/ylucas/Bureau/expe_postdoc/data_postdoc/287.846.PATRIC.features.tab', sep='\t')
a=  df.describe(include='all')

dfbis=pd.read_csv('/home/ylucas/Bureau/expe_postdoc/data_postdoc/287.847.PATRIC.features.tab', sep='\t')
b=  dfbis.describe(include='all')

dfter=df+dfbis
c=dfter.describe(include='all')
"""

# with open('/home/ylucas/toydata_pseudomonas_levofloxacin/traindata/Resistant287.8519.fa','r') as f:
#     lines=f.readlines()
# ls_lengths=[]
#
# len_contig=0
# for line in lines:
#     if(len(line))>10 & (len(line))<100:
#         len_contig+= len(line)
#     elif (len(line))<10:
#         ls_lengths.append(len_contig)
#         len_contig=0
# print(ls_lengths)

#TODO doc louvic example
# POUR APPELER:
# .. automodule:: deep_search.loaders.query_title_loader
#    :members: build_left_right_matches
"""
   Split training dataset in left and right part. Build the corresponding matches.
   Args:
       data (pd.Dataframe): Dataframe holding positive examples
       left_cols (array-like): Array of columns specific to the left part
       right_cols (array-like): Array of columns specific to the right part
   Returns:
       left (pd.DataFrame): The left part (without duplicates).
       right (pd.DataFrame): The right part (without duplicates).
       matches (np.NDArray): Array holding le the left (0-th column) and right (1-th column) couples.
   Examples:
       >>> data = pd.DataFrame(
       ...    [
       ...        [1, 2, 3, 4],
       ...        [2, 2, 4, 4],
       ...        [1, 2, 5, 4],
       ...        [4, 2, 6, 4],
       ...    ],
       ...    columns=[
       ...        "left_col1",
       ...        "left_col2",
       ...        "right_col1",
       ...        "right_col2",
       ...    ]
       ... )
       >>> left_cols = ["left_col1", "left_col2"]
       >>> right_cols = ["right_col1", "right_col2"]
       >>> build_left_right_matches(data, left_cols, right_cols)
       (
           left_col1  left_col2
       0          1          2
       1          2          2
       3          4          2,
           right_col1  right_col2
       0           3           4
       1           4           4
       2           5           4
       3           6           4,
           array([[0, 0],
                  [1, 1],
                  [0, 2],
                  [2, 3]])
       )
   """
# import pandas as pd
# import numpy as np
# import os
# from collections import Counter
# path='/home/ylucas/Bureau/Data/levofloxacin/'
# dfs= []
# for file in os.listdir(path):
#     dfs.append(pd.read_csv(os.path.join(path, file), sep='\t'))
#
# dfs=pd.concat(dfs, ignore_index=True)
#

# dfs['size']=dfs['end']-dfs['start']
#
# for plfam in dfs['plfam_id'][:50]:
#     a=dfs['size'][dfs['plfam_id']==plfam]
#     b,c = np.unique(a, return_counts=True)
#     print(b)
#     print(c)
#     print(np.sum(c))
#     print('---')

# from BCBio import GFF
# from sgt import SGT
#
# infile='/home/ylucas/Bureau/annotated_fastas/1313.5461.gff'
# inhandle=open(infile)
# for rec in GFF.parse(inhandle):
#     print(rec)


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
        import pandas as pd
        dfs=[]
        for key in ls_dic[0]:
            print(key)
            print(len(ls_dic[0][key]))
        for dic in ls_dic:
            dfs.append(pd.DataFrame.from_dict(dic))
        print(dfs)
        self.dfs=pd.concat(dfs, ignore_index=True)
        return dfs

    def run(self):
        ls_path=['/home/ylucas/Bureau/annotated_fastas/1313.5461.gff','/home/ylucas/Bureau/annotated_fastas/1313.5465.gff','/home/ylucas/Bureau/annotated_fastas/1313.5466.gff']

        ls_dic=[]

        for path in ls_path:
            ls_dic.append(self.gffs_to_dic(path))

        dfs=self.lsdic_to_pandas(ls_dic)
        self.dic=ls_dic

a=gff_to_pandas()
a.run()

