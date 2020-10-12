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
import pandas as pd
import numpy as np
df=pd.read_excel('/home/ylucas/Bureau/SP_strains_metadata.xlsx')

df = pd.read_excel('/home/ylucas/Bureau/SP_strains_metadata.xlsx')
row = np.where(df.values == 'SP1')[0]
print(df['chloramphenicol'].values[row])

