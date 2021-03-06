    MAGITICS is a research project using machine learning for the
    prediction of antimicrobial resistance (AMR). More precisely, we
    leverage DNA sequence data in the form of k-mers to predict the
    resistance of the strains and potentially identify new resistance
    genes.


INSTALLATION

Python3 packages required

   .. code:: shell

       $ Pandas
       $ Scikit-learn
       $ Numpy
       $ Tensorflow

Clone

.. code:: shell

    $ git clone https://github.com/yvanlucas/magitics_project

Paths configuration


Before starting, you must specify the path to your datas in the
*config.py* file.

USAGE

Kmer parsing into sparse matrix (scipy.sparse.coo\_matrix):

::

    from data import Kmercount_to_matrix
    parser = Kmercount_to_matrix()
    parser.run() 

Learning & testing


You can specify before in *config.py* the model you want to use, or pass
one as argument

::

    model = 'gradient'  # can be ['rf','SCM', 'gradient', 'Ada']

Learning and testing all at once (RAM consuming process):

::

    from learning import Train_kmer_clf
    learner = Train_kmer_clf(clf= *model*)
    learner.run(evaluate=True)

Learning and testing in streaming mode (RAM efficient):

::

    from learning import Train_kmer_clf,  Test_streaming
    learner = Train_kmer_clf(clf= *model*)
    learner.run(evaluate=False)

    tester = Test_streaming(batchsize=1, kmer_to_index=learner.kmer_to_index, clf=learner.cv_clf)
    tester.run()


INPUT FORMAT

-  Fasta files
   (`Example <https://docs.google.com/document/d/1YcsyUQMT7aTPcqjS0t5VULW1UO9O5rHCiLeukZ47jjc/edit?usp=sharing>`__)

OUTPUT FORMAT

-  Report file
   (`Example <https://docs.google.com/document/d/1_gKsC3LW8TkwoGkbpJm-ubarPMmP7lEJoW3HCsvKzZo/edit?usp=sharing>`__)
-  Pickle file of the form:

::
    output = {"classifier": clf, "features": kmer_to_index, "y_pred": y_preds,
     "y_pruned": y_pruned, "y_true": y_test, "score":score}


CONTACT
`Dr. Yvan Lucas <mailto:yvanlucas44@gmail.com>`__
