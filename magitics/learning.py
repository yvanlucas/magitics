import os
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns
from sklearn import (ensemble, tree, feature_selection, metrics, model_selection, preprocessing)
import config as cfg


class Train_kmer_clf(object):
    """
    Train in batch version and optionnally test in batch settings
    Also optimize the treshold for computing accuracies
    """
    def __init__(self):
        if cfg.model == "rf":
            self.clf = ensemble.RandomForestClassifier()
            self.param_grid = cfg.rf_grid
        elif cfg.model == "gradient":
            #self.clf = ensemble.GradientBoostingClassifier(max_depth=4, max_features=None)
            self.clf = ensemble.GradientBoostingRegressor(max_features=None)
            self.param_grid = cfg.gradient_grid
        elif cfg.model == 'Ada':
            self.clf = ensemble.AdaBoostClassifier()
            self.param_grid = cfg.ada_grid

        if cfg.dtype == 'sparse':
            with open(os.path.join(cfg.pathtoxp, cfg.xp_name, cfg.id, "kmers_mats.pkl"), "rb") as f:
                [self.mat, self.labels, self.strain_to_index, self.kmer_to_index] = pickle.load(f)
            self.testratio = 0.0
        elif cfg.dtype == 'df':
            with open(os.path.join(cfg.pathtoxp, cfg.xp_name, cfg.id, "kmers_DF.pkl"), "rb") as f:
                [self.mat, self.labels] = pickle.load(f)
            self.kmer_to_index = self.mat.columns
            self.testratio = 0.3

    def preprocess_y(self):
        """
        Transform y into a binary vector
        """
        # to_drop = ["label", "strain"]
        # X = self.mat.drop(to_drop, axis=1)
        # self.columns = X.columns
        if cfg.MIC==False:
            le = preprocessing.LabelEncoder()
            self.y = le.fit_transform(self.labels)
        else:
            print(self.labels)
            self.y=self.labels


    def split_train_test(self):
        """
        Split the data matrix and target vector into train and test matrices and vector

        Returns:
            X_train, X_test, y_train, y_test
        """
        if self.testratio > 0:
            X_train, X_test, y_train, y_test = model_selection.train_test_split(self.mat, self.y, test_size=self.testratio)
        else:
            X_train = self.mat
            y_train = self.y
            X_test = None
            y_test = None
        print(y_train)
        del self.mat, self.y
        return X_train, X_test, y_train, y_test

    def chi2_feature_selection(self, X_train, X_test, y_train):
        """
        Refactor X_train and X_test only keeping features that are correlated with the target y_train

        Args:
            X_train: train numpy array
            X_test: test numpy array
            y_train: train target variable vector

        Returns:
            refined X_train and X_test
        """

        chi2_selector = feature_selection.SelectKBest(feature_selection.chi2, k=1000000)
        X_train = chi2_selector.fit_transform(X_train, y_train)
        X_test = chi2_selector.transform(X_test)
        return X_train, X_test

    def fit(self, X_train, y_train):
        """
        Fit the chosen classifier using gridsearch cross validation

        Args:
            X_train: input data matrix
            y_train: ground truth vector

        Returns:

        """
        self.cv_clf = model_selection.GridSearchCV(estimator=self.clf, param_grid=self.param_grid, cv=2,
                                                   scoring="roc_auc", n_jobs=1)
        self.cv_clf.fit(X_train, y_train)
        self.y_pred = self.cv_clf.predict_proba(X_train)
        with open(os.path.join(cfg.pathtoxp, cfg.xp_name, cfg.id, f'{cfg.model}_CVresults.pkl'), 'wb') as f:
            pickle.dump({"classifier": self.cv_clf, "features": self.kmer_to_index}, f, protocol=4)

    def predict(self, X_test):
        """
        Batch predict of X_test labels

        Args:
            X_test: input data matrix

        Returns:
            y_predict: prediction vector
        """
        y_predict = self.cv_clf.predict_proba(X_test)
        return y_predict



    def get_accuracy_treshold(self, X_train, y_train):
        """
        Calculate the treshold to obtain the best accuracy on the train set

        Args:
            X_train: input data matrix
            y_train: ground truth vector

        Returns:
            treshold value
        """

        train_predict = self.cv_clf.best_estimator_.predict_proba(X_train)[:, -1]
        accuracies = []
        nsteps = 100
        for i in range(nsteps):
            tres = i / nsteps
            tresd_predict = []
            for pred in train_predict:
                if pred > tres:
                    tresd_predict.append(1)
                else:
                    tresd_predict.append(0)
            accuracies.append(metrics.accuracy_score(y_train, tresd_predict))
        ind = accuracies.index(max(accuracies))
        treshold = float(ind) / float(nsteps)
        with open(os.path.join(cfg.pathtoxp, cfg.xp_name, cfg.id, f"{cfg.model}_tres_value.txt"), "w") as f:
            f.write(str(treshold))
        return treshold

    def adapted_accuracy(self, y_test, y_pred, treshold):
        """
        Predict accuracy with respect to a calculated treshold (inherited from Test_streaming class)

        Args:
            y_test: ground truth vector
            y_pred: prediction vector
            treshold

        Returns:
            adapted accuracy score
        """
        return Test_streaming.adapted_accuracy(self, y_test, y_pred, treshold)

    def evaluate_and_write_report(self, y_pred, y_test, treshold):
        """
        Predict scores and write report

        Args:
            y_test: ground truth vector
            y_pred: prediction vector
            treshold
        """
        Test_streaming.evaluate_and_write_report(self, y_pred, y_test, treshold)

    def plot_CV_heatmap(self):
        """
        Plot gridsearchCV heatmap
        """
        ls_params = list(self.param_grid.keys())
        self.pvt = pd.pivot_table(pd.DataFrame(self.cv_clf.cv_results_), values="mean_test_score",
                                  index="param_" + ls_params[0], columns="param_" + ls_params[1])
        ax = sns.heatmap(self.pvt)
        ax.set(ylabel=ls_params[0], xlabel=ls_params[1])
        ax.figure.savefig(os.path.join(cfg.pathtoxp, cfg.xp_name, cfg.id, f"{cfg.model}_gridCV_heatmap.png"))

    def plot_boosting_learning_curve(self, X_test, y_test):
        """
        Plot boosting learning curve (test/train deviance = f (n_estimator) )

        Args:
            X_test: input data matrix
            y_test: ground truth vector
        """
        if cfg.model == "gradient":
            test_score = np.zeros((self.cv_clf.best_params_["n_estimators"],), dtype=np.float64)
            for i, y_pred in enumerate(self.cv_clf.best_estimator_.staged_predict(X_test)):
                test_score[i] = self.cv_clf.best_estimator_.loss_(y_test, y_pred)
            sns.set()
            fig = plt.figure(figsize=(6, 6))
            plt.subplot(1, 1, 1)
            plt.title(cfg.xp_name + "  //  ROC AUC = " + str(self.score))
            plt.plot(np.arange(self.cv_clf.best_params_["n_estimators"]) + 1, self.cv_clf.best_estimator_.train_score_,
                     "b-", label="Training Set Deviance")
            plt.plot(np.arange(self.cv_clf.best_params_["n_estimators"]) + 1, test_score, "r-",
                     label="Test Set Deviance")
            plt.legend(loc="upper right")
            plt.xlabel("Boosting Iterations")
            plt.ylabel("Deviance")
            fig.tight_layout()
            plt.savefig(os.path.join(cfg.pathtoxp, cfg.xp_name, cfg.id, f"{cfg.model}boosting_learning_curve.png"))

    def run(self):
        """
        Run method to wrap the class methods, train and optionnaly test the model in batch settings
        """
        self.preprocess_y()
        X_train, X_test, y_train, y_test = self.split_train_test()
        # X_train, X_test = self.chi2_feature_selection(X_train, X_test, y_train)
        self.fit(X_train, y_train)
        tres = self.get_accuracy_treshold(X_train, y_train)
        if cfg.dtype == 'df':
            y_predict = self.predict(X_test)
            self.evaluate_and_write_report(y_predict, y_test, tres)
            self.plot_CV_heatmap()
            self.plot_boosting_learning_curve(X_test, y_test)

            with open(os.path.join(cfg.pathtoxp, cfg.xp_name, cfg.id, f"{cfg.model}_CVresults.pkl"), "wb") as f:
                pickle.dump(
                    {"classifier": self.cv_clf, "features": self.kmer_to_index, "y_pred": y_predict, "y_true": y_test},
                    f,
                    protocol=4)


class NN_model(object):
    """
    Work in progress: neural network for end classification, because it's nice, and also in order to wrap different input datas
    of different sources in an architecture.
    """
    def __init__(self):
        return

    def create_model(self, X_train):
        import tensorflow as tf

        # Complete version
        model = tf.keras.Sequential([
            tf.keras.layers.InputLayer(input_shape=(X_train.shape[1],)),  # Input Layer
            tf.keras.layers.Dense(64, activation='relu'),  # Hidden Layer
            tf.keras.layers.Dense(1)])  # Output

        # Optional (easy)
        # Here the InputLayer will be generated automatically
        # model = tf.keras.Sequential([
        #     tf.keras.layers.Dense(64, activation='relu', input_shape=(X_train.shape[1],)), # Hidden Layer
        #     tf.keras.layers.Dense(1)]) # Output

        # Compiling our model
        model.compile(loss='mse',
                      optimizer=tf.keras.optimizers.RMSprop(),
                      metrics=['mae'])

        # Summary
        model.summary()
        return model

    def fit(self, model, X_train, y_train):
        model.fit(X_train, y_train, epochs=110, verbose=2)

        loss, mae = model.evaluate(X_val, y_val, verbose=2)
        y_pred = model.predict(X_test).flatten() #what is flatten
        return model


class Test_streaming(object):
    """
    Test in stream settings
    Also optionally prune the classifier discarding trees that use
    """
    def __init__(self, kmer_to_index=None, clf=None, batchsize=10):
        self.batchsize = 1

        self.testdir = os.path.join(cfg.pathtodata, cfg.testdir)
        self.kmer_to_index = kmer_to_index
        try:
            self.clf = clf.best_estimator_
        except Exception as e:
            print(e)
            self.clf = clf
        self.pathtotemp = os.path.join(cfg.pathtoxp, cfg.xp_name, cfg.id, "test-temp")
        self.pathtosave = os.path.join(cfg.pathtoxp, cfg.xp_name, cfg.id, "test-output")

        if not (os.path.isdir(self.pathtotemp) and os.path.isdir(self.pathtosave)):
            mkdirCmd = "mkdir %s" % (self.pathtotemp)
            os.system(mkdirCmd)
            mkdirCmd = "mkdir %s" % (self.pathtosave)
            os.system(mkdirCmd)

    def create_sparse_coos(self, cols, rows, datas, col, row, data):
        """
        Extend existing list of coordinates with new coordinates corresponding to the kmers counts of the current file
        """
        cols.extend(col)
        rows.extend(row)
        datas.extend(data)
        return cols, rows, datas

    def prune_boosting(self):
        """
        Prune boosting according to the kmer redundancy:
        if two very redundant kmers are used for the prediction:
        prune one of the two trees of depth=1 with the least importance
        """
        import difflib
        # Select index of redundant kmers with lower importances
        ls_index = []
        featimp = self.clf.feature_importances_
        kmers = [list(self.kmer_to_index.keys())[i] for i in np.nonzero(featimp)[0]]
        imps = [featimp[i] for i in np.nonzero(featimp)[0]]
        index = [i for i in np.nonzero(featimp)[0]]
        for kmer1, imp1, ind1 in zip(kmers, imps, index):
            for kmer2, imp2, ind2 in zip(kmers, imps, index):
                similarity = difflib.SequenceMatcher(None, kmer1, kmer2).ratio()
                if similarity > cfg.pruning_tresh and kmer1 != kmer2:
                    if imp1 > imp2:
                        ls_index.append(ind2)
                    elif imp2 > imp1:
                        ls_index.append(ind1)
        return list(set(ls_index))  # list of redundant kmer indexes

    def predict_pruned(self, X_test, ls_index):
        """
        Predict using only trees retained by the method prune boosting

        Args:
            X_test: input data matrix
            ls_index: list of trees index to use for pruned prediction

        Returns:
            prediction vector
        """

        cumpred = np.array([x for x in self.clf.staged_decision_function(X_test)])[:, :, 0]
        preds_out = cumpred[-1, :]
        for i in ls_index:  # i can't be 0 but who would prune first tree of boosting
            preds_out = preds_out - (cumpred[i - 1, :] - cumpred[i, :])
        return preds_out

    def populate_sparse_matrix(self, cols, rows, datas, batchsize=1):
        """
        Create a sparse matrix using coordinate vectors

        Args:
            cols: coordinate vector
            rows: coordinate vector
            datas: data vector
            batchsize
        """
        X_test = sp.csr_matrix((datas, (rows, cols)), shape=(batchsize, len(self.kmer_to_index)), dtype=np.int16)
        return X_test

    def append_prediction(self, X_test, y_preds, y_pruned, y_test, ls_index, y):
        """
        Predict and append predictions to prediction vectors

        Args:
            X_test: input data matrix
            y_preds: prediction list
            y_pruned: pruned prediction list
            y_test: ground truth list
            ls_index: list of trees index to use for pruned prediction
            y: labels of the current batch

        Returns:
            y_preds, y_pruned, y_test

        """
        y_preds.extend(self.clf.predict_proba(X_test))
        # y_pruned.extend(self.predict_pruned(X_test, ls_index))
        y_test.append(y)
        return y_preds, y_pruned, y_test

    def parse_and_map_kmers(self, fastaname, batchnumber):
        """
        Parse kmers_count using dsk and "get_kmer_counts" method
        and map them to the correct coordinates using the "map_data_to_coords" method

        Args:
            fastaname: path to the fasta file
            batchnumber: iteration counter of the current batch

        Returns:
            coordinates lists, label
        """
        self.parse_kmers_dsk(fastaname)
        y = fastaname[:5]
        y = self.get_label_from_csv_metadata(fastaname[:-3])
        kmer_count = self.get_kmer_counts(fastaname)
        col, row, data = self.map_data_to_coords(kmer_count, batchnumber)
        return col, row, data, y

    def get_label_from_csv_metadata(self, strain):
        """
	Get label from csv metadata
        TODO: pass parameters for col, metadata_path

        Args:
            strain: strainID

        Returns:
            label
        """
        # df = pd.read_excel('/home/ylucas/Bureau/SP_strains_metadata.xlsx')
        df = pd.read_excel('/scratch/MAGITICS_data/Streptococcus_pneumoniae/SP_strains_metadata.xlsx')
        row=np.where(df.values==strain)[0]
        return df['chloramphenicol'].values[row]


    def get_kmer_counts(self, kmercountfile):
        """
        Parse kmer_counts from the kmercountfile generated using dsk

        Args:
            kmercountfile: pathtokmercountfile

        Returns:
            kmer_count dict: {kmer1:count, kmer2: count, ...}

        """
        kmer_count = {}
        with open(os.path.join(self.pathtosave, kmercountfile), "r") as fasta:
            lines = fasta.readlines()
            for line in lines:
                try:
                    [ID, count] = line.split(" ")
                    kmer_count[str(ID)] = int(count)
                except:
                    print('line = {0}'.format(line))
        return kmer_count

    def map_data_to_coords(self, kmer_count, batchnumber):
        """
        Iterate over kmer_count dict to fill coordinates vector

        Args:
            kmer_count: dict {kmer1:count, kmer2: count, ...}
            batchnumber: iteration counter of the current batch

        Returns:
            coordinate lists
        """
        row = []
        data = []
        column = []
        for kmer in kmer_count:
            try:
                column.append(self.kmer_to_index[kmer])
                row.append(batchnumber)
                if cfg.kmer_count == 1:
                    data.append(kmer_count[kmer])
                else:
                    data.append(1)
            except:
                self.missing_kmers.append(kmer)
        return column, row, data

    def parse_kmers_dsk(self, fastaname):
        """
        Create a kmercount file using the command line software dsk

        Args:
            fastaname: path to the fasta file
        """
        kmerCmd = "dsk -file %s -out %s -kmer-size %d -abundance-min 1 -verbose 0" % (
            os.path.join(self.testdir, fastaname), os.path.join(self.pathtotemp, fastaname), cfg.len_kmers)
        os.system(kmerCmd)
        outputCmd = "dsk2ascii -file %s -out  %s" % (
            os.path.join(self.pathtotemp, fastaname), os.path.join(self.pathtosave, fastaname))
        os.system(outputCmd)

    def adapted_accuracy(self, y_test, y_preds, tres):
        """
        Predict accuracy score with respect to a pre-calculated treshold

        Args:
            y_test: ground truth vector
            y_preds: prediction vector
            tres: treshold

        Returns:
            adapted accuracy score
        """

        y_preds_adapted = []
        for pred in y_preds[:, -1]:
            if pred > float(tres):
                y_preds_adapted.append(1.0)
            else:
                y_preds_adapted.append(0.0)
        score = metrics.accuracy_score(y_test, y_preds_adapted)
        print(score)
        return score

    def evaluate_and_write_report(self, y_preds, y_test, tres=None, pruned=False):
        """
        Compute several evaluation metrics and save them in a score dictionary,
        then write report containing prediction parameters and scores.
        Args:
            y_preds: prediction vector
            y_test: ground truth vector
            tres: adated accuracy treshold
            pruned:  indicator if we are in a pruned setting or not
        """

        score = {}
        print(np.shape(y_test))
        print(np.shape(y_preds))
        score["ROC_AUC"] = metrics.roc_auc_score(y_test, y_preds[:, -1])
        if tres==None:
            with open(os.path.join(cfg.pathtoxp, cfg.xp_name, cfg.id, f"{cfg.model}_tres_value.txt"), "r") as f:
                tres = f.readlines()[0]
        score["Accuracy"] = self.adapted_accuracy(y_test, y_preds, tres)
        score["Accuracy"] = metrics.accuracy_score(y_test, y_preds[:, -1].round())
        score["MAE"] = metrics.mean_absolute_error(y_test, y_preds[:, -1])
        score["MSE"] = metrics.mean_squared_error(y_test, y_preds[:, -1])
        score["MAPE"] = metrics.mean_absolute_percentage_error(y_test, y_preds[:,-1])
        print("*** ROC AUC = ***")
        print(score["ROC_AUC"])
        self.write_report(pruned, score)
        return score

    def write_report(self,score, pruned=False):
        """
        Write report containing prediction parameters and scores

        Args:
            score: dict containing evaluation metrics values
            pruned: indicator if we are in a pruned setting or not
        """
        with open(os.path.join(cfg.pathtoxp, cfg.xp_name, cfg.id, f"{cfg.model}_report.txt"), "a") as txt:
            if pruned:
                txt.write('PRUNED' + "\n\n")
            txt.write(cfg.xp_name + "/" + cfg.id + "\n\n")
            txt.write(str(score) + "\n")
            txt.write("Len_kmers = " + str(cfg.len_kmers) + "\n")
            txt.write("Model = " + str(self.clf) + "\n")
            # txt.write("Best_params = "+str(self.clf.best_params_)+"\n")
            # txt.write("Param_grid = " + str(self.param_grid) + "\n")
            # txt.write("best params = " + str(self.clf.best_params_)+'\n')
            txt.write("\n Relevant kmers : \n")
            if cfg.model == "rf" or cfg.model == "gradient":
                featimp = self.clf.feature_importances_
                kmers = [list(self.kmer_to_index.keys())[i] for i in np.nonzero(featimp)[0]]
                for kmer in kmers:
                    txt.write(str(kmer) + "\n")

    def clean_temp_directories(self):
        """
        Delete the intermediate file that were created during the command line call of the kmer-parser
        """
        cleankmertempcmd = "rm -rf %s" % (self.pathtotemp)
        os.system(cleankmertempcmd)
        cleantempcmd = "rm -rf %s" % (self.pathtosave)
        os.system(cleantempcmd)

    def run(self):
        """
        Run method to wrap the class methods and test in stream settings
        """
        self.missing_kmers = []
        files = [file for file in os.listdir(self.testdir)]
        remaining = len(files)
        fileindex = 0
        y_test = []
        y_preds = []
        y_pruned = []
        ls_index = self.prune_boosting()
        while remaining > 0:
            batchiter = 0
            batch = min(remaining, self.batchsize)
            for file in files[fileindex: fileindex + batch]:
                cols = []
                rows = []
                datas = []
                col, row, data, y = self.parse_and_map_kmers(file, batchiter)
                cols, rows, datas = self.create_sparse_coos(cols, rows, datas, col, row, data)
                y = file[:5]
                batchiter += 1
                remaining -= 1
                X_test = self.populate_sparse_matrix(cols, rows, datas, batchiter)
                try:
                    y_preds, y_pruned, y_test = self.append_prediction(X_test, y_preds, y_pruned, y_test, ls_index, y)
                except Exception as e:
                    print('exception')
                    print(e)
            fileindex += batch
        y_preds = np.vstack(y_preds)
        le = preprocessing.LabelEncoder()
        y_test = le.fit_transform(y_test)
        score=self.evaluate_and_write_report(y_preds, y_test)
        print(ls_index)
        with open(os.path.join(cfg.pathtoxp, cfg.xp_name, cfg.id, f"{cfg.model}_CVresults.pkl"), "wb") as f:
            pickle.dump({"classifier": self.clf, "features": self.kmer_to_index, "y_pred": y_preds, "y_pruned": y_pruned,
                 "y_true": y_test, "score": score}, f, protocol=4)
        self.clean_temp_directories()


# train=Train_kmer_clf()
# train.run()
#
#
# #with open(os.path.join(cfg.pathtoxp, cfg.xp_name, cfg.id, f'{cfg.model}_CVresults.pkl'), 'rb') as f:
# #    dic=pickle.load(f)
# #test=Test_streaming(batchsize=1, kmer_to_index=dic['features'], clf=dic['classifier'])
# test = Test_streaming(batchsize=1, kmer_to_index=train.kmer_to_index, clf=train.cv_clf)
# test.run()
