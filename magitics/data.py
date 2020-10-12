import os
import pickle
import numpy as np
import pandas as pd
import scipy.sparse as sp
# from Bio import SeqIO
import config as cfg


class Kmer_parser(object):
    """
    Parse kmers from a fasta file using command-line softwares (DSK or Gerbil)
    DSK: https://github.com/GATB/dsk
    Gerbil: https://github.com/uni-halle/gerbil
    class called in KmersCount2Dataframe class
    """

    def __init__(self, fastaname):
        """
        Init class with the name of the file from which kmers are parsed

        Args:
            fastaname:
        """
        self.pathtofasta = os.path.join(cfg.pathtodata, cfg.data, fastaname)
        self.strainnumber = self.pathtofasta.split('/')[-1][:-3]
        self.label = fastaname[:5]
        self.available_commands = ["parse_kmers_dsk", "parse_kmers_gerbil"]
        self.len_kmers = cfg.len_kmers
        self.pathtotemp = os.path.join(cfg.pathtoxp, cfg.xp_name, "temp", self.label + self.strainnumber)
        self.pathtosave = os.path.join(cfg.pathtoxp, cfg.xp_name, "kmers", self.label + self.strainnumber)

    def parse_kmers_gerbil(self):
        """
        Use the gerbil software to create a kmer-count file from the fasta file containing contigs
        """
        kmerCmd = "gerbil -k %d   %s %s %s" % (
            self.len_kmers, self.pathtofasta, self.pathtotemp, self.pathtosavetemp)
        os.system(kmerCmd)
        tofastaCmd = "toFasta %s %d %s" % (self.pathtosavetemp, self.len_kmers, self.pathtosave)
        os.system(tofastaCmd)

    def count_kmers_gerbil(self):
        """
        Parse the kmer-count file generated by gerbil in a dictionnary of the form {kmer: count} for the considered strain
        """
        self.kmer_counts = {}
        with open(self.pathtosave, "r") as fasta:
            lines = fasta.readlines()
            self.kmer_counts["strain"] = self.strainnumber
            for kmercount, kmerID in zip(*[iter(lines)] * 2):
                try:
                    count = int(kmercount[1:])
                    ID = kmerID[:-1]
                    self.kmer_counts[ID] = count
                except Exception as e:
                    print(e)
                    print(kmercount, kmerID)

    def parse_kmers_dsk(self):
        """
        Use the DSK software to create a kmer-count file from the fasta file containing contigs
        """
        print(self.pathtofasta)
        kmerCmd = "dsk -file %s -out %s -kmer-size %d -abundance-min 1 -verbose 0" % (
            self.pathtofasta, self.pathtotemp, self.len_kmers)
        os.system(kmerCmd)
        outputCmd = "dsk2ascii -file %s -out  %s" % (self.pathtotemp, self.pathtosave)
        os.system(outputCmd)

    def count_kmers_dsk(self):
        """
        Parse the kmer-count file generated by DSK in a dictionnary of the form {kmer: count} for the considered strain
        """
        self.kmer_counts = {}
        with open(self.pathtosave, "r") as fasta:
            lines = fasta.readlines()
            for line in lines:
                try:
                    [ID, count] = line.split(" ")
                    self.kmer_counts[str(ID)] = int(count)
                except Exception as e:
                    print(e)
                    print("line = " + line)


class plfam_parser(object):
    """
    Parse the plfam list (list of genes identifiers) for a given strain using the PATRIC.features.tab sheet.
    """

    def __init__(self, strainID):
        self.strainID = strainID
        self.pathtofile = os.path.join(cfg.pathtodata, cfg.data, strainID + '.PATRIC.features.tab')

    def run(self):
        """
        Iterate through the lines of the PATRIC.features.tab file

        Returns:
            list of plfams in the file (some may appear twice since genes are sometimes redundant), y: prediction target
        """
        # plfams = []
        # with open(self.pathtofile, 'r') as f:
        #     lines = f.readlines()
        # for line in lines:
        #     try:
        #         plfams.append(line.split('\t')[16])
        #     except Exception as e:
        #         print(e)
        #         print(line.split('\t'))

        df=pd.read_csv(self.pathtofile, sep='\t')
        plfams=df['plfam_id'].to_list()
        y = self.strainID[:9]
        return plfams, y


class plfam_to_matrix(object):
    """
    Use the class plfam_parser to parse plfams and create a pandas dataframe of plfams for AMR prediction
    """

    def __init__(self):
        return

    def get_list_plfams(self, dic_data):
        """
        get the list of unique plfams for initiating the pandas dataframe

        Args:
            dic_data: {strain1:[plfam1, plfam2, ...]}

        Returns:
            list of unique plfams
        """

        ls_plfams = []
        for strain in dic_data:
            ls_plfams.extend(dic_data[strain])
        ls_plfams = list(set(ls_plfams))
        return ls_plfams

    def fill_dic(self, dicdata, ls_plfams):
        """
        Use dicdata and ls_plfams to create a dictionary (dic_df) from which a Pandas dataframe can be created


        Args:
            dicdata: dic {strainID: [plfams contained in the strain genome]
            ls_plfams: list of unique plfams

        Returns:
            dic_df: dictionary from which a Pandas dataframe can be created
        """
        dic_df = {}
        for strainID in dicdata:
            dic_df[strainID] = {}
            for plfam in ls_plfams:
                dic_df[strainID][plfam] = 0
            for plfam in dicdata[strainID]:
                dic_df[strainID][plfam] += 1
        return dic_df

    def run(self):
        """
        Run method to wrap the class methods and create the plfams Pandas dataframe
        """
        dic_data = {}
        targets = []
        for filename in os.listdir(os.path.join(cfg.pathtodata, cfg.data)):
            print(filename)
            strainID = filename.split('/')[-1][:-20]
            parser = plfam_parser(strainID)
            plfams, y = parser.run()
            targets.append(y)
            dic_data[strainID] = plfams
        ls_plfam = self.get_list_plfams(dic_data)
        dic_df = self.fill_dic(dic_data, ls_plfam)
        self.df = pd.DataFrame(dic_df).T
        # self.df['target']=targets
        with open(os.path.join(cfg.pathtoxp, cfg.xp_name, cfg.id, "kmers_DF.pkl"), "wb") as f:
            pickle.dump([self.df, targets], f, protocol=4)


class Kmercount_to_matrix(object):
    """
    This class allows us to iterate over fastas file, count kmers using the class KmerExtractionandCount and create a
    scipy sparse csr matrix or a pandas dataframe using these kmers counts
    """

    def __init__(self):
        #self.kmerdicts=0
        #if not self.kmerdicts:
        #    print("loading kmers dictionary")
        #    with open(os.path.join(cfg.pathtoxp, cfg.xp_name, "kmerdicts.pkl"), "rb") as f:
        #        self.kmerdicts = pickle.load(f)
        return


    def extend_kmerdicts(self, kmer):
        """
        Merge all individual kmer_count dictionaries created with kmers_parser into a dictionary of the form
        {kmer1:{strain_a:count, strain_b:count}, kmer2:{strain_a:count, strain_c:count...},...}
        Args:
            kmer: object of the class kmers_parser containing a dict {kmer1:count, kmer2:count,...}

        Returns:
            self
        """

        for key in kmer.kmer_counts.keys():
            if key in self.kmerdicts: #todo pass kmerdict as param instead
                self.kmerdicts[key][kmer.strainnumber] = int(kmer.kmer_counts[key])
            else:
                self.kmerdicts[key] = {kmer.strainnumber: int(kmer.kmer_counts[key])}

    def create_sparse_coos(self):
        """
        Iterate over the keys of the kmerdict containing the kmer_counts of all the strains in order to gather coordinates
        and datas informations for creating the sparse_csr matrix

        Returns:
            rows, cols, datas: coordinates and values of kmercount for the sparse matrix
        """
        print("*** Creating matrix ***")
        rows = []
        columns = []
        data = []
        for kmer in self.kmerdicts: #todo pass kmerdict as param instead
            for strain in self.kmerdicts[kmer]:
                rows.append(self.strain_to_index[strain])
                columns.append(self.kmer_to_index[kmer])
                if cfg.kmer_count == 1:
                    data.append(self.kmerdicts[kmer][strain])
                else:
                    data.append(1)
        del self.kmerdicts
        return rows, columns, data

    def populate_sparse_matrix(self, rows, cols, data):
        """
        Create a sparse_csr matrix with the coordinates contained in rows cols, data

        Args:
            rows: list of row coordinates
            cols: list of cols coordinates
            data: kmer count values for the coordinates

        Returns:
            self

        """
        n_strains = len(self.strains)
        self.mat = sp.csr_matrix((data, (rows, cols)), shape=(n_strains, len(self.kmer_to_index)), dtype=np.int8)

        mkdircmd = "mkdir %s" % (os.path.join(cfg.pathtoxp, cfg.xp_name, cfg.id))
        os.system(mkdircmd)

        with open(os.path.join(cfg.pathtoxp, cfg.xp_name, cfg.id, "kmers_mats.pkl"), "wb") as f:
            pickle.dump([self.mat, self.labels, self.strain_to_index, self.kmer_to_index], f, protocol=4)

    def create_dataframe(self):
        import pyarrow as pa
        import pyarrow.parquet as pq
        print("*** Creating dataframe ***")
        if not self.kmerdicts:
            with open(os.path.join(cfg.pathtoxp, cfg.xp_name, "kmerdicts.pkl"), "rb") as f:
                self.kmerdicts = pickle.load(f)
        self.kmerdb = pd.DataFrame(self.kmerdicts)
        table = pa.Table.from_pandas(self.kmerdicts.transpose)
        pq.write_table(table, os.path.join(cfg.pathtoxp, cfg.xp_name, "kmers_DF.parquet"))

    def clean_temp_directories(self, kmer):
        """
        delete the intermediate file that were created during the command line call of the kmer-parser

        Args:
            kmer: kmer object from which the path are extracted

        Returns:
            self
        """

        cleankmertempcmd = "rm -rf %s" % (kmer.pathtotemp)
        os.system(cleankmertempcmd)
        # cleantempcmd="rm -rf %s" % (self.kmer.pathtosavetemp)
        # os.system(cleantempcmd)


    def get_label_from_csv_metadata(self, strain):
        """
        Get label from csv metadata
        TODO: pass parameters for col, metadata_path

        Args:
            strain: strainID

        Returns:
            label
        """
        import pandas as pd
        # df = pd.read_excel('/home/ylucas/Bureau/SP_strains_metadata.xlsx')
        df = pd.read_excel('/scratch/MAGITICS_data/Streptococcus_pneumoniae/SP_strains_metadata.xlsx')
        row=np.where(df.values==strain)[0]
        print('label')
        print(int(df['chloramphenicol'].values[row]))
        return int(df['chloramphenicol'].values[row])

    def run(self):
        """
        Run method to wrap the class methods and create the kmer count sparse matrix
        """
        self.kmerdicts = {}
        self.labels = []
        self.strains = []
        for filename in os.listdir(os.path.join(cfg.pathtodata, cfg.data)):
            kmer = Kmer_parser(os.path.join(filename))
            kmer.parse_kmers_dsk()
            kmer.count_kmers_dsk()
            self.strains.append(kmer.strainnumber)
            #self.labels.append(kmer.label)
            self.labels.append(self.get_label_from_csv_metadata(kmer.strainnumber))
            print(self.labels)
            self.extend_kmerdicts(kmer)
            # self.clean_temp_directories(kmer)
        self.strain_to_index = {strain: i for i, strain in zip(range(len(self.strains)), self.strains)}
        self.kmer_to_index = {kmer: i for i, kmer in enumerate(self.kmerdicts)}
        rows, cols, data = self.create_sparse_coos()
        self.populate_sparse_matrix(rows, cols, data)


class kmer_gene_correspondance(object):
    """
    This class allows us to pairwise iterate over fastas and PATRIC.features.tab files, parse the genes sequences contained
    in the DNA sequence and associate the kmers with the genes they appear in in a dictionary.
    """

    def __init__(self):
        return

    def get_genes_limit(self, pathtofile):
        """
        Create a dictionary of the form {'plfam':[(limit1, limit2), numcontig],..}


        Args:
            pathtofile: path to the .PATRIC.features.tab

        Returns:
            dic_limits: {'plfam':[(limit1, limit2), numcontig],..}

        """
        with open(pathtofile, 'r') as f:
            lines = f.readlines()[1:]
        count_contigs = 0
        dic_limits = {}
        for line in lines:
            line = line.split('\t')
            if line[9] == '1':
                count_contigs += 1
            else:
                pgfam = line[16]
                limits = [(int(line[9]), int(line[10])), count_contigs]
                if pgfam in dic_limits:
                    dic_limits[pgfam].append(limits)
                else:
                    dic_limits[pgfam] = [limits]
        return dic_limits

    def extract_sequence(self, pathtofile):
        """
        return a list of contigs DNA sequences from the fasta file

        Args:
            pathtofile: path to the fasta file

        Returns:
            contigs: list of assembled DNA sequence (contigs)
        """
        with open(pathtofile, 'r') as f:
            file = f.readlines()
        contigs = ['']
        skipfirst = 1
        for line in file:
            if skipfirst == 0:
                contigs[-1] = contigs[-1] + str(line)[:-1]
            skipfirst = 0
            if len(line) == 1:
                contigs.append('')
                skipfirst = 1
        return contigs

    def extract_genes_from_seq(self, contigs, dic_limits):
        """
        Use the genes limits in dic_limits to extract genes sequences from the fasta

        Args:
            contigs: list of assembled DNA sequence (contigs)
            dic_limits: {'plfam':[(limit1, limit2), numcontig],..}

        Returns:
            dic_genes: {'plfam':[geneseq1, geneseq2],..}
        """

        dic_genes = {}
        for pgfam in dic_limits:
            dic_genes[pgfam] = []
            for limits in dic_limits[pgfam]:
                geneseq = contigs[limits[1] - 1][limits[0][0]:limits[0][1]]
                dic_genes[pgfam].append(geneseq)
        return dic_genes

    def kmers_within_gene(self, geneseq, len_kmers=31):
        """
        extract kmers from geneseq


        Args:
            geneseq: DNA sequence
            len_kmers: len of the kmer extracted

        Returns:
            ls_kmers: liste of kmers contained in geneseq
        """
        ls_kmers = []
        for i in range(len(geneseq) - len_kmers):
            ls_kmers.append(geneseq[i:i + len_kmers])
        return ls_kmers

    def build_kmer_gene_dict(self, dic_kmer_gene, dic_gene):
        """
        For each geneseq in dic_genes, parse kmers and fill dic_kmer_gene


        Args:
            dic_kmer_gene: {kmer:[pgfam1, pgfam2, ..]}
            dic_gene: {pgfam:[seq1, seq2, ..]}

        Returns:
            filled dic_kmer_gene: {kmer:[pgfam1, pgfam2, ..]}
        """

        for pgfam in dic_gene.keys():
            for seq in dic_gene[pgfam]:
                ls_kmers = self.kmers_within_gene(seq)
                for kmers in ls_kmers:
                    if kmers in dic_kmer_gene:
                        if pgfam not in dic_kmer_gene[kmers]:
                            dic_kmer_gene[kmers].append(pgfam)
                    else:
                        dic_kmer_gene[kmers] = [pgfam]
        return dic_kmer_gene

    def get_genes_from_kmers(self):
        """
        Method used to understand results
        """
        return

    def run(self):
        """
        Run method to wrap the class methods and create the dic_kmer_gene
        dic_kmer_gene: {kmer:[pgfam1, pgfam2, ..]}

        Returns:
            self
        """
        dic_kmer_gene = {}
        ls_file = os.listdir(os.path.join(cfg.pathtodata, cfg.data))
        ls_file.sort()
        for filename in ls_file[::2]:
            print(filename)
            limitfile = os.path.join(cfg.pathtodata, cfg.data, filename)
            seqfile = os.path.join(cfg.pathtodata, cfg.data, filename[:-20] + '.fna')
            dic_limits = self.get_genes_limit(limitfile)
            contigs = self.extract_sequence(seqfile)
            dic_genes = self.extract_genes_from_seq(contigs, dic_limits)
            dic_kmer_gene = self.build_kmer_gene_dict(dic_kmer_gene, dic_genes)
        with open(os.path.join(cfg.pathtoxp, cfg.xp_name, cfg.id, 'dic_kmer_gene.pkl'), 'wb') as f:
            pickle.dump(dic_kmer_gene, f)

# genedf=plfam_to_matrix()

# genedf.run()
#
# dic_kmer=kmer_gene_correspondance()
# dic_kmer.run()
#
#
# for kmer in list(dic_kmer.dic_kmer_gene.keys())[:25]:
#     print(dic_kmer.dic_kmer_gene[kmer])
#     print(kmer)
