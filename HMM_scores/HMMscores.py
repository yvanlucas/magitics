#TODO:
# For each strain:
# 1) Gather PLFAM ids and corresponding nuc sequence
# 2) HMMscan each locus against the corresponding HMM flatfile
# 3) get bitscore (susceptible - resistant) or only the one existings
# 4) put bitscore into a dic for each strain: strain: {plfam1:score1, ...}
# 5) Create dataframe