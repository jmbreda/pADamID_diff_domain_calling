import numpy as np
import pandas as pd
from sklearn.metrics.cluster import normalized_mutual_info_score

DATASETS = []
#DATASETS.append(('SCII','Top2B',True))
#DATASETS.append(('SC','Top1',True))
# reproducibility control:


# get data
Diff = {}
for i in range(8,11):
	data = ('SCII_r'+str(i),'Top2B_r'+str(i))
	# load input file
	infolder = f'data/significant_diff_regions/{data[1]}-{data[0]}'
	infile = f'{infolder}/D_{data[1]}_{data[0]}_norm_counts_and_domains.csv'
	Diff[str(i)] = pd.read_csv(infile,index_col=0)

Corr_regions = []
NMI = []
for i in range(8,10):
	for j in range(i+1,11):
		Corr_regions.append( np.corrcoef(Diff[str(i)].domain,Diff[str(j)].domain)[0,1] )
		NMI.append( normalized_mutual_info_score(Diff[str(i)].domain,Diff[str(j)].domain) )



