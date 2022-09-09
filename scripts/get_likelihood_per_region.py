import numpy as np
import pandas as pd
from scipy.stats import beta
from write_to_bigwig import *

# get data
for i in range(9,11):
	data = ('SCII_r8','SCII_r'+str(i))
	print(data)

	# load input file
	datafolder = f'data/significant_diff_regions/{data[1]}-{data[0]}'
	infile = f'{datafolder}/D_{data[1]}_{data[0]}_norm_counts_and_domains.csv'
	Diff = pd.read_csv(infile,index_col=0)
	infile = f'{datafolder}/D_{data[1]}_{data[0]}_sign_regions_no_small_extended.csv'
	Diff_domain = pd.read_csv(infile,index_col=0)

	Diff_domain['likelihood'] = np.zeros(len(Diff_domain))
	for i in Diff_domain.index:
		idx_chr = np.array( Diff['chr']==Diff_domain.loc[i,'chr'] )
		idx_start = np.array( Diff['start']>=Diff_domain.loc[i,'start'] )
		idx_end = np.array( Diff['end']<=Diff_domain.loc[i,'end'] )
		n_nan = np.sum(np.isnan(Diff.loc[idx_chr*idx_start*idx_end ,'value']))
		p = np.sum(Diff.loc[idx_chr*idx_start*idx_end ,'value'] > 0) + n_nan/2
		n = np.sum(Diff.loc[idx_chr*idx_start*idx_end ,'value'] < 0) + n_nan/2
		if Diff_domain.loc[i,'domain'] == 1:
			Diff_domain.loc[i,'likelihood'] = 1 - beta.cdf(0.5,p+1,n+1)
		elif Diff_domain.loc[i,'domain'] == -1:
			Diff_domain.loc[i,'likelihood'] = - beta.cdf(0.5,p+1,n+1)

	bw_file = f'{datafolder}/D_{data[1]}_{data[0]}_domain_likelihood.bw'
	write_bw(Diff_domain[['chr','start','end','likelihood']],bw_file)


