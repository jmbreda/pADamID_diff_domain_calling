import numpy as np
import pyBigWig as bw
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
sys.path.insert(1, 'scripts')
from get_pADamID_datasets import *
from write_to_bigwig import *
from domain_calling import *
import os
import pandas as pd
from scipy.stats import beta, norm

# data
OUTFOLDER = 'results/significant_diff_regions'
bin_size = 20 # kb
Datasets = get_pADamID_datasets()
# more stringent
eps=0.001
percertiles = [eps,1-eps]
Win = 40 # window size (nr. of bins) to compute diff significance
min_region_size = 5*bin_size # kb
n_bin_bridge = 3
# less stringent
eps=0.01
percertiles = [eps,1-eps]
Win = 50 # window size (nr. of bins) to compute diff significance
min_region_size = 10*bin_size # kb
n_bin_bridge = 10


# loop over in files
for data in ['Top2B_2']:#Datasets.keys():
	print(data)

	# get tables
	print('\tget data')
	[Treat,Ctrl] = get_pADamID_tables(Datasets[data])

	# tables zscore normalization, and standard error rescaling
	Ctrl.value = (Ctrl.value - Ctrl.value.mean())/Ctrl.value.std()
	Ctrl.value_std_err = Ctrl.value_std_err/Ctrl.value.std()
	Treat.value = (Treat.value - Treat.value.mean())/Treat.value.std()
	Treat.value_std_err = Treat.value_std_err/Treat.value.std()

	# Compute differencial track
	Diff = Treat.copy(deep=True)
	# remove chrY
	Diff.drop(index=Diff[Diff.chr=='chrY'].index,inplace=True)
	Diff['value'] = Treat['value'] - Ctrl['value']

	# check if value have standard error
	std_err='value_std_err' in Diff.columns

	# Compute Diff standard error
	if std_err:
		Diff['value_std_err'] = np.sqrt( Treat['value_std_err']**2 + Ctrl['value_std_err']**2 )
		# Compute normal cdf from 0 to inf
		mu = Diff.value.values
		sig = Diff.value_std_err.values
		mu[np.isnan(mu)] = 0
		sig[np.isnan(sig)] = 1
		Diff['norm_cdf'] = 1-norm.cdf(0,mu,sig)
		outfolder = f'{OUTFOLDER}/{data}'
	else:
		outfolder = f'{OUTFOLDER}/{data}_no_std_err'
	if not os.path.exists(outfolder):
		os.mkdir(outfolder)

	# Write bw tracks
	bw_file = f'{outfolder}/{data}_ctrl.bw'
	write_bw(Ctrl,bw_file)
	bw_file = f'{outfolder}/{data}_treat.bw'
	write_bw(Treat,bw_file)
	bw_file = f'{outfolder}/D_{data}.bw'
	write_bw(Diff[['chr','start','end','value']],bw_file)

	# get first estimate of domains
	print('\tget domains')
	Diff = get_significant_regions_in_win(Diff,Win,percertiles,std_err)

	# collapse in domains
	Diff_domain = collapse_domains(Diff[['chr','start','end','domain']])

	if False:
		print('\tremove small domains')
		Diff_domain = remove_small_domains(Diff_domain,min_region_size,bin_size)

		# update domains in Diff
		for i in Diff_domain.index:
			idx_chr = np.array( Diff['chr']==Diff_domain.loc[i,'chr'] )
			idx_start = np.array( Diff['start']>=Diff_domain.loc[i,'start'] )
			idx_end = np.array( Diff['end']<=Diff_domain.loc[i,'end'] )
			Diff.loc[ idx_chr*idx_start*idx_end ,'domain'] = Diff_domain.loc[i,'domain']

	# refine regions
	print('\textend domains')
	# for each pos/neg region try to extend if it increases q_{.01} for pos or decreases q_{.99} for neg.
	extend_domains(Diff,Diff_domain,1,percertiles,std_err,n_bin_bridge)
	extend_domains(Diff,Diff_domain,-1,percertiles,std_err,n_bin_bridge)

	# remove 0-length regions
	idx_out = Diff_domain[Diff_domain.start==Diff_domain.end].index
	Diff_domain.drop(index=idx_out,inplace=True)

	# collapse
	Diff_domain = collapse_domains(Diff_domain)

	# remove small domains
	print('\tremove small domains')
	Diff_domain = remove_small_domains(Diff_domain,min_region_size,bin_size)

	# update domains in Diff
	for i in Diff_domain.index:
		idx_chr = np.array( Diff['chr']==Diff_domain.loc[i,'chr'] )
		idx_start = np.array( Diff['start']>=Diff_domain.loc[i,'start'] )
		idx_end = np.array( Diff['end']<=Diff_domain.loc[i,'end'] )
		Diff.loc[ idx_chr*idx_start*idx_end ,'domain'] = Diff_domain.loc[i,'domain']

	# write outputs
	print('\twrite outputs')
	bw_file = f'{outfolder}/D_sign_regions.bw'
	write_bw(Diff_domain,bw_file)
	csv_file = f'{outfolder}/D_sign_regions.csv'
	Diff_domain.to_csv(csv_file)
	csv_file = f'{outfolder}/D_norm_counts_and_domains.csv'
	Diff.to_csv(csv_file)

	# also write regions in bed format
	bed_file = f'{outfolder}/D_sign_regions'
	write_bed(Diff_domain,bed_file)
