import numpy as np
import pandas as pd
from scipy.stats import beta

# Get beta posterior distribution of p in windows of 'win'
def get_significant_regions_in_win(Diff,Win,percentiles,std_err):
	win = int(Win/2)
	CHR = np.unique( Diff['chr'] )
	Diff['domain'] = np.zeros(len(Diff))

	for c in CHR:
		idx_c = Diff[Diff.chr==c].index
		for i in idx_c:
			i_low = max([idx_c[0],i-win])
			i_high = min([idx_c[-1],i+win])
			if std_err:
				p = np.sum(Diff.loc[i_low:i_high,'norm_cdf'])
				n = np.sum(1-Diff.loc[i_low:i_high,'norm_cdf'])
			else:
				n_nan = np.sum(np.isnan(Diff.loc[i_low:i_high,'value']))
				p = np.sum(Diff.loc[i_low:i_high,'value'] > 0) + n_nan/2
				n = np.sum(Diff.loc[i_low:i_high,'value'] < 0) + n_nan/2

			q = beta.ppf(percentiles,p+1,n+1)
			if q[0] > .5:
				Diff.loc[i,'domain'] = 1
			elif q[1] < .5:
				Diff.loc[i,'domain'] = -1
			#else: #if q[0] <= .5 and .5 <= q[1]:
			#	Diff_regions.iloc[i,3] = 0
	return Diff

def remove_small_domains(Diff,min_region_size,bin_size):

	CHR = np.unique(Diff['chr'])
	for c in CHR:
		# get first and last index of the chr
		[idx_first,idx_last] = Diff[Diff['chr']==c].index[[0,-1]]
		# loop by domains to remove (+,-,neutral)
		for dom in [1,-1,0]:
			# loop my size
			for min_size in range(bin_size*1000,min_region_size*1000,bin_size*1000):
				idx = Diff[(Diff['chr']==c)&(Diff.end-Diff.start==min_size)&(Diff.domain==dom)].index
				for i in idx:
					# ignore windows at chr borders
					if (i != idx_first) & (i != idx_last):
						same_chr = Diff.loc[i-1,'chr'] == Diff.loc[i+1,'chr']
						same_value = Diff.loc[i-1,'domain'] == Diff.loc[i+1,'domain']
						if all([same_chr,same_value]):
							Diff.loc[i,'domain'] = Diff.loc[i+1,'domain']
	Diff = collapse_domains(Diff)

	return Diff

def collapse_domains(df):

	out = []
	last_chr = None
	last_value = np.nan
	for l in df.values:

		# ignore nans
		if np.isnan(l[3]):
			continue

		# if new chromosome write if not nan and reinitialise start/end/values
		if l[0] != last_chr:
			if not np.isnan(last_value):
				out.append([last_chr,start,end,last_value])

			last_chr = l[0]
			start = l[1]
			end = l[2]
			last_value = l[3]
			continue

		if l[3] == last_value:
			# if value same as last, extend end
			end = l[2]
			last_value = l[3]
		# if new value and not nan, write and update start
		else:
			if not np.isnan(last_value):
				out.append([last_chr,start,end,last_value])

			last_chr = l[0]
			start = l[1]
			end = l[2]
			last_value = l[3]
	# write
	if not np.isnan(last_value):
		out.append([last_chr,start,end,last_value])

	out = pd.DataFrame(out)
	out.columns = df.columns

	return out


# refine regions
# for each pos/neg region try to extend if it increases q_{.01} for pos or decreases q_{.99} for neg.
def extend_domains(Diff,Diff_domain,domain_sign,percentiles,std_err,n_bin_bridge):

	bin_size = Diff.loc[0,'end'] - Diff.loc[0,'start']

	if not (domain_sign == 1 or domain_sign == -1):
		return

	if domain_sign == 1:
		percentile = percentiles[0]
	else:
		percentile = percentiles[1]

	# loop over domains
	for domain_idx in Diff_domain[Diff_domain['domain']==domain_sign].index:
		d = Diff_domain.loc[domain_idx,:]
		d_new = d.copy(deep=True)

		# get Diff indices of my region
		idx_chr = np.array( Diff['chr']==d.chr )
		idx_start = np.array( Diff['start']>=d.start )
		idx_end = np.array( Diff['end']<=d.end )
		my_idx = Diff.loc[ idx_chr*idx_start*idx_end ].index

		# skip empty idx
		if len(my_idx)==0:
			continue

		# compute cdf at percentile
		if std_err:
			p = np.sum(Diff.loc[my_idx,'norm_cdf'])
			n = np.sum(1-Diff.loc[my_idx,'norm_cdf'])
		else:
			n_nan = np.sum(np.isnan(Diff.loc[my_idx,'value']))
			p = np.sum(Diff.loc[my_idx,'value'] > 0) + n_nan/2
			n = np.sum(Diff.loc[my_idx,'value'] < 0) + n_nan/2
		q0 = beta.ppf(percentile,p+1,n+1)

		# extend
		for direction in ['l','r']:
			counter = 0
			n_bin = 0
			idx = my_idx
			while True:
				if direction == 'l':
					# avoid going out of dataframe
					if idx[0]-1 == -1:
						break
					# avoid changing chromosome
					if Diff.loc[idx[0]-1,'chr'] != Diff.loc[idx[0],'chr']:
						break
					# only extend in neutral regions
					if Diff.loc[idx[0]-1,'domain'] != 0:
						break

					idx = np.concatenate(([idx[0]-1],idx))

				else:
					# avoid going out of dataframe
					if idx[-1]+1 == len(Diff):
						break
					# avoid changing chromosome
					if Diff.loc[idx[-1]+1].chr != Diff.loc[idx[-1]].chr:
						break
					# only extend in neutral regions
					if Diff.loc[idx[-1]+1,'domain'] != 0:
						break

					idx = np.concatenate((idx,[idx[-1]+1]))

				# Compute updated q1
				if std_err:
					p = np.sum(Diff.loc[idx,'norm_cdf'])
					n = np.sum(1-Diff.loc[idx,'norm_cdf'])
				else:
					n_nan = np.sum(np.isnan(Diff.loc[idx,'value']))
					p = np.sum(Diff.loc[idx,'value'] > 0) + n_nan/2
					n = np.sum(Diff.loc[idx,'value'] < 0) + n_nan/2
				q1 = beta.ppf(percentile,p+1,n+1)

				if domain_sign == 1:
					cond = (q1>=q0)
				else:
					cond = (q1<=q0)

				if cond:
					n_bin += 1
					if direction == 'r':
						d_new.end += n_bin*bin_size
					else:
						d_new.start -= n_bin*bin_size
					n_bin = 0
					counter = 0
					q0=q1
				else:
					counter += 1
					n_bin += 1

				if counter >= n_bin_bridge:
					break
			# end while (extend domain)

			# update Diff_domain
			if direction == 'l':
				# updade start
				Diff_domain.loc[domain_idx,'start'] = d_new.start
				# update previous domain only if same chr
				if domain_idx > Diff_domain.index[0]:
					if Diff_domain.loc[domain_idx,'chr'] == Diff_domain.loc[domain_idx-1,'chr']:
						Diff_domain.loc[domain_idx-1,'end'] = d_new.start

			else:
				# update end
				Diff_domain.loc[domain_idx,'end'] = d_new.end
				# update next domain only if not end and same chr
				if domain_idx < Diff_domain.index[-1]:
					if Diff_domain.loc[domain_idx,'chr'] == Diff_domain.loc[domain_idx+1,'chr']:
						Diff_domain.loc[domain_idx+1,'start'] = d_new.end

			# update Diff
			idx_chr = np.array( Diff['chr'] == d_new.chr )
			idx_start = np.array( Diff['start'] >= d_new.start )
			idx_end = np.array( Diff['end'] <= d_new.end )
			Diff.loc[ idx_chr*idx_start*idx_end ,'domain'] = d_new.domain

		#end direction loop
	# end loop on domains

