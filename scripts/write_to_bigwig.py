import numpy as np
import pyBigWig as bw
import pandas as pd


def write_bw(df,outfile):
	CHR = np.unique(df.chr)

	# Create diff bw file
	bw_fid = bw.open(outfile,'w')

	my_head = [(c, df[df['chr']==c].iloc[-1,2]) for c in CHR]
	bw_fid.addHeader(my_head)
	for c in CHR:
		my_chr = df[df.chr==c].chr.values.astype(str)
		start = df[df.chr==c].start.values
		end = df[df.chr==c].end.values-1
		vals = df[df.chr==c].iloc[:,-1].values
		bw_fid.addEntries(my_chr, start, ends=end, values=vals)
	bw_fid.close()

def write_bed(df,outfile):

	fid_neg = open(f'{outfile}_negative_regions.bed','w')
	fid_0 = open(f'{outfile}_neutral_regions.bed','w')
	fid_pos = open(f'{outfile}_positive_regions.bed','w')
	for l in df.values:
		if l[3] == -1:
			fid_neg.write('%s\t%u\t%u\n'%(l[0],l[1],l[2]))
		elif l[3] == 0:
			fid_0.write('%s\t%u\t%u\n'%(l[0],l[1],l[2]))
		elif l[3] == 1:
			fid_pos.write('%s\t%u\t%u\n'%(l[0],l[1],l[2]))
	fid_neg.close()
	fid_0.close()
	fid_pos.close()






# plot diff
if False:
	fig = plt.figure()
	rows = len(CHR)
	cols = 1
	f=0

	for i,c in enumerate(CHR[:-1]):
		print(c)
		f+=1
		ax = fig.add_subplot(rows,cols,f)
		ax.bar(np.arange(len(E['diff'][c])),E['diff'][c],1)
		ax.set_title(c)

	fig.set_size_inches([12*cols,1*rows])
	plt.tight_layout()


def running_mean(x,f,win):
  _x_ = np.array([ np.mean(x[i:i+win]) for i in range(int(len(x)-win+1)) ])
  _f_ = np.array([ np.nanmean(f[i:i+win]) for i in range(int(len(f)-win+1)) ])
  _df_ = np.array([ np.nanstd(f[i:i+win]) for i in range(int(len(f)-win+1)) ])

  return _x_, _f_, _df_

if False:
	fig = plt.figure()
	rows = len(CHR)
	cols = 1
	f=0
	# get running mean
	for c in CHR:
		print(c)
		bins = np.arange(len(x))*bin_size

		[x,Delta,dDelta] = running_mean(bins,E['diff'][c],100)

		f+=1
		ax = fig.add_subplot(rows,cols,f)

		ax.plot(bins,Delta,color='r')
		ax.fill_between(x, Delta - dDelta, Delta + dDelta, color='r',alpha=0.2)
		ax.grid('on')
		ax.set_yticks([0])

	fig.set_size_inches([20*cols,3*rows])
	plt.tight_layout()
	fig.savefig('Fig/moving_average.pdf')



# plot histogramms
if False:
	fig = plt.figure()
	rows = 3
	cols = 1
	f=0
	cm = plt.get_cmap('gist_ncar',len(CHR)+1)

	for cond in COND:
		f+=1
		ax = fig.add_subplot(rows,cols,f)
		for i,c in enumerate(CHR[:-1]):
			[h,x] = np.histogram(E[cond][c][~np.isnan(E[cond][c])],bins=30,density=1)
			ax.plot((x[:-1] + x[1:])/2,h,color=cm(i))

		if cond==COND[0]:
			ax.legend(CHR[:-1])
		ax.set_title(cond)
		ax.grid('on')

	fig.set_size_inches([12*cols,8*rows])
	plt.tight_layout()
	fig.savefig('Fig/hist_expr.pdf')





