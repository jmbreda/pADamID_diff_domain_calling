import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import beta, norm

x=np.arange(0,1,.001)
data = {}
data['neut'] = [20,30]
data['pos'] = [48,2]
data['neg'] = [8,42]
percentiles = [0.001,0.999]

fig = plt.figure()
rows = 1
cols = 3
f=0
for i in data.keys():
	f+=1
	ax = fig.add_subplot(rows,cols,f)
	plt.plot(x,beta.pdf(x,data[i][0]+1,data[i][1]+1))
	q = beta.ppf(percentiles,data[i][0]+1,data[i][1]+1)
	plt.plot(q,[0,0],'o-')
	plt.plot(0.5,0,'ro')
	ax.set_xlabel('$p$')
	ax.set_ylabel('$P(p|k,n)$')
	ax.set_title(f'{i}. region')
	ax.set_yticks([])
	ax.set_xticks(np.arange(0,1.01,.25))

fig.set_size_inches([3*cols,2.5*rows])
plt.tight_layout()
fig.savefig(f'Fig/beta_cartoon.pdf')
plt.close(fig)
