import pandas as pd
import numpy as np


def get_pADamID_datasets():
	Datasets = {}

	Datasets['Top1'] = {}
	Datasets['Top1']['datafolder'] = '/DATA/scratch/usr/t.v.schaik/proj/sManzo_pADamID/ts190515_pADamID_RPE_Top1_DRB/results/normalized/bin-20kb'
	Datasets['Top1']['treatment'] = [f'pADamID-RPE_Top1_r{i}_Lmnb2-20kb.norm.txt.gz' for i in [2,3,4,6]]
	Datasets['Top1']['control'] = [f'pADamID-RPE_SC_r{i}_Lmnb2-20kb.norm.txt.gz' for i in [2,3,4,6]]

	Datasets['Top2B'] = {}
	Datasets['Top2B']['datafolder'] = '/DATA/scratch/usr/t.v.schaik/proj/sManzo_pADamID/ts190515_pADamID_RPE_Top1_DRB/results/normalized/bin-20kb'
	Datasets['Top2B']['treatment'] = [f'pADamID-RPE_Top2B_r{i}_Lmnb2-20kb.norm.txt.gz' for i in range(8,13)]
	Datasets['Top2B']['control'] = [f'pADamID-RPE_SCII_r{i}_Lmnb2-20kb.norm.txt.gz' for i in range(8,13)]

	Datasets['TPL'] = {}
	Datasets['TPL']['datafolder'] = '/DATA/scratch/usr/t.v.schaik/proj/sManzo_pADamID/ts190515_pADamID_RPE_Top1_DRB/results/normalized/bin-20kb'
	Datasets['TPL']['treatment'] = [f'pADamID-RPE_TPL2_r{i}_Lmnb2-20kb.norm.txt.gz' for i in [6,8,9]]+[f'pADamID-RPE_TPL2_CAL_r{i}_Lmnb2-20kb.norm.txt.gz' for i in [8,9]]
	Datasets['TPL']['control'] = [f'pADamID-RPE_CT_r{i}_Lmnb2-20kb.norm.txt.gz' for i in [6,8,9]]+[f'pADamID-RPE_CT_CAL_r{i}_Lmnb2-20kb.norm.txt.gz' for i in [8,9]]

	Datasets['FLV'] = {}
	Datasets['FLV']['datafolder'] = '/DATA/scratch/usr/t.v.schaik/proj/sManzo_pADamID/ts190515_pADamID_RPE_Top1_DRB/results/normalized/bin-20kb'
	Datasets['FLV']['treatment'] = ['pADamID-RPE_FLV2-3h_r6_Lmnb2-20kb.norm.txt.gz','pADamID-RPE_FLV_r4_Lmnb2-20kb.norm.txt.gz']
	Datasets['FLV']['control'] = ['pADamID-RPE_CT_r6_Lmnb2-20kb.norm.txt.gz','pADamID-RPE_CT_r4_Lmnb2-20kb.norm.txt.gz']

	Datasets['NegCtrl'] = {}
	Datasets['NegCtrl']['datafolder'] = '/DATA/scratch/usr/t.v.schaik/proj/sManzo_pADamID/ts190515_pADamID_RPE_Top1_DRB/results/normalized/bin-20kb'
	Datasets['NegCtrl']['treatment'] = [f'pADamID-RPE_SCII_r{i}_Lmnb2-20kb.norm.txt.gz' for i in range(8,11)]
	Datasets['NegCtrl']['control'] = [f'pADamID-RPE_SCII_r{i}_Lmnb2-20kb.norm.txt.gz' for i in range(11,13)]

	Datasets['Top2B_2'] = {}
	Datasets['Top2B_2']['datafolder'] = '/DATA/usr/j.breda/Workspace/jb20220118_Stefano_pADamID/data/input_pADamID'
	Datasets['Top2B_2']['treatment'] = [f'pADamID-RPE_WT_Top2B_r{i}_Lmnb2-20kb.norm.txt.gz' for i in range(13,17)]
	Datasets['Top2B_2']['control'] = [f'pADamID-RPE_WT_SCII_r{i}_Lmnb2-20kb.norm.txt.gz' for i in range(13,17)]


	return Datasets


def get_pADamID_tables(Datasets):

	for i,f in enumerate(Datasets['treatment']):
		if i==0:
			Treat = pd.read_csv(Datasets['datafolder']+'/'+f,sep='\t',header=None)
			X = Treat.iloc[:,3].values[:,None]
			Treat.drop(columns=3,inplace=True)
		else:
			tmp = pd.read_csv(Datasets['datafolder']+'/'+f,sep='\t',header=None)
			if not np.all(Treat.values == tmp.iloc[:,:3].values):
				print('error: chr-start-end unequal in replicates')
				break
			X = np.concatenate((X,tmp.iloc[:,3].values[:,None]),axis=1)

	n = len(Datasets['treatment'])
	if n>1:
		Treat['value'] = np.zeros(X.shape[0])
		Treat['value_std_err'] = np.zeros(X.shape[0])
		idx_nan = np.sum(np.isnan(X),axis=1)>=(n-1)
		Treat.iloc[idx_nan,3:] = np.nan
		Treat.iloc[~idx_nan,3] = np.nanmean(X[~idx_nan,:],axis=1)
		Treat.iloc[~idx_nan,4] = np.nanstd(X[~idx_nan,:],axis=1)/np.sqrt(n)
		Treat.columns = ['chr','start','end','value','value_std_err']
	else:
		Treat['value'] = X
		Treat.columns = ['chr','start','end','value']

	for i,f in enumerate(Datasets['control']):
		if i==0:
			Ctrl = pd.read_csv(Datasets['datafolder']+'/'+f,sep='\t',header=None)
			X = Ctrl.iloc[:,3].values[:,None]
			Ctrl.drop(columns=3,inplace=True)
		else:
			tmp = pd.read_csv(Datasets['datafolder']+'/'+f,sep='\t',header=None)
			if not np.all(Ctrl.values == tmp.iloc[:,:3].values):
				print('error: chr-start-end unequal in replicates')
				break
			X = np.concatenate((X,tmp.iloc[:,3].values[:,None]),axis=1)

	n = len(Datasets['control'])
	if n>1:
		Ctrl['value'] = np.zeros(X.shape[0])
		Ctrl['value_std_err'] = np.zeros(X.shape[0])
		idx_nan = np.sum(np.isnan(X),axis=1)>=(n-1)
		Ctrl.iloc[idx_nan,3:] = np.nan
		Ctrl.iloc[~idx_nan,3] = np.nanmean(X[~idx_nan,:],axis=1)
		Ctrl.iloc[~idx_nan,4] = np.nanstd(X[~idx_nan,:],axis=1)/np.sqrt(n)
		Ctrl.columns = ['chr','start','end','value','value_std_err']
	else:
		Ctrl['value'] = X
		Ctrl.columns = ['chr','start','end','value']


	return Treat, Ctrl



