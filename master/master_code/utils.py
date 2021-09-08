#!usr/bin/env/python3
#author: @rayezh
import pandas as pd
import numpy as np
import numpy.ma as ma
from collections import defaultdict

class encode_categorical_feature():
	""" Categorical feature encoding object
	
	Params
	------
	data: '../../feature/data/synergy_responses_with_monotherapy.tsv' pandas frame
	
	Yields
	------
	batch: dict
	cell_line: dict
	moa: dict
	remark: dict
	cancer_type: dict
	cancer_subtype: dict
	"""
	data = pd.read_csv('../../feature/data/synergy_responses_with_monotherapy.tsv', sep = '\t')
	
	def __init__(self):
		self.batch = defaultdict(lambda:np.nan, {j:i for i, j in enumerate(sorted(set(self.data['.identifier_batch'])))})
		self.cell_line = defaultdict(lambda:np.nan, {j:i for i, j in enumerate(sorted(set(self.data['.identifier_sample_name'])))})
		self.moa = defaultdict(lambda:np.nan, {j:i for i, j in enumerate(sorted(set([str(x) for x in self.data['.metadata_moa_1'].to_list()+self.data['.metadata_moa_2'].to_list()])))})
		self.treatment = defaultdict(lambda:np.nan, {j:i for i, j in enumerate(sorted(set([str(x) for x in self.data['.metadata_treatment_1'].to_list()+self.data['.metadata_treatment_2'].to_list()])))})
		self.remark = defaultdict(lambda:np.nan, {j:i for i, j in enumerate(sorted(set([str(x) for x in self.data['.metadata_treatment_remarks']])))})
		self.cancer_type = defaultdict(lambda:np.nan, {j:i for i, j in enumerate(sorted(set(self.data['.metadata_cancer_type'])))})
		self.cancer_subtype = defaultdict(lambda:np.nan, {j:i for i, j in enumerate(sorted(set(self.data['.metadata_cancer_subtype'])))})

def boostrapping_confidence_interval(pred_all, ci):
	""" Boostrapping to get a 95 confidence interval for prediction performance
	
	Params
	------
	pred_all: Numpy array
	ci: confidence interval
	
	Yields
	------
	cor_mean: float 
		middle bound
	lb: float
		lower bound
	ub: float
	upper bound
	
	"""
	import random
	# set random seed
	random.seed(0)
	
	# calculate overall correlation
	cor_mean = ma.corrcoef(ma.masked_invalid(pred_all[:,0]), ma.masked_invalid(pred_all[:,1]))[0,1]
	print("Overall prediction/gold standard correlation is %.4f" % cor_mean)
	# start boostrapping ...
	cor_all = [] 
	for i in range(100):
		pred_new = random.choices(pred_all, k = len(pred_all))
		pred_new = np.array(pred_new)
		cor = ma.corrcoef(ma.masked_invalid(pred_new[:,0]), ma.masked_invalid(pred_new[:,1]))[0,1]
		cor_all.append(cor)
	cor_all = sorted(cor_all)
	
	lb = cor_all[round(100*(0.5-ci*0.5))]
	ub = cor_all[round(100*(0.5+ci*0.5))]
	print("%d%s Confidence interval is: (%.4f, %.4f)" % (int(ci*100), '%', lb, ub))
	
	return cor_mean, lb, ub





