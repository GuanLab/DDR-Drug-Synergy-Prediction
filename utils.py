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
	data: '../../feature/data/synergy_responses.tsv' pandas frame
	
	Yields
	------
	batch: dict
	cell_line: dict
	moa: dict
	remark: dict
	cancer_type: dict
	cancer_subtype: dict
	"""
	data = pd.read_csv('../../feature/data/synergy_responses.tsv', sep = '\t')

	
	def __init__(self):
		def one_hot_encode(col):
			""" One hot encode categorical feature
		
			Params
			------
			col: pandas series
		
			Yields
			------
			one_hot: list of encoded categorical feature
			fnames: list
			"""
			tmp = pd.get_dummies(col.unique(), dtype =int)
			fnames = tmp.columns
			one_hot = {col.unique()[i]:r.tolist() for i, r in enumerate(tmp.values)}

			return one_hot, fnames

		self.batch, self.batch_name = one_hot_encode(self.data['.identifier_batch'])
		self.cell_line, self.cell_line_name = one_hot_encode(self.data['.identifier_sample_name'])
		# combine moa_1 and moa_2
		tmp = pd.concat([self.data[['.metadata_moa_1']].rename(columns = {'.metadata_moa_1':'.metadata_moa'}), self.data[['.metadata_moa_2']].rename(columns = {'.metadata_moa_2':'.metadata_moa'})])
		self.moa, self.moa_name = one_hot_encode(tmp['.metadata_moa'])
		# combine treatment_1 and treatment_2
		tmp = pd.concat([self.data[['.metadata_treatment_1']].rename(columns = {'.metadata_treatment_1':'.metadata_treatment'}), self.data[['.metadata_treatment_2']].rename(columns = {'.metadata_treatment_2':'.metadata_treatment'})])
		self.treatment, self.treatment_name = one_hot_encode(tmp['.metadata_treatment']) 
		self.remark, self.remark_name = one_hot_encode(self.data['.metadata_treatment_remarks'])
		self.cancer_type, self.cancer_type_name = one_hot_encode(self.data['.metadata_cancer_type'])
		self.cancer_subtype, self.cancer_subtype_name = one_hot_encode(self.data['.metadata_cancer_subtype'])

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





