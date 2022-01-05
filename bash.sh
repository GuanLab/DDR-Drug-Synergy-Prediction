cd feature/QC
python QC.py
python pull_drug_target.py
python split_train_test.py
cd ../../

for i in cell_line indication
do
	mkdir experiment_test_by_${i}
	cd experiment_test_by_${i}

	cp -r  ../master/master_code features_monotherapy_moa_molecular_target_gene_geneset_chemical_structure_drug_name
	cd features_monotherapy_moa_molecular_target_gene_geneset_chemical_structure_drug_name
	python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f monotherapy moa molecular target_gene geneset chemical_structure drug_name --evaluate_shap Train Test
	python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f monotherapy moa molecular target_gene geneset chemical_structure drug_name --evaluate_shap Train Test
	cd ../../
done

