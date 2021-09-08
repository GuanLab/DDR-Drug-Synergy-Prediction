for i in cell_line indication
do
	mkdir experiment_test_by_${i}
	cd experiment_test_by_${i}

	############################ total set ot features ############################
	cp -r  ../master/master_code features_monotherapy_moa_molecular_target_gene_geneset_chemical_structure
	cd features_monotherapy_moa_molecular_target_gene_geneset_chemical_structure
	python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f monotherapy moa molecular target_gene geneset chemical_structure --evaluate_shap Train Test
	python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f monotherapy moa molecular target_gene geneset chemical_structure --evaluate_shap Train Test
	cd downstream_analysis
	python statistical_analysis.py
	python select_surrogate_features.py # get top genes/geneset for surrogate models
	cd ../../
	
	cp -r  ../master/master_code features_monotherapy_moa_molecular_target_gene_geneset_chemical_structure_drug_name
	cd features_monotherapy_moa_molecular_target_gene_geneset_chemical_structure_drug_name
	python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f monotherapy moa molecular target_gene geneset chemical_structure drug_name --evaluate_shap Train Test
	python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f monotherapy moa molecular target_gene geneset chemical_structure drug_name --evaluate_shap Train Test
	cd ..
	
	############################ total set ot features (excluding monotherapy) ############################
	cp -r  ../master/master_code features_moa_molecular_target_gene_geneset_chemical_structure
	cd features_moa_molecular_target_gene_geneset_chemical_structure
	python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f moa molecular target_gene geneset chemical_structure --evaluate_shap Train Test
	python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f moa molecular target_gene geneset chemical_structure --evaluate_shap Train Test
	cd ..

    cp -r  ../master/master_code features_moa_molecular_target_gene_geneset_chemical_structure_drug_name
    cd features_moa_molecular_target_gene_geneset_chemical_structure_drug_name
    python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f moa molecular target_gene geneset chemical_structure drug_name --evaluate_shap Train Test
    python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f moa molecular target_gene geneset chemical_structure drug_name --evaluate_shap Train Test
    cd ..
	cd ..
done

