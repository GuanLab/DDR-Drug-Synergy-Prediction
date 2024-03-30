for i in cell_line indication
do
	cd experiment_test_by_${i}

	#monotherapy+mode-of-action+chemical_structure+drug_name_geneset+molecular
	#cp -r  ../master/master_code/ features_monotherapy_moa_chemical_structure_drug_name_geneset_molecular
	#cd features_monotherapy_moa_chemical_structure_drug_name_geneset_molecular
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f monotherapy moa chemical_structure drug_name geneset molecular --hold_out_test 
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f monotherapy moa chemical_structure dru_name geneset molecular --hold_out_test 
	#cd ..

	#monotherapy+mode-of-action+chemical_structure+drug_name_geneset+molecular_target_gene
	#cp -r  ../master/master_code/ features_monotherapy_moa_chemical_structure_drug_name_geneset_molecular_target_gene
	cd features_monotherapy_moa_chemical_structure_drug_name_geneset_molecular_target_gene
	screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f monotherapy moa chemical_structure drug_name geneset molecular target_gene --hold_out_test 
	screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f monotherapy moa chemical_structure dru_name geneset molecular target_gene --hold_out_test 
	cd ..
	
	#monotherapy+mode-of-action+chemical_structure+drug_name_geneset+molecular_target_gene_network
	#cp -r  ../master/master_code/ features_monotherapy_moa_chemical_structure_drug_name_geneset_molecular_target_gene_network
	#cd features_monotherapy_moa_chemical_structure_drug_name_geneset_molecular_target_gene_network
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f monotherapy moa chemical_structure drug_name geneset molecular target_gene network --hold_out_test 
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f monotherapy moa chemical_structure dru_name geneset molecular target_gene network --hold_out_test 
	#cd ..

	# TODO: synthetic lethality
#	cp -r  ../master/master_code/ features_synthetic_lethality
#	cd features_synthetic_lethality
#	screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f synthetic_lethality --hold_out_test
#	screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f synthetic_lethality --hold_out_test
#	cd ..
	
	cd ..
done