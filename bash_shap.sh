# shap analysis on the final model
for i in cell_line indication
do
    cd experiment_test_by_${i}
	#monotherapy+mode-of-action+geneset+molecular+synthetic_lethality+target_gene+network+chemical_structure+drug_name done
	cp -r  ../master/master_code/*py features_monotherapy_moa_geneset_molecular_synthetic_lethality_target_gene_network_chemical_structure_drug_name
	cd features_monotherapy_moa_geneset_molecular_synthetic_lethality_target_gene_network_chemical_structure_drug_name
	screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type --exclude_cancer_type -s aoc -f monotherapy moa geneset molecular synthetic_lethality target_gene network chemical_structure drug_name --evaluate_shap Train Test --skip_training
	screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type --exclude_cancer_type -s bliss -f monotherapy moa geneset molecular synthetic_lethality target_gene network chemical_structure drug_name --evaluate_shap Train Test --skip_training
	cd ..
    
    cd ..
done
