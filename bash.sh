for i in cell_line indication
do
	mkdir experiment_test_by_${i}
	cd experiment_test_by_${i}

	# part 1: each feature

	#monotherapy done
	#cp -r  ../master/master_code/ features_monotherapy
	#cd features_monotherapy
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f monotherapy --hold_out_test 
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f monotherapy --hold_out_test 
	#cd ..

	#monotherapy+mode-of-action done
	#cp -r  ../master/master_code/ features_monotherapy_moa
	#cd features_monotherapy_moa
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f monotherapy moa --hold_out_test 
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f monotherapy moa --hold_out_test 
	#cd ..

	#monotherapy+mode-of-action+geneset done
	#cp -r  ../master/master_code/ features_monotherapy_moa_geneset
	#cd features_monotherapy_moa_geneset
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f monotherapy moa geneset --hold_out_test 
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f monotherapy moa geneset --hold_out_test 
	#cd ..

	#monotherapy+mode-of-action+geneset+chemical_structure
	#cp -r  ../master/master_code/ features_monotherapy_moa_geneset_chemical_structure
	#cd features_monotherapy_moa_geneset_chemical_structure
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f monotherapy moa geneset chemical_structure --hold_out_test 
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f monotherapy moa geneset chemical_structure --hold_out_test 
	#cd ..

	#monotherapy+mode-of-action+geneset+chemical_structure+drug_name
	#cp -r  ../master/master_code/ features_monotherapy_moa_geneset_chemical_structure_drug_name
	#cd features_monotherapy_moa_geneset_chemical_structure_drug_name
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f monotherapy moa geneset chemical_structure drug_name --hold_out_test 
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f monotherapy moa geneset chemical_structure drug_name --hold_out_test 
	#cd ..

	#drug_name done
	#cp -r  ../master/master_code/ features_drug_name
	#cd features_drug_name
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f drug_name --hold_out_test --evaluate_shap Train Test
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f drug_name --hold_out_test --evaluate_shap Train Test
	#cd ..

	#molecular done
	#cp -r  ../master/master_code/ features_molecular
	#cd features_molecular
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f molecular --hold_out_test 
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f molecular --hold_out_test 
	#cd ..
	
	#molecular+target_gene done
	#cp -r  ../master/master_code/ features_molecular_target_gene
	#cd features_molecular_target_gene
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f molecular target_gene --hold_out_test 
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f molecular target_gene --hold_out_test 
	#cd ..

	#molecular+target_gene+network done
	#cp -r  ../master/master_code/ features_molecular_target_gene_network
	#cd features_molecular_target_gene_network
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f molecular target_gene network --hold_out_test
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f molecular target_gene network --hold_out_test
	#cd ..

	#synthetic_lethality done
	#cp -r  ../master/master_code/ features_synthetic_lethality
	#cd features_synthetic_lethality
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f synthetic_lethality --hold_out_test 
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f synthetic_lethality --hold_out_test 
	#cd ..
	
	#synthetic_lethality+target_gene done
	#cp -r  ../master/master_code/ features_synthetic_lethality_target_gene
	#cd features_synthetic_lethality_target_gene
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f synthetic_lethality target_gene --hold_out_test 
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f synthetic_lethality target_gene --hold_out_test 
	#cd ..

	#synthetic_lethality+target_gene+network done
	#cp -r  ../master/master_code/ features_synthetic_lethality_target_gene_network
	#cd features_synthetic_lethality_target_gene_network
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f synthetic_lethality target_gene network --hold_out_test
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f synthetic_lethality target_gene network --hold_out_test
	#cd ..

	#molecular+synthetic_lethality+target_gene+network done
	#cp -r  ../master/master_code/ features_molecular_synthetic_lethality_target_gene_network
	#cd features_molecular_synthetic_lethality_target_gene_network
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f molecular synthetic_lethality target_gene network --hold_out_test
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f molecular synthetic_lethality target_gene network --hold_out_test
	#cd ..

	#molecular+synthetic_lethality+target_gene+network+drug_name done
	#cp -r  ../master/master_code/ features_molecular_synthetic_lethality_target_gene_network_drug_name
	#cd features_molecular_synthetic_lethality_target_gene_network_drug_name
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f molecular synthetic_lethality target_gene network drug_name --hold_out_test
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f molecular synthetic_lethality target_gene network drug_name --hold_out_test
	#cd ..

	#molecular+synthetic_lethality+target_gene+network+drug_name+chemical_structure done
	#cp -r  ../master/master_code/ features_molecular_synthetic_lethality_target_gene_network_drug_name_chemical_structure
	#cd features_molecular_synthetic_lethality_target_gene_network_drug_name_chemical_structure
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f molecular synthetic_lethality target_gene network drug_name chemical_structure --hold_out_test
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f molecular synthetic_lethality target_gene network drug_name chemical_structure --hold_out_test
	#cd ..

	#molecular+synthetic_lethality+target_gene+network+drug_name+chemical_structure+monotherapy done
	#cp -r  ../master/master_code/ features_molecular_synthetic_lethality_target_gene_network_drug_name_chemical_structure_monotherapy
	#cd features_molecular_synthetic_lethality_target_gene_network_drug_name_chemical_structure_monotherapy
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f molecular synthetic_lethality target_gene network drug_name chemical_structure monotherapy --hold_out_test
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f molecular synthetic_lethality target_gene network drug_name chemical_structure monotherapy --hold_out_test
	#cd ..

	#molecular+synthetic_lethality+target_gene+network+drug_name+chemical_structure+monotherapy+mode-of-action done
	#cp -r  ../master/master_code/ features_molecular_synthetic_lethality_target_gene_network_drug_name_chemical_structure_monotherapy_moa
	#cd features_molecular_synthetic_lethality_target_gene_network_drug_name_chemical_structure_monotherapy_moa
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f molecular synthetic_lethality target_gene network drug_name chemical_structure monotherapy moa --hold_out_test
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f molecular synthetic_lethality target_gene network drug_name chemical_structure monotherapy moa --hold_out_test
	#cd ..


	#monotherapy+mode-of-action+geneset+molecular done
	#cp -r  ../master/master_code/ features_monotherapy_moa_geneset_molecular
	#cd features_monotherapy_moa_geneset_molecular
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f monotherapy moa geneset molecular --hold_out_test
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f monotherapy moa geneset molecular --hold_out_test
	#cd ..

	#monotherapy+mode-of-action+geneset+molecular+target_gene 
	#cp -r  ../master/master_code/ features_monotherapy_moa_geneset_molecular_target_gene
	#cd features_monotherapy_moa_geneset_molecular_target_gene
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f monotherapy moa geneset molecular target_gene --hold_out_test
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f monotherapy moa geneset molecular target_gene --hold_out_test
	#cd ..

	#monotherapy+mode-of-action+geneset+molecular+target_gene+network done
	#cp -r  ../master/master_code/ features_monotherapy_moa_geneset_molecular_target_gene_network
	#cd features_monotherapy_moa_geneset_molecular_target_gene_network
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f monotherapy moa geneset molecular target_gene network --hold_out_test
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f monotherapy moa geneset molecular target_gene network --hold_out_test
	#cd ..

	#monotherapy+mode-of-action+geneset+molecular+synthetic_lethality+target_gene+network done
	#cp -r  ../master/master_code/ features_monotherapy_moa_geneset_molecular_synthetic_lethality_target_gene_network
	#cd features_monotherapy_moa_geneset_molecular_synthetic_lethality_target_gene_network
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f monotherapy moa geneset molecular synthetic_lethality target_gene network --hold_out_test
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f monotherapy moa geneset molecular synthetic_lethality target_gene network --hold_out_test
	#cd ..

	#monotherapy+mode-of-action+geneset+molecular+synthetic_lethality+target_gene+network+chemical_structure done
	#cp -r  ../master/master_code/ features_monotherapy_moa_geneset_molecular_synthetic_lethality_target_gene_network_chemical_structure
	#cd features_monotherapy_moa_geneset_molecular_synthetic_lethality_target_gene_network_chemical_structure
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type --exclude_cancer_type -s aoc -f monotherapy moa geneset molecular synthetic_lethality target_gene network chemical_structure --hold_out_test
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type --exclude_cancer_type -s bliss -f monotherapy moa geneset molecular synthetic_lethality target_gene network chemical_structure --hold_out_test
	#cd ..

	#monotherapy+mode-of-action+geneset+molecular+synthetic_lethality+target_gene+network+chemical_structure+drug_name done
	#cp -r  ../master/master_code/ features_monotherapy_moa_geneset_molecular_synthetic_lethality_target_gene_network_chemical_structure_drug_name
	#cd features_monotherapy_moa_geneset_molecular_synthetic_lethality_target_gene_network_chemical_structure_drug_name
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type --exclude_cancer_type -s aoc -f monotherapy moa geneset molecular synthetic_lethality target_gene network chemical_structure drug_name --hold_out_test --evaluate_shap Train Test
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type --exclude_cancer_type -s bliss -f monotherapy moa geneset molecular synthetic_lethality target_gene network chemical_structure drug_name --hold_out_test --evaluate_shap Train Test
	#cd ..

	
	cd ..
done

