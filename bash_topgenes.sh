for i in cell_line indication
do
	# part 2: top 100 genes
	cd experiment_test_by_${i}

	# top genes
	#for j in 5 10 20 30 40 50 60 70 80 90 100 150 250 500 1000
	#do

	#cp -r  ../master/master_code/ features_monotherapy_moa_molecular_target_gene_network_geneset_chemical_structure_drug_name_topgene${j}
	#cd features_monotherapy_moa_molecular_target_gene_network_geneset_chemical_structure_drug_name_topgene${j}
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f monotherapy moa molecular target_gene network geneset chemical_structure drug_name --surrogate_gene ${j} --hold_out_test 
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f monotherapy moa molecular target_gene network geneset chemical_structure drug_name --surrogate_gene ${j} --hold_out_test
	#cd ..

	#done

	#cp -r  ../master/master_code/ features_monotherapy_moa_molecular_target_gene_network_geneset_chemical_structure_drug_name_topgene0
	#cd features_monotherapy_moa_molecular_target_gene_network_geneset_chemical_structure_drug_name_topgene0
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f monotherapy moa geneset chemical_structure drug_name --hold_out_test 
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f monotherapy moa geneset chemical_structure drug_name  --hold_out_test
	#cd ..

	#cp -r  ../master/master_code/ features_monotherapy_moa_molecular_target_gene_network_geneset_chemical_structure_drug_name_topgene2725
	#cd features_monotherapy_moa_molecular_target_gene_network_geneset_chemical_structure_drug_name_topgene2725
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f monotherapy moa molecular target_gene network geneset chemical_structure drug_name --hold_out_test 
	#screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f monotherapy moa molecular target_gene network geneset chemical_structure drug_name  --hold_out_test
	#cd ..

	# top synthetic lethality
	for j in 30 40 50 #5 10 15 20 25 30 40 50 80 100 200 500
	do
	cp -r  ../master/master_code/ features_monotherapy_moa_synthetic_lethality_target_gene_network_geneset_chemical_structure_drug_name_topsynleth${j}
	cd features_monotherapy_moa_synthetic_lethality_target_gene_network_geneset_chemical_structure_drug_name_topsynleth${j}
	screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f monotherapy moa geneset chemical_structure drug_name synthetic_lethality target_gene network --surrogate_synleth ${j} --hold_out_test 
	screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f monotherapy moa geneset chemical_structure drug_name synthetic_lethality target_gene network --surrogate_synleth ${j} --hold_out_test
	cd ..

	done

	# features_monotherapy_moa_geneset_molecular_synthetic_lethality_target_gene_network_chemical_structure_drug_name

	cd ..
done


#for i in cell_line indication
#do
	# part 2: top synthetic lethality
#	cd experiment_test_by_${i}
#
#	for j in 10 # 10 20 30 40 50 60 70 80 90 100 150 250 500 1000
#	do
#		for k in 10 #5 10 15 20 25 30 40 50 80 100 200
#		do
#			cp -r  ../master/master_code/ features_monotherapy_moa_molecular_target_gene_network_geneset_chemical_structure_drug_name_synthetic_lethality_topgene${j}_topsynleth${k}
#			cd features_monotherapy_moa_molecular_target_gene_network_geneset_chemical_structure_drug_name_synthetic_lethality_topgene${j}_topsynleth${k}
#			screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f monotherapy moa molecular target_gene network geneset chemical_structure drug_name synthetic_lethality --surrogate_gene ${j} --surrogate_synleth ${k} --hold_out_test 
#			screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f monotherapy moa molecular target_gene network geneset chemical_structure drug_name synthetic_lethality --surrogate_gene ${j} --surrogate_synleth ${k} --hold_out_test
#			cd ..
#		done
#	done

#	cd ..
#done

