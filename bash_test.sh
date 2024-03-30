# shap analysis on the final model
for i in cell_line indication
do
    cd experiment_test_by_${i}
    #molecular+target_gene+network done
    for j in exp cnv snv lof coh_pat lof_pat ddr
    do
	    cp -r  ../master/master_code/ features_molecular_target_gene_network_${j}
	    cd features_molecular_target_gene_network_${j}
	    screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s aoc -f molecular target_gene network --mol_type ${j} --hold_out_test
	    screen -d -m python main.py -p ../../test_by_${i} --exclude_synergy_batch --exclude_cell_line --exclude_cancer_type -s bliss -f molecular target_gene network --mol_type ${j} --hold_out_test
	    cd ..
    done

    
    cd ..
done
