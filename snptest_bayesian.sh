PATH=$PATH://plink-1.07-mac-intel/

plink --file sim --make-bed --out sim

./snptest_v2.5.6 \
-data sim.bed sim.sample \
-bayesian 1 \
-method score \
-pheno phenotype \
-o snptest_out
