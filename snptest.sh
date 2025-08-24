PATH=$PATH://plink-1.07-mac-intel/

plink --file sim --make-bed --out sim

./snptest_v2.5.6 \
-data sim.bed sim.sample \
-frequentist 1 \
-pheno phenotype \
-method newml \
-o snptest.out
