PATH=$PATH://plink-1.07-mac-intel/

plink --file sim --make-bed --out sim

plink --bfile sim --assoc --adjust  --out simout

