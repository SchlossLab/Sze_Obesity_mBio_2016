## master script to submit all prephasing jobs in parallel (per chromosome)
## will vary walltime request depending on scaling from chromosome length
## GWAS model.. to be saved in directory /g01/home/nfishb/LHS/GWAS_9
## slope ~ geno + PC1-5 + sex + age + fev0 + o'connor slope + BMI
## in continuous smokers only
## FEV0 is first non-missing observation 

cd ~
hours=($(<Jobs/SNPTEST_times.txt))
for i in `seq 1 22`; do
	let "ix=$i-1"
	qsub -l walltime=${hours[$ix]},nodes=1:ppn=6,pmem=2gb -q nestor -N gwas9_chr${i} -v chr=${i},dir=/g01/home/nfishb/LHS/GWAS_9,samp=phenotypes_GWAS9.sample,pheno=slope Jobs/GWAS_snptest.sh
done
