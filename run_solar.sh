#!/bin/bash

# target=intensity_asy_total_ca
target=${1} 
wd=../hcp/solar/${target}

rm -rf ${wd}/*INORM
rm -rf ${wd}/header
rm -rf ${wd}/*inorm
rm -rf ${wd}/*result*
rm -rf ${wd}/*.mod
rm -rf ${wd}/phenotypes.info

cp ../hcp/solar/pedi_shell/* ${wd}/

if [[ ${target} == *"ca"* ]]; then
  List=$(seq 1 12)
else
  List=$(seq 1 180)
fi

for i in ${List}
do
  echo node_${i} >> ${wd}/header
done

pheno_list=`more ${wd}/header`
phenotype="1_phen.csv"

# Lets make solar_run file
out_file=${wd}/solar.run.inorm
echo pheno load ${phenotype} > ${out_file}
echo model new >> ${out_file}
echo covar age^1,2#sex >> ${out_file}

for cur_dir in $pheno_list; do
    echo define ${cur_dir}'_'INORM = inorm_${cur_dir} >> ${out_file}
    echo trait ${cur_dir}'_'INORM >> ${out_file}
    echo polyg -all -s >> ${out_file}
done

# activate solar SOLARECL

cd $wd
solar < solar.run.inorm

pheno_list=`find . -type d -name "node*"`

# Lets make solar_run file

echo node","name","H2r","SE","rp >> "1_heri_result.csv"

for cur_dir  in  $pheno_list
do
echo $cur_dir

XX=`echo $cur_dir | awk -F "_" '{print $2}'`
echo $XX

X=`echo $cur_dir`
echo $X

line=`more $cur_dir/polygenic.out | grep "H2r is"`
H2r=`echo $line |  awk '{print $3}'`
echo $H2r

rp=`echo $line |  awk '{print $6}'`
echo $rp

line=`more $cur_dir/polygenic.out | grep "H2r Std. Error: "`
SE=`echo $line |  awk '{print $4}'`
echo $SE

echo $XX","$X","$H2r","$SE","$rp >> "1_heri_result.csv"

done
