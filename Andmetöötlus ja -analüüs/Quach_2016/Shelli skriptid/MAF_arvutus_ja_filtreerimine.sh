#!/bin/bash

bcftools view ${inputname} -S population_1_names.txt > ${tempName}
bcftools +fill-tags ${tempname} -Ov -o ${mafname} -- -t MAF
bcftools filter ${mafname} -i 'MAF[0] > 0.01' > ${outputname}
bgzip ${outputname}
bcftools index ${outputname}.gz
bcftools view ${outputname}.gz -r 9:20318520-21318520 > ${outputname}_IFNB1.vcf
