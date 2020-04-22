#!/bin/bash

args=("$@")
ALIGN_DIR=${args[0]}
PEAK_DIR=${args[1]}
output=${args[2]}

sample_number=$(($#-1))
name_array=(" ")
score_array=("frip")

for i in $(seq 3 $sample_number)
do
   name_array+=(${args[i]})
   sample_peak=$PEAK_DIR"/"${args[i]}"_peaks.filtered.narrowPeak"
   sample_bam=$ALIGN_DIR"/"${args[i]}".filtered.dedup.sorted.bam"
   NR1=$(sambamba view $sample_bam | wc -l)
   RIP1_SELF=$(sambamba view -L $sample_peak $sample_bam | wc -l)
   score=$(awk -v var1=$RIP1_SELF -v var2=$NR1 'BEGIN { print  ( var1 / var2 ) }')
   score_array+=($score)
done

printf "%s\t" "${name_array[@]}" > $output
printf "\n" >> $output
printf "%s\t" "${score_array[@]}" >> $output
