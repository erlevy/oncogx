#!/bin/bash

# Arg 1: bam file
# Arg 2: bed file
# Arg3: output directory
# Arg4: output directory intermedia files
# sample usage:
#   ./Bam_to_Fastq-fromTCGA-Pairend.sh ../NomenclatureTRBgenes/BAM_files/Bam_Files_TCGA/TCGA-B6-A0IJ-01A-11W-A050-09_IlluminaGA-DNASeq_exome_1.sorted.bam  /mnt/oncogxA/Valentina/NomenclatureTRBgenes/SureSelect/Bait-exome_SureSelect/TRA-B-D-G/Target-genesTR.ABDG-GenomicReg.bed outdirectory1 outdirectory2
region_file="$2" #"/mnt/oncogxA/Valentina/NomenclatureTRBgenes/SureSelect/Bait-exome_SureSelect/TRA-B-D-G/Target-genesTR.ABDG-GenomicReg.bed"
name=`basename $1 ".bam"`
bam_inpath=`dirname $1`
bam_outpath=/scratch
outdir=$3
outdir2=$4

#cp $1 /scratch
#cp $1.bai /scratch
#samtools index $bam_outpath/$name.bam

touch $bam_outpath/$name.summary.txt

# getting the flagstat from the raw bam file 
echo 'FLAGSTAT:' >> $bam_outpath/$name.summary.txt
samtools flagstat $bam_inpath/$name.bam >> $bam_outpath/$name.summary.txt 

# extracting overlapping+unmapped TCR reads
#mv $bam_outpath/$name.bam $bam_outpath/$name.TCRreg.bam
#samtools view -b -L $region_file $bam_outpath/$name.bam > $bam_outpath/$name.TCRreg.bam

 # double check proper assembly file was chosen
 # samtools view header | awk print AS:<ref assem> from second line |
 #   sed remove "AS:"
 #echo "Double checking reference alignment:"
 bam_align=$(samtools view -H $bam_inpath/$name.bam | awk '{if(NR==2) print}' | awk 'match($0, /AS:.*/) {split(substr($0, RSTART, RLENGTH), a, " ")
 print a[1]}' | sed 's/^...//')
 if [ -z "${bam_align}" ] ; then # files doesn't contain "AS:". Try "UR:"
 	bam_align=$(samtools view -H $bam_inpath/$name.bam | awk '{if(NR==2) print}' | awk 'match($0, /UR:.*/) {split(substr($0, RSTART, RLENGTH), a, " ")
 	print a[1]}' | sed 's/^...//')
 fi
 if [ -z "${bam_align}" ] ; then # files doesn't contain "UR:". Try "SN:"
 	bam_align=$(samtools view -H $bam_inpath/$name.bam | awk '{if(NR==2) print}' | awk 'match($0, /SN:.*/) {split(substr($0, RSTART, RLENGTH), a, " ")
 	print a[1]}' | sed 's/^...//')
 	if [[ "$bam_align" =~ "chr" ]] ; then # only hg19 contains "chr"
 		bam_align="hg19"
 	else
 		bam_align="grch37"
 	fi
 fi
 
# find which assembly the bam file says it is
if [[ "$bam_align" =~ .*([Gg][Rr][Cc][Hh]37)+.* ]] ; then
	bam_align="grch37"
elif [[ "$bam_align" =~ .*([Hh][Gg]19)+.* ]] ; then
	bam_align="hg19"
elif [[ "$bam_align" =~ .*(Homo_sapiens_assembly19.fasta)+.* ]] ; then
	bam_align="hg19"
fi
 
echo $bam_align

if [ "$bam_align" == "hg19" ] ; then
	samtools view -b -F 4 $bam_inpath/$name.bam chr7:38295938-38407399 chr7:142000817-142510993 chr9:33618203-33662661 chr14:22090036-23014042 > $bam_outpath/$name.TCR.bam
	wait
fi
if [ "$bam_align" == "grch37" ] ; then
	samtools view -b -F 4 $bam_inpath/$name.bam 7:38295938-38407399 7:142000817-142510993 9:33618203-33662661 14:22090036-23014042 > $bam_outpath/$name.TCR.bam
	wait
fi
wait
	
samtools view -b -f 4 $bam_inpath/$name.bam > $bam_outpath/$name.unmapped.bam
wait
samtools index $bam_outpath/$name.TCR.bam
wait
samtools index $bam_outpath/$name.unmapped.bam
wait
samtools merge $bam_outpath/$name.TCRreg.bam $bam_outpath/$name.TCR.bam $bam_outpath/$name.unmapped.bam
wait
echo 'Overlap-Unmapped-count:' >> $bam_outpath/$name.summary.txt
samtools view -c $bam_outpath/$name.TCRreg.bam >> $bam_outpath/$name.summary.txt # unmapped-overlapping-count
wait
samtools sort -n  $bam_outpath/$name.TCRreg.bam $bam_outpath/$name.TCRreg.sorted
wait
bedtools bamtofastq -i $bam_outpath/$name.TCRreg.sorted.bam -fq $bam_outpath/$name.TCRreg.fastq
wait
bgzip $bam_outpath/$name.TCRreg.fastq
wait


# move back
mv $bam_outpath/$name.TCRreg.fastq.gz $outdir
mv $bam_outpath/$name.summary.txt $outdir
mv $bam_outpath/$name.TCRreg.sorted.bam $outdir
wait

# clean up
rm $bam_outpath/$name*

#rm $bam_outpath/$name.TCRreg.bam.bai
#rm $bam_outpath/$name.unmapped.bam
#rm $bam_outpath/$name.unmapped.bam.bai
#rm $bam_outpath/$name.TCRreg.merged.bam
#rm $bam_outpath/$name.TCRreg.merged.bam.bai
#rm $bam_outpath/$name.TCRreg.merged.sorted.bam
