#!/bin/sh
## Script to run for variant quality, variant annotation filtering
###########################################################################################################################

annotationResults=$1   
genotypingVariants=$2   
parameter=$3    
outdir=$4

HOME="/rfs-storageservice/GIS/Projects/HG6/mmlian/PD_Intl_WES/manuscript_scripts/annot"

if [ $# -lt 4 ]
then
        echo "Usage:"
        echo "	filter_annotation.pl [arguments] <variantannot_samplegeno.tsv> <samplegeno.tsv> <parameters.txt> <outputFolder/>"
        echo "	Function  : To do snp filtering based on user-defined criteria for annotated variant list"
        echo "	Example   : perl filter_annotation.pl variantannot_samplegeno.tsv samplegeno.tsv v2.3_ns_MAF1_del_can_params.txt v2.3_ns_MAF1_del_can/"
	echo "	Argument1 : Annotated Main TSV file"
        echo "	Argument2 : Genotypes file"
	echo "	Argument3 : Parameter/criteria file"
	echo "	Argument4 : OutputFolder ended with slash"
	echo ""
        exit
fi


isFileEmpty()
{
	if [[ ! -s $1 ]]; then
		echo "$1 file is empty"
		exit 1;
	fi
 }
isFileExist()
{
	if [[ ! -f $1 ]]; then
		echo "Cannot open $1"
		exit 1;
	fi;
}
isDirExist()
{
	if [[ ! -d $1 ]]; then
		 mkdir "$outdir" 
	fi;
}
isFileExist $annotationResults
isFileEmpty $annotationResults
isFileExist $genotypingVariants
isFileEmpty $genotypingVariants
isFileExist $parameter
isFileEmpty $parameter
isDirExist $outdir

count=`grep "####" $parameter |wc -l`
#echo "count:$count"
head=`head -1 $genotypingVariants |  awk -v OFS="\t" '$1=$1' | sed -e 's/:1//g;s/:0//g'`
#IFS='  ' 
IFS=$'\t\n:'
g_head=($head)

intSize=${#g_head[@]};
if [[ $count == 0 ]] ; then
	arr=$parameter
fi
if [[ $count == 1 ]] ; then
	csplit -f $outdir\/param_split $parameter "/^##/"
	arr=`ls -rt $outdir\/param_split*`
fi
if [[ $count > 1 ]] ; then
	let "ss=$count - 1"
	csplit -f $outdir\/param_split $parameter "/^##/" {$ss}
	arr=`ls -rt $outdir\/param_split*`
fi
unset IFS

IFS=$'\n'
param_array=($arr)

#read -a param_array <<<$arr

for (( w=0; w <= $count; w++ )); do 
	#echo "param array:${param_array[$w]}"
	sed -i '/####/d' ${param_array[$w]}

	unset IFS

	perl $HOME\/filtering_annotations.pl $annotationResults $genotypingVariants ${param_array[$w]} $outdir$w
	head -n1 $annotationResults > $outdir\/header_annotation_results
	head -n1 $genotypingVariants > $outdir\/header_genotyping_results
done


cat $outdir\/*final_all_results_filtered* > $outdir\/temp_final_tsv_results_filtered
cat $outdir\/*samplegeno_filtered* > $outdir\/temp_seq_results_filtered

rm $outdir\/*tsv
cat $outdir\/header_annotation_results $outdir\/temp_final_tsv_results_filtered |awk 'NR==1; NR > 1 {print $0 | "sort -k1,1 -k2,2n"}' |uniq > $outdir\/variantannot_samplegeno_final_all_results_filtered.tsv
cat $outdir\/header_genotyping_results $outdir\/temp_seq_results_filtered |awk 'NR==1; NR > 1 {print $0 | "sort -k1,1 -k2,2n"}' |uniq > $outdir\/samplegeno_filtered.tsv

rm $outdir\/header_annotation_results
rm $outdir\/header_genotyping_results
rm $outdir\/temp_final_tsv_results_filtered
rm $outdir\/temp_seq_results_filtered

