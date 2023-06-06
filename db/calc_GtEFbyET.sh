#!/bin/bash
#------------------------------------------------
#author: Ramiro A. Espada. April, 2023.
#------------------------------------------------
export LC_NUMERIC="en_US.UTF-8"
#============================================
#Inputs:
     EFfile="./orig/EFv210806.csv"           #Emission factor of each VegId
SpecGTfiles=("./orig/SpeciationCrop210806.csv"  "./orig/SpeciationHerb210806.csv"  "./orig/SpeciationShrub210806.csv"  "./orig/SpeciationTree210725.csv")      #"Speciation" files for each GrowthType.
         GT=("crop" "grass" "shrub" "tree") #GrowthType (same order than above!)
	
    outFile='GtEFbyEcotype.csv'
#============================================
#Reescribo los archivos de forma mas facil de trabajar:
echo "Parsing Vegetation EFs and LDFs .. "
if [ -f 'tmp_vegFrac.csv' ] 
then
	rm 'tmp_vegFrac.csv'
fi;
if [ -f 'tmp_specGT.csv' ] 
then
	rm 'tmp_specGT.csv'
fi;

#tmp_vegFrac.csv:
#vegid, EF1,EF2,....,LDF19,LDF1,..,LDF4
awk -F"," 'BEGIN {OFS = ","} NR>1{print $1,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29}' $EFfile  > tmp_vegFrac.csv

echo "Parsing Speciation files .. "
# tmp_specGT.csv:
# gtype,ecotype,vegid,vegfrac
#echo "gtype,ecotype,vegid,vegfrac" > tmp_specGT.csv
for k in $(seq 0 $(( ${#GT[@]}-1)) )
do
        echo $k ${GT[$k]}
        gt=${GT[$k]}
        gtFile=${SpecGTfiles[$k]}
        awk -F"," -v gt=$gt 'BEGIN {OFS=","}NR>1{print gt,$0}' $gtFile >> tmp_specGT.csv
done
#============================================
#Calculo el EF de cada "growth-type" (GT) agrupado por Ecotype
echo " Computing EF for each Growth-Type grouped by Ecotype.. "

if [ -f $outFile ] 
then
	rm $outFile
fi

awk -F ',' '
BEGIN {OFS = ","} 
NR == FNR {
   #VegId, EF1,EF2,....,LDF1,...,LDF4
   #$1     $2, $3,......,...,...,$(n)
   n=NF
   for (i = 2; i <= n; i++){
        EF[$1,i-1] = $i
   };next
} 
NR != FNR {
   # gtype,ecotype,vegid,vegfrac
   #$1     $2      $3     $4
   for (i = 2; i <= n; i++) {
       EFe[$1,$2,i-1] += $4 * EF[$3,i-1]
   }
} 
END {
   for (i in EFe) {                     #Es posible que acá este haciendo más loops de lo necesario.
       split(i,a,SUBSEP);
       printf("%s,%s", a[1],a[2])
       for (j = 1; j < n; j++) {
          printf(",%.3f",EFe[a[1],a[2],j]);
       }
       printf("\n");
   }
}' tmp_vegFrac.csv tmp_specGT.csv | sort -u > $outFile

#Limpio archivos temprales:
rm tmp_vegFrac.csv tmp_specGT.csv

