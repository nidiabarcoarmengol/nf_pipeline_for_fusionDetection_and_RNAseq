#!/bin/bash

while getopts u:a:f:b:d:x:w:o: flag
do
    case "${flag}" in
        a) Arriba=${OPTARG};;
        f) FusionCatcher=${OPTARG};;
        u) STAR=${OPTARG};;
        b) Sample=${OPTARG};;
        d) Sample_type=${OPTARG};;
        x) Disease_name=${OPTARG};;
        w) Cicero=${OPTARG};;
        o) Output_dir=${OPTARG};;
    esac
done

##### Arriba Fusions

if [[ -n "$Arriba" ]]; then

	lines=$(($(wc -l < "$Arriba") - 1))

	echo "Fusionses de Arriba introducidas"
	
	paste <(tail -n +2 $Arriba | awk '{print $5}' | awk -F':' '{print $1}') \
	<(tail -n +2 $Arriba | awk '{print $5}' | awk -F':' '{print $2}') \
	<(tail -n +2 $Arriba | awk '{print $3}' | awk -F'/' '{print $2}') \
	<(tail -n +2 $Arriba | awk '{print $6}' | awk -F':' '{print $1}') \
	<(tail -n +2 $Arriba | awk '{print $6}' | awk -F':' '{print $2}') \
	<(tail -n +2 $Arriba | awk '{print $4}' | awk -F'/' '{print $2}') \
	<(yes "RNA" | head -n "$lines") \
	<(yes "$Sample" | head -n "$lines") \
	<(yes "$Sample_type" | head -n "$lines") \
	<(yes "$Disease_name" | head -n "$lines") \
	<(yes "Arriba" | head -n "$lines") \
	<(tail -n +2 $Arriba | awk '{print $1}') \
	<(tail -n +2 $Arriba | awk '{print $7}') \
	<(tail -n +2 $Arriba | awk '{print $2}') \
	<(tail -n +2 $Arriba | awk '{print $8}') \
	<(tail -n +2 $Arriba | awk '{print $23}') \
	<(tail -n +2 $Arriba | awk '{print $24}') \
	<(tail -n +2 $Arriba | awk '{print $9}') \
	<(tail -n +2 $Arriba | awk '{print $15}') \
	<(yes "NA" | head -n "$lines") \
	<(tail -n +2 $Arriba | awk '{print $16}') \
	<(yes "NA" | head -n "$lines") \
	<(tail -n +2 $Arriba | awk '{print $28}') \
	<(tail -n +2 $Arriba | awk '{print $17}') \
	<(tail -n +2 $Arriba | awk '{print $10}') \
	<(tail -n +2 $Arriba | awk '{print $11}') \
	<(tail -n +2 $Arriba | awk '{print $11}') > "$Output_dir"/"$Sample"_Arriba_fusions.txt
fi

##### FusionCatcher Fusions

if [[ -n "$FusionCatcher" ]]; then
	
	sed -i 's/\t\t/\tNA\t/g' $FusionCatcher
	lines=$(($(wc -l < "$FusionCatcher") - 1))
	
	echo "Fusionses de FusionCatcher introducidas"
	
	paste <(tail -n +2 $FusionCatcher | awk '{print $9}' | awk -F':' '{print $1}') \
	<(tail -n +2 $FusionCatcher | awk '{print $9}' | awk -F':' '{print $2}') \
	<(tail -n +2 $FusionCatcher | awk '{print $9}' | awk -F':' '{print $3}') \
	<(tail -n +2 $FusionCatcher | awk '{print $10}' | awk -F':' '{print $1}') \
	<(tail -n +2 $FusionCatcher | awk '{print $10}' | awk -F':' '{print $2}') \
	<(tail -n +2 $FusionCatcher | awk '{print $10}' | awk -F':' '{print $3}') \
	<(yes "RNA" | head -n "$lines") \
	<(yes "$Sample" | head -n "$lines") \
	<(yes "$Sample_type" | head -n "$lines") \
	<(yes "$Disease_name" | head -n "$lines") \
	<(yes "FusionCatcher" | head -n "$lines") \
	<(tail -n +2 $FusionCatcher | awk '{print $1}') \
	<(yes "NA" | head -n "$lines") \
	<(tail -n +2 $FusionCatcher | awk '{print $2}') \
	<(yes "NA" | head -n "$lines") \
	<(yes "NA" | head -n "$lines") \
	<(yes "NA" | head -n "$lines") \
	<(yes "NA" | head -n "$lines") \
	<(yes "NA" | head -n "$lines") \
	<(yes "NA" | head -n "$lines") \
	<(tail -n +2 $FusionCatcher | awk '{print $16}') \
	<(tail -n +2 $FusionCatcher | awk '{print $3}') \
	<(tail -n +2 $FusionCatcher | awk '{print $15}') \
	<(yes "NA" | head -n "$lines") \
	<(tail -n +2 $FusionCatcher | awk '{print $6}') \
	<(yes "NA" | head -n "$lines") \
	<(tail -n +2 $FusionCatcher | awk '{print $5}') > "$Output_dir"/"$Sample"_FusionCatcher_fusions.txt
fi

##### STAR Fusions

if [[ -n "$STAR" ]]; then
	
	lines=$(($(wc -l < "$STAR") - 1))
	
	echo "Fusionses de STAR introducidas"
	
	paste <(tail -n +2 $STAR | awk '{print $6}' | awk -F':' '{gsub(/^chr/, "", $1); print $1}') \
	<(tail -n +2 $STAR | awk '{print $6}' | awk -F':' '{print $2}') \
	<(tail -n +2 $STAR | awk '{print $6}' | awk -F':' '{print $3}') \
	<(tail -n +2 $STAR | awk '{print $8}' | awk -F':' '{gsub(/^chr/, "", $1); print $1}') \
	<(tail -n +2 $STAR | awk '{print $8}' | awk -F':' '{print $2}') \
	<(tail -n +2 $STAR | awk '{print $8}' | awk -F':' '{print $3}') \
	<(yes "RNA" | head -n "$lines") \
	<(yes "$Sample" | head -n "$lines") \
	<(yes "$Sample_type" | head -n "$lines") \
	<(yes "$Disease_name" | head -n "$lines") \
	<(yes "STAR_Fusion" | head -n "$lines") \
	<(tail -n +2 $STAR | awk '{print $1}' | awk -F'--' '{print $1}') \
	<(yes "NA" | head -n "$lines") \
	<(tail -n +2 $STAR | awk '{print $1}' | awk -F'--' '{print $2}') \
	<(yes "NA" | head -n "$lines") \
	<(yes "NA" | head -n "$lines") \
	<(yes "NA" | head -n "$lines") \
	<(yes "NA" | head -n "$lines") \
	<(yes "NA" | head -n "$lines") \
	<(yes "NA" | head -n "$lines") \
	<(yes "NA" | head -n "$lines") \
	<(yes "NA" | head -n "$lines") \
	<(yes "NA" | head -n "$lines") \
	<(tail -n +2 $STAR | awk '{print $17}') \
	<(tail -n +2 $STAR | awk '{print $2}') \
	<(yes "NA" | head -n "$lines") \
	<(tail -n +2 $STAR | awk '{print $3}') > "$Output_dir"/"$Sample"_STAR_fusions.txt
fi

##### CICERO Fusions

if [[ -n "$Cicero" ]]; then
	
	lines=$(($(wc -l < "$Cicero") - 1))
	
	echo "Fusionses de CICERO introducidas"
	
	paste <(tail -n +2 $Cicero | awk '{print $3}' | awk -F':' '{gsub(/^chr/, "", $1); print $1}') \
	<(tail -n +2 $Cicero | awk '{print $4}') \
	<(tail -n +2 $Cicero | awk '{print $5}') \
	<(tail -n +2 $Cicero | awk '{print $8}' | awk -F':' '{gsub(/^chr/, "", $1); print $1}') \
	<(tail -n +2 $Cicero | awk '{print $9}') \
	<(tail -n +2 $Cicero | awk '{print $10}') \
	<(yes "RNA" | head -n "$lines") \
	<(yes "$Sample" | head -n "$lines") \
	<(yes "$Sample_type" | head -n "$lines") \
	<(yes "$Disease_name" | head -n "$lines") \
	<(yes "CICERO" | head -n "$lines") \
	<(tail -n +2 $Cicero | awk '{print $2}') \
	<(tail -n +2 $Cicero | awk '{print $6}') \
	<(tail -n +2 $Cicero | awk '{print $7}') \
	<(tail -n +2 $Cicero | awk '{print $11}') \
	<(yes "NA" | head -n "$lines") \
	<(yes "NA" | head -n "$lines") \
	<(tail -n +2 $Cicero | awk '{print $32}') \
	<(tail -n +2 $Cicero | awk '{print $30}') \
	<(tail -n +2 $Cicero | awk '{print $31}') \
	<(yes "NA" | head -n "$lines") \
	<(yes "NA" | head -n "$lines") \
	<(tail -n +2 $Cicero | awk '{print $27}') \
	<(yes "NA" | head -n "$lines") \
	<(tail -n +2 $Cicero | awk '{print $13}') \
	<(tail -n +2 $Cicero | awk '{print $14}') \
	<(yes "NA" | head -n "$lines") > "$Output_dir"/"$Sample"_CICERO_fusions.txt
fi

cat "$Output_dir"/"$Sample"_Arriba_fusions.txt "$Output_dir"/"$Sample"_FusionCatcher_fusions.txt "$Output_dir"/"$Sample"_STAR_fusions.txt "$Output_dir"/"$Sample"_CICERO_fusions.txt > "$Output_dir"/Fusions.txt

if [[ -n "$Arriba" ]]; then
   awk -F'\t' '{print $1,$2,$3,$12,$4,$5,$6,$14}' "$Output_dir"/"$Sample"_Arriba_fusions.txt > "$Output_dir"/Arriba_fusions_2.txt
fi
if [[ -n "$FusionCatcher" ]]; then
   awk -F'\t' '{print $1,$2,$3,$12,$4,$5,$6,$14}' "$Output_dir"/"$Sample"_FusionCatcher_fusions.txt > "$Output_dir"/FusionCatcher_fusions_2.txt
fi
if [[ -n "$STAR" ]]; then
   awk -F'\t' '{print $1,$2,$3,$12,$4,$5,$6,$14}' "$Output_dir"/"$Sample"_STAR_fusions.txt > "$Output_dir"/STAR_fusions_2.txt
fi

if [[ -n "$Cicero" ]]; then
   awk -F'\t' '{print $1,$2,$3,$12,$4,$5,$6,$14}' "$Output_dir"/"$Sample"_CICERO_fusions.txt > "$Output_dir"/CICERO_fusions_2.txt
fi

echo -e "Chr1\tPos1\tStrand1\tChr2\tPos2\tStrand2\tRNA\tSample\tSample_type\tDisease_name\tTool\tGene1\tGene_location1\tGene2\tGene_location2\tTranscript_id1\tTranscript_id2\tType\tRating\tMedal\tReading_frame\tFusion_description\tFusion_sequence\tAnnotation\tSplit_reads1\tSplit_reads2\tSpanning_pairs" > "$Output_dir"/"$Sample"_Fusions_all.csv
sed -i 's/ /\t/g' "$Output_dir"/Fusions.txt
cat "$Output_dir"/Fusions.txt >> "$Output_dir"/"$Sample"_Fusions_all.csv

awk '!seen[$0]++' "$Output_dir"/Fusions.txt > "$Output_dir"/Unique_fusions.txt

awk -F'\t' '{print $1,$2,$3,$12,$4,$5,$6,$14}' "$Output_dir"/Unique_fusions.txt > "$Output_dir"/Unique_fusions_filtered.txt

echo -e "Chr1\tPos1\tStrand1\tGene1\tChr2\tPos2\tStrand2\tGene2\tNcaller\tSet" > "$Output_dir"/Fusions_set_no_filt.csv

while IFS= read -r line 
do
    count=0
    fusion_tools=()
     if [[ -n "$Arriba" ]]; then
    	if grep -F "$line" "$Output_dir"/Arriba_fusions_2.txt > /dev/null; then
    		fusion_tools+=("Arriba")
    		((count++))
    	fi
    fi
    
    if [[ -n "$FusionCatcher" ]]; then
    	if grep -F "$line" "$Output_dir"/FusionCatcher_fusions_2.txt > /dev/null; then
    		fusion_tools+=("FusionCatcher")
    		((count++))
    	fi
    fi
    
    if [[ -n "$STAR" ]]; then
        	if grep -F "$line" "$Output_dir"/STAR_fusions_2.txt > /dev/null; then
    		fusion_tools+=("STAR_fusions")
    		((count++))
    	fi
    fi
    
    if [[ -n "$Cicero" ]]; then
    	if grep -F "$line" "$Output_dir"/CICERO_fusions_2.txt > /dev/null; then
    		fusion_tools+=("Cicero")
    		((count++))
    	fi
    fi
   
    fusion_tools_string=$(IFS=/; echo "${fusion_tools[*]}")

    paste <(echo $line) <(echo $count) <(echo $fusion_tools_string) > "$Output_dir"/temp.txt
    sed -i 's/ /\t/g' "$Output_dir"/temp.txt
    
    cat "$Output_dir"/temp.txt >> "$Output_dir"/Fusions_set_no_filt.csv
done < "$Output_dir"/Unique_fusions_filtered.txt

awk '!seen[$0]++' "$Output_dir"/Fusions_set_no_filt.csv > "$Output_dir"/"$Sample"_Fusions_set.csv

if [[ -n "$Arriba" ]]; then
rm "$Output_dir"/Arriba_fusions_2.txt
fi

if [[ -n "$FusionCatcher" ]]; then
rm "$Output_dir"/FusionCatcher_fusions_2.txt
fi

if [[ -n "$STAR" ]]; then
rm "$Output_dir"/STAR_fusions_2.txt
fi

if [[ -n "$Cicero" ]]; then
rm "$Output_dir"/CICERO_fusions_2.txt
fi

rm "$Output_dir"/Fusions.txt
rm "$Output_dir"/Unique_fusions.txt
rm "$Output_dir"/Unique_fusions_filtered.txt
rm "$Output_dir"/Fusions_set_no_filt.csv
rm "$Output_dir"/temp.txt

