for indices in {1..2}
do
	for chrs in {1..22}
	do
	  #if [ ! -f `printf "simulation_%03d.csv" $i` ]; then
	    #echo $i
	    bsub -q long_normal -n 1 -R "span[hosts=1]" -o job_%J.out -e job_%J.err -W 24:00 -M 20G "Rscript run_coloc.R $indices $chrs"
	  #fi
	done
done
