for i in *
do

obabel -imol $i/*.mol -oxyz --addfilename | obabel -ixyz -ocan | sed 's/^/'"$i"'    /g' >> active.smi

done
