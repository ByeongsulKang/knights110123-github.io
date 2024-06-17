for i in {1..100}; do
       obabel -ipdb FMM2.pdb -omol --gen3d -O "FMM2_$[i].mol"
done       
