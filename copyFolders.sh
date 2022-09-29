mkdir Results
for i in {1..10}
do
for j in M1_$i M2_$i M1_05_$i M1D_$i M2D_$i M1R_$i M2R_$i Ped_$i M1D_Est_$i
do
echo $j
cp ${j}/${j}_info.txt ${j}/${j}_p_markers.txt ${j}/${j}_p_neutral.txt ${j}/${j}_p_qtl.txt ${j}/${j}_qtl_effects.txt Results
done
done
