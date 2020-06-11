
cd build
mkdir -p out

for CSIZE in 0.1 0.092 0.084 0.076 0.068 0.061 0.053 0.045 0.037 0.03
do

FSIZE=$CSIZE

	./mesh_generation $CSIZE 0
done

cd ..
mkdir -p build/out/plots

mv build/out/*.meshfile build/out/plots/
mv build/out/meshes.txt build/out/plots/
# mv build/out/origin_moves.csv build/out/plots/

FOLDERNAME="11_doublesphere" #_r${REP}"

echo $FOLDERNAME

mkdir build/out/$FOLDERNAME
mv build/out/plots build/out/$FOLDERNAME
mv build/out/*.png build/out/$FOLDERNAME
mv build/out/*.csv build/out/$FOLDERNAME

# ----- HEAT DIFF EXP ---------
cd build
# ./tet_experiment out/$FOLDERNAME/plots/
#cd ..
#rm build/out/$FOLDERNAME/plots/*.meshfile
