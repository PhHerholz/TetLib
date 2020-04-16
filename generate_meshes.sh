cd build

mkdir -p out

for CSIZE in 0.1 0.05 #0.3 0.1 0.05 #0.092 0.084 0.076 0.068 0.061 0.053 0.045 0.037 0.03
do
for CERATIO in 2
do
for FSIZE in 0.05 0.1 0.25 0.5 0.75 1. 
do
for NFLIPS in 0 #250 500 750 1000
do
for OPTS in "000"
do
	./mesh_generation $CSIZE $CERATIO $FSIZE ${OPTS:0:1} ${OPTS:1:1} ${OPTS:2:1} $NFLIPS 
done
done
done
done
done

cd ..
python3 plot_mesh_quality.py build/out
mv build/out/*.meshfile build/out/plots/
mv build/out/meshes.txt build/out/plots/

FOLDERNAME="00_facetsize"

echo $FOLDERNAME

mkdir build/out/$FOLDERNAME
mv build/out/plots build/out/$FOLDERNAME
mv build/out/*.png build/out/$FOLDERNAME
mv build/out/*.csv build/out/$FOLDERNAME

# ----- HEAT DIFF EXP ---------
cd build
./tet_experiment out/$FOLDERNAME/plots/
