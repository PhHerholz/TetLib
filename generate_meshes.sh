cd build

mkdir -p out

for CVOL in 0.5 0.4 #0.092 0.084 0.076 0.068 0.061 0.053 0.045 0.037 0.03
do
for CERATIO in 2
do
for NFLIPS in 0 #1000
do
for OPTS in "000"
do
	./mesh_generation $CVOL $CERATIO ${OPTS:0:1} ${OPTS:1:1} ${OPTS:2:1} $NFLIPS 
done
done
done
done

cd ..
python3 plot_mesh_quality.py build/out
mv build/out/*.meshfile build/out/plots/
mv build/out/meshes.txt build/out/plots/

FOLDERNAME="run_00"

echo $FOLDERNAME

mkdir build/out/$FOLDERNAME
mv build/out/plots build/out/$FOLDERNAME
mv build/out/*.png build/out/$FOLDERNAME
mv build/out/*.csv build/out/$FOLDERNAME
