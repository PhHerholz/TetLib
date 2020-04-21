ADD_PLOTS=false
ADD_ORIGIN=true

if $ADD_PLOTS
  then
	plot_flag=1
  else
	plot_flag=0
fi
if $ADD_ORIGIN
  then
	origin_flag=1
  else
	origin_flag=0
fi

cd build
mkdir -p out

for CSIZE in 0.1 0.03 #0.092 0.084 0.076 0.068 0.061 0.053 0.045 0.037 0.03

	#0.06 #0.1 0.0975862068966 0.0951724137931 0.0927586206897 0.0903448275862 0.0879310344828 0.0855172413793 0.0831034482759 0.0806896551724 0.078275862069 0.0758620689655 0.0734482758621 0.0710344827586 0.0686206896552 0.0662068965517 0.0637931034483 0.0613793103448 0.0589655172414 0.0565517241379 0.0541379310345 0.051724137931 0.0493103448276 0.0468965517241 0.0444827586207 0.0420689655172 0.0396551724138 0.0372413793103 0.0348275862069 0.0324137931034 0.03
do
for NORBITPOINTS in 0 #in 0 20 40 60 80 100 120 140 160 180 200
do
for CERATIO in 2
do
#for FSIZE in 0.1 #0.05 0.1 0.25 0.5 0.75 1. 
#do
for VAR in 0 0.2 0.5 # 0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 #250 500 750 1000 1250 1500 1750 2000 2250 2500 2750 3000 3250 3500 3750 4000 # 150 500 750 1000 1250 1500 1750 2000 2250 2500 2750 3000 5000 10000  #250 500 750 1000
do
for OPTS in "000"
do

FSIZE=$CSIZE

	./mesh_generation $CSIZE $CERATIO $FSIZE ${OPTS:0:1} ${OPTS:1:1} ${OPTS:2:1} $VAR $plot_flag $origin_flag
done
done
done
done
done
#done

cd ..
mkdir -p build/out/plots

if $ADD_PLOTS
then
	python3 plot_mesh_quality.py build/out
fi

mv build/out/*.meshfile build/out/plots/
mv build/out/meshes.txt build/out/plots/
mv build/out/origin_moves.csv build/out/plots/

FOLDERNAME="17r2_withorigin" #_r${REP}"

echo $FOLDERNAME

mkdir build/out/$FOLDERNAME
mv build/out/plots build/out/$FOLDERNAME
mv build/out/*.png build/out/$FOLDERNAME
mv build/out/*.csv build/out/$FOLDERNAME

# ----- HEAT DIFF EXP ---------
#cd build
#./tet_experiment out/$FOLDERNAME/plots/
#cd ..
#rm build/out/$FOLDERNAME/plots/*.meshfile
