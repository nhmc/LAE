cp find_par.py $1
cp model.cfg $1
cp model.py $1
cp printpar.sh $1
cp run_emcee.py $1
cp run.sh $1
d=$1/comp1
mkdir $d
cp comp1/dont_use $d
cp comp1/priors $d 
cp comp1/observed_logN $d 
cd $d 
ln -s ../model.cfg .
ln -s ../model.py .


