cd comp1/final
python ../../find_par.py --header
cd ../..

for n in 2 3 4 5 6 7 8 
do
    cd comp$n/final
    python ../../find_par.py
    cd ../..
done

