for comp in all
#for comp in comp?
do 
    echo "Component $comp"
    cd $comp
    for n in uvb_k??
    do
	echo "    $n"
	cd $n
	#echo "dry run"
	qsub qsub.sh
	cd ..
    done
    cd ..
done
