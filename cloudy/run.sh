for n in 1 2 3 4 5 6 7 8
do
    cd comp$n
    #mv fig/model.pdf final/
    #rm -rf fig
    #cp final/model.png ~/Projects/MPIA_QSO_LBG/LAE/letter/emulateapj/model$n.png
    #cp final/posterior_mcmc.png ~/Projects/MPIA_QSO_LBG/LAE/letter/emulateapj/posterior$n.png
    #pngquant -f final/par2.png
    #cp final/par2-fs8.png ~/Projects/MPIA_QSO_LBG/LAE/letter/emulateapj/pos$n.png
    cp final/model.pdf ~/Projects/MPIA_QSO_LBG/LAE/letter/emulateapj/model$n.pdf
    cd ..
done

# for n in 1 2 3 4 5 6 7 8
# do
#     cd comp$n
#     plot_mcmc final/samples_mcmc.sav.gz
#     mv fig/* final/
#     rmdir fig
#     cd ..
# done


# for n in comp[2345678]
# do
#     cd $n
#     python ../run_emcee.py final
#     cd ..
# done

# for n in comp[12345678]
# do
#     cd $n/final
#     python ../../find_par.py
#     cd ../..
# done


# for n in comp[345678]
# do
#     cd $n
#     rm observed_logN
#     ln -s ../../J0004/$n/observed_logN .
#     cd ..
# done
