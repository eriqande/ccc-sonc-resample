fi=simCoNe.log; rm -f  $fi; for (( i=10; i<=400; i+=2 )); do echo "Starting simulations for n0 = $i on $(date)" >> $fi; simCoNeprob $i 1000000 > ${i}pr.txt; echo "Done with simulations for $i on $(date)" >> $fi; echo >> $fi;  done 