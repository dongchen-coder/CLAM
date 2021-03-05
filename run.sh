for i in 2 4 8 128
do
    make run_all WAYS=$i
    mv ./leases/*_fnsl_leases ./leases/poly_all/128blocks_${i}ways
    mv ./leases/*_fsl_leases ./leases/poly_all/128blocks_${i}ways_phases
done
