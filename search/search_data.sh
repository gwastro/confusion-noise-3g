#!/bin/bash


for i in {0..500}; do

    sleep 0.1
#    step=191
    step=305
    start=`echo "$step * ($i + 0)" | bc`
    end=`echo "$step * ($i + 1) - 1" | bc`
    echo $i $start $end

    OMP_NUM_THREADS=1 condor_run -a accounting_group=cbc.prod.search -a request_memory=80000 python find_triggers.py \
    --confusion-noise-path /mnt/d/confusion_noise_3g/dataset_sim/bns_median_ETD_norotation_6h/confusion_noise_ET-D_1_BNS_21600s.gwf \
    --det-noise-path /mnt/d/confusion_noise_3g/dataset_sim/ETD_6h/det_noise_ET-D_21600s.gwf \
    --bank-path /mnt/d/confusion_noise_3g/template_banks/BNS_ETD_10Hz_IMRPhenomD/nospinbank.hdf \
    --template-id-start $start \
    --template-id-end $end \
    --snr-threshold 6.0 \
    --fmin 10 \
    --save-fig False \
    --output-folder .  >  job_$i-$start-$end.txt &

done
