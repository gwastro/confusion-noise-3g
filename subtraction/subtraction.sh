#!/bin/bash

OMP_NUM_THREADS=1 condor_run -a accounting_group=cbc.prod.search -a request_memory=80000 python subtract_triggers_waveform.py \
--confusion-noise-path /mnt/d/confusion_noise_3g/dataset_sim/bns_median_CE_norotation_6h/confusion_noise_CE_1_BNS_21600s.gwf \
--triggers-path /mnt/d/confusion_noise_3g/search/results_bns_CE_medianR0nospin_norotation_snr6p0/ \
--psd-path /mnt/d/confusion_noise_3g/search/data_segments_bns_CE_medianR0nospin_norotation/psd_0.00s_512.00s.txt \
--snr-threshold 10.0 \
--f-low 5.0 \
--time-window 1 \
--output-path ./subtraction_CE_medianR0_snr10p0 >  job_subtraction_CE_medianR0_snr10p0.txt &
