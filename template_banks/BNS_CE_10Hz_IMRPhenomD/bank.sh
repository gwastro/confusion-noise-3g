 #!/bin/bash


if [ $# -eq 0 ]
then
	echo "Provide the psd file"
	exit 1
else
	echo "Computing aligned bank using psd $1"


OMP_NUM_THREADS=1 python -m cProfile -o log /home/anaconda3/envs/pycbc_dev/bin/pycbc_brute_bank \
--verbose \
--output-file nospinbank.hdf \
--minimal-match 0.97 \
--tolerance .001 \
--buffer-length 4 \
--sample-rate 4096 \
--tau0-threshold 0.5 \
--approximant IMRPhenomD \
--tau0-crawl 10 \
--tau0-start 0 \
--tau0-end 400 \
--asd-file $1 \
--input-config CE_bns_bank.ini \
--seed 1 \
--low-frequency-cutoff 10.0

fi