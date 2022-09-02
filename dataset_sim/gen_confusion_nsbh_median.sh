 #!/bin/bash

 python confusion_noise.py --signal-type NSBH --detector-name ET-D \
                           --rotation False --duration 6 \
                           --average-time-interval 46.512696517894504 \
                           --seed 0 --input-config ./population_files/prior_files/o1o2o3_lvk_nsbh_median_2Hz.ini \
                           --output-folder ./nsbh_median_ETD_norotation_6h \
                           --thread-number 2
