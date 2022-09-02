 #!/bin/bash

 python confusion_noise.py --signal-type BNS --detector-name CE \
                           --rotation False --duration 6 \
                           --average-time-interval 4.16166232002214 \
                           --seed 0 --input-config ./population_files/prior_files/o1o2o3_lvk_bns_upper_5Hz.ini \
                           --output-folder ./bns_upper_CE_norotation_6h \
                           --thread-number 2
