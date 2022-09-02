 #!/bin/bash

 python confusion_noise.py --signal-type BBH --detector-name CE \
                           --rotation False --duration 6 \
                           --average-time-interval 359.4162912746394 \
                           --seed 0 --input-config ./population_files/prior_files/o1o2o3_lvk_bbh_median_5Hz.ini \
                           --output-folder ./bbh_median_CE_norotation_6h \
                           --thread-number 1
