 #!/bin/bash

 python confusion_noise.py --signal-type BBH --detector-name ET-D \
                           --rotation False --duration 6 \
                           --average-time-interval 175.71463128982361 \
                           --seed 0 --input-config ./population_files/prior_files/o1o2o3_lvk_bbh_upper_2Hz.ini \
                           --output-folder ./bbh_upper_ETD_norotation_6h \
                           --thread-number 1
