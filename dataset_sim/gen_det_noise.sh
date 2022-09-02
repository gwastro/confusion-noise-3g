 #!/bin/bash

python det_noise_pycbc.py --detector-name CE \
                          --fs 4096 --flow 5 --duration 6 \
                          --seed 0 --output-folder ./CE_6h
