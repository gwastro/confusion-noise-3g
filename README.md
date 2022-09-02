## A mock data study for 3G ground-based detectors: the performance loss of matched filtering due to the correlated confusion noise

Shichao Wu <sup>1,2</sup> and Alexander H. Nitz <sup>1,2</sup>

<sub>1. [Albert-Einstein-Institut, Max-Planck-Institut for Gravitationsphysik, D-30167 Hannover, Germany](http://www.aei.mpg.de/obs-rel-cos)</sub>  
<sub>2. Leibniz Universitat Hannover, D-30167 Hannover, Germany</sub>

## Introduction

The next-generation (3G/XG) ground-based gravitational-wave (GW) detectors such as Einstein Telescope (ET) and Cosmic Explorer (CE) will begin observing in the next decade. Due to the extremely high sensitivity of these detectors, the majority of stellar-mass compact-binary mergers in the entire Universe will be observed. It is also expected that 3G detectors will have significant sensitivity down to 2-7 Hz; the observed duration of binary neutron star signals could increase to several hours or days. The abundance and duration of signals will cause them to overlap in time, which may form a confusion noise that could affect the detection of individual GW sources when using naive matched filtering; Matched filtering is only optimal for stationary Gaussian noise. We create mock data for CE and ET using the latest population models informed by GWTC-3 catalog and investigate the performance loss of matched filtering due to overlapping signals. We find the performance loss mainly comes from a deviation in the noise's measured amplitude spectral density. The redshift reach of CE (ET) can be reduced by 15-38 (8-21) % if we use the median (upper) merger rate estimate. The direct contribution of confusion noise to the total SNR is generally negligible compared to the contribution from instrumental noise. We also find that correlated confusion noise has a negligible effect on the quadrature summation rule of network SNR for ET, but might reduce the network SNR of high detector-frame mass signals for detector networks including CE if no mitigation is applied. For ET, the null stream can mitigate the astrophysical foreground. For CE, we demonstrate that a computationally efficient, straightforward single-detector signal identification and subtraction method successfully suppresses the total noise to almost the instrument noise level.

## Folders and Scripts

This repository contain the main scripts used in [our paper](). Some key codes or functions called by these scripts have been written into `PyCBC`, such as [`population.population_models`](https://github.com/gwastro/pycbc/blob/master/pycbc/population/population_models.py), [`pycbc_brute_bank`](https://github.com/gwastro/pycbc/blob/master/bin/bank/pycbc_brute_bank), [`distributions.external`](https://github.com/gwastro/pycbc/blob/master/pycbc/distributions/external.py), [`distributions.utils`](https://github.com/gwastro/pycbc/blob/master/pycbc/distributions/utils.py), and [`psd.read`](https://github.com/gwastro/pycbc/blob/master/pycbc/psd/read.py) modules.
 * `dataset_sim`: scripts to generate CE and ET mock datasets.
 * `figure_notebooks`: notebooks to generate all figures in our paper (BE PUBLIC ONCE PAPER IS ACCEPTED).
 * `search`: scripts to search our mock datasets.
 * `subtraction`: scripts to do the first-stage foreground cleaning.
 * `template_banks`: scripts to generate the BNS template banks for CE and ET.

## License and Citation

This repository is licensed under [GNU General Public License v3.0](https://github.com/gwastro/confusion_noise_3g/blob/main/LICENSE).
We encourage use of these scripts in derivative works. If you use the material provided here, please cite the paper using the reference:

```
@article{
}
```