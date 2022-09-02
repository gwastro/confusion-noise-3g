## A mock data study for 3G ground-based detectors: the performance loss of matched filtering due to the correlated confusion noise

Shichao Wu <sup>1,2</sup> and Alexander H. Nitz <sup>1,2</sup>

<sub>1. [Albert-Einstein-Institut, Max-Planck-Institut for Gravitationsphysik, D-30167 Hannover, Germany](http://www.aei.mpg.de/obs-rel-cos)</sub>  
<sub>2. Leibniz Universitat Hannover, D-30167 Hannover, Germany</sub>

## Introduction

In the next decade, the next-generation (3G/XG) ground-based GW detectors such as Einstein Telescope and Cosmic Explorer will begin their observation. Due to the extremely high sensitivity of these nearly all-sky detectors, they will detect the majority of stellar-mass compact-binary mergers in the entire Universe. It is expected that 3G detectors will have significant sensitivity down to 2-7 Hz; the observed duration of binary neutron star signals could increase to several hours or days. The abundance and duration of signals will cause them to overlap in time, forming a confusion noise, it might affect the detection of an individual GW signal by using matched filtering, which is only optimal in stationary and Gaussian noise. In this paper, we implement a mock data search for Cosmic Explorer and Einstein Telescope using the latest population models, and we investigate the performance loss of matched filtering due to overlapping signals. We find that the performance loss mainly comes from a deviation in the noise's measured amplitude spectral density. Depending on the merger rate, the reach in the redshift of Cosmic Explorer can be reduced by 15-38\%, or up to 8-21\% for Einstein Telescope. The contribution of confusion noise to the total SNR is generally negligible compared to the contribution from instrumental noise. Furthermore, the correlated confusion noise among the detector network has a negligible effect on the quadrature summation rule of the network SNR for Einstein Telescope, and might reduce the network SNR of high detector-frame mass signals when including Cosmic Explorer. We demonstrate that a straightforward signal identification and subtraction method can suppress the total noise almost to the instrument noise level for Cosmic Explorer; we can expect Einstein Telescope to use the null stream method to achieve a similar result.

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