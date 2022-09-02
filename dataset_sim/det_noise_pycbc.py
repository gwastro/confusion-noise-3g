# Copyright (C) 2022 Shichao Wu
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; If not, see <http://www.gnu.org/licenses/>.


import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import lal
import pycbc.psd
from pycbc import init_logging
from pycbc.frame import write_frame
from pycbc.types.timeseries import TimeSeries
from pycbc.noise.reproduceable import colored_noise
from tqdm import tqdm
import logging


# setup log  # TODO: Not work yet.
logging_level = logging.INFO
logging.basicConfig(format='%(asctime)s : %(message)s', level=logging_level)

parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('--detector-name', required=True,
                    help="The name of the detector to be injected, "
                    "choose from 'ET-D', 'CE'.")
parser.add_argument('--fs', required=True, default=4096,
                    help="The sampling frequency of simulated data. In the unit of Hz.")
parser.add_argument('--flow', required=True, default=5,
                    help="The lower cut-off frequency of simulated data. In the unit of Hz.")
parser.add_argument('--duration', type=float, required=True,
                    help="The duration of simulated data. In the unit of hour.")
parser.add_argument('--seed', type=int, default=0,
                    help="The seed to be used for the random number generator. Default is 0.")
parser.add_argument('--output-folder', required=True,
                    help="The output folder to save simulated noise and paramters file.")
parser.add_argument("--force", action="store_true", default=False,
                    help="If the output-file already exists, overwrite it. "
                         "Otherwise, an OSError is raised.")
parser.add_argument("--verbose", action="store_true", default=False,
                    help="Print logging messages.")

opts = parser.parse_args()
init_logging(opts.verbose)

if os.path.exists(opts.output_folder) and not opts.force:
    raise OSError("output-file already exists; use --force if you wish to "
                  "overwrite it.")
os.makedirs(opts.output_folder)

detector = opts.detector_name
fs = float(opts.fs)
f_low = float(opts.flow)
duration = 3600*opts.duration #  In seconds.
random_seed = opts.seed
output_path = opts.output_folder

ASD_DIR = "./design_asd" 

# Mapping the input detector's name to the detector ASD file name.
det_name_map = {'ET-D': 'ET_D', 'CE': 'CE'}

delta_f = 1.0/8
length = int(fs/2/delta_f) + 1
filename = "%s/%s_asd.txt" % (ASD_DIR, det_name_map[detector])
det_psd = pycbc.psd.from_txt(filename, length, delta_f, f_low, is_asd_file=True)


# Generate stationary and gaussian noise.
print("Simulation start...")
detector_noise = colored_noise(psd=det_psd, start_time=0, end_time=duration,
                               seed=random_seed, sample_rate=fs,
                               low_frequency_cutoff=f_low,
                               filter_duration=128)

# Save the detector noise to the .gwf file.
print("Writing the simulated noise into the .gwf file...")
epoch = lal.LIGOTimeGPS(detector_noise.sample_times[0])
data_noise = TimeSeries(np.array(detector_noise), delta_t=1.0/fs, epoch=epoch)
write_frame("%s/det_noise_%s_%ds.gwf" % (output_path, detector, duration), 
            "H1:LDAS-STRAIN", data_noise)

# Plot the detector noise.
print("Writing is finished, start to plot the timeseries and ASD figure...")
plt.figure(dpi=150)
plt.plot(detector_noise.sample_times, detector_noise)
plt.xlabel("Time (s)")
plt.ylabel("Strain")
plt.title("det_noise_%s_%ds" % (detector, duration))
plt.savefig("%s/det_noise_%s_%ds.png" % (output_path, detector, duration))
plt.close()

# Estimate the PSD as a check of the simulated noise. NOT for the production usage.
## We choose 16 seconds samples that are overlapped 50 percent.
delta_t = 1.0/fs
seg_len = int(16/delta_t)
seg_stride = int(seg_len/2)
estimated_psd = pycbc.psd.welch(detector_noise, seg_len=seg_len, seg_stride=seg_stride)

plt.figure(dpi=150)
plt.loglog(estimated_psd.sample_frequencies, np.sqrt(estimated_psd), label="%s ASD (Welch)" % detector)
plt.loglog(det_psd.sample_frequencies, np.sqrt(det_psd), label="%s ASD (Design)" % detector)
plt.xlim(f_low, fs/2)
plt.ylim(10**(-25), )
plt.xlabel("Frequency (Hz)")
plt.ylabel(r"Strain noise [1/$\sqrt{\mathrm{Hz}}$]")
plt.legend()
plt.savefig("%s/ASD_check_%s_%ds.png" % (output_path, detector, duration))
plt.close()

print("Finished.")
