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
from pycbc import frame
from pycbc.detector import Detector
from pycbc.filter import matched_filter
from pycbc.filter.matchedfilter import sigma
from pycbc.psd import from_txt
from pycbc.types.timeseries import TimeSeries
from pycbc.types.frequencyseries import load_frequencyseries
from pycbc.waveform import get_td_waveform
from pycbc.events.coinc import cluster_over_time
from tqdm import tqdm
from pycbc import init_logging
import logging



# setup log  # TODO: Not work yet.
logging_level = logging.INFO
logging.basicConfig(format='%(asctime)s : %(message)s', level=logging_level)

parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('--confusion-noise-path', required=True,
                    help="The confusion noise to be loaded.")
parser.add_argument('--triggers-path', required=False,
                    help="The triggers which found by 'find_triggers.py'.")
parser.add_argument('--psd-path', required=False,
                    help="The pre-generated PSD.")
parser.add_argument('--snr-threshold', required=True,
                    help="The SNR threshold to be used.")
parser.add_argument('--f-low', required=True,
                    help="The low frequency cutoff to be used.")
parser.add_argument('--f-final', required=False, default=2048,
                    help="The final frequency cutoff to be used in the 2nd waveform.")
parser.add_argument('--time-window', required=True,
                    help="The length of time window (in seconds) to cluster triggers.")
parser.add_argument('--output-path', required=True,
                    help="The output folder to save triggers' files.")
parser.add_argument("--force", action="store_true", default=False,
                    help="If the output-file already exists, overwrite it. "
                         "Otherwise, an OSError is raised.")
parser.add_argument("--verbose", action="store_true", default=False,
                    help="Print logging messages.")

opts = parser.parse_args()
init_logging(opts.verbose)

confusion_noise_path = opts.confusion_noise_path
triggers_path = opts.triggers_path
psd_path = opts.psd_path
snr_threshold = np.float64(opts.snr_threshold)
f_low = np.float64(opts.f_low)
f_final = np.float64(opts.f_final)
time_window = np.float64(opts.time_window)
output_path = opts.output_path

os.makedirs(output_path)
confusion_noise = frame.read_frame(location=confusion_noise_path, channels='H1:LDAS-STRAIN')

triggers_dirs = os.listdir(triggers_path)
if os.path.isfile("%s/psd_0.00s_3600.00s.txt" % triggers_path):
    triggers_dirs.remove("psd_0.00s_3600.00s.txt")

print("Loading all triggers_info files...")
for index in tqdm(range(len(triggers_dirs))):
    if index == 0:
        triggers_info = np.loadtxt("%s/%s/triggers_info.txt" % (triggers_path, triggers_dirs[index]), skiprows=1, unpack=True)
    else:
        triggers_info_new = np.loadtxt("%s/%s/triggers_info.txt" % (triggers_path, triggers_dirs[index]), skiprows=1, unpack=True)
        if np.shape(triggers_info_new) == (12,):
            triggers_info_new = np.expand_dims(triggers_info_new, axis=1)
        triggers_info = np.concatenate((triggers_info, triggers_info_new), axis=1)

triggers_time = triggers_info[0]
triggers_win_start = triggers_info[1]
triggers_win_end = triggers_info[2]
template_masse1 = triggers_info[3]
template_masse2 = triggers_info[4]
triggers_tot_snr = triggers_info[6] + 1j*triggers_info[7]
triggers_op_snr = triggers_info[8] + 1j*triggers_info[9]
triggers_noise_snr = triggers_info[10] + 1j*triggers_info[11]
reduced_triggers_index = []
fs = confusion_noise.sample_rate
flow_1 = 10
flow_2 = f_low

epoch = lal.LIGOTimeGPS(confusion_noise.start_time)
triggers_waveform = TimeSeries(np.zeros(int(confusion_noise.duration*fs+0.5)), delta_t=1.0/fs, epoch=epoch)
reduced_triggers_file = open("%s/reduced_triggers_info.txt" % output_path, 'w')
reduced_triggers_file.write("# time\t mass1\t mass2\t snr_tot_real\t snr_tot_imag\t snr_op_real\t snr_op_imag\t snr_noise_real\t snr_noise_imag\n")

print("Selecting triggers above the given SNR threshold...")
l = np.abs(triggers_tot_snr) >= snr_threshold
print("Clustering all triggers within the given time window, and find the true trigger...")
reduced_triggers_index = cluster_over_time(stat=np.abs(triggers_tot_snr)[l], time=triggers_time[l], window=time_window)

print("Subtracting triggers' rescaled templates...")
for i in tqdm(reduced_triggers_index):
    reduced_triggers_time = triggers_time[l][i]
    reduced_triggers_snr = triggers_tot_snr[l][i]
    reduced_triggers_mass1 = template_masse1[l][i]
    reduced_triggers_mass2 = template_masse2[l][i]
    print("tc: %f tot_SNR: %f mass: (%f, %f)" % (reduced_triggers_time, np.abs(reduced_triggers_snr), reduced_triggers_mass1, reduced_triggers_mass2))
    print("confusion_noise_SNR: %f det_noise_SNR: %f" % (np.abs(triggers_op_snr[l][i]), np.abs(triggers_noise_snr[l][i])))

    hp_10Hz, _ = get_td_waveform(approximant="IMRPhenomD", 
                                 mass1=reduced_triggers_mass1, mass2=reduced_triggers_mass2,
                                 distance=440, f_lower=flow_1, f_ref=flow_1, delta_t=1.0/fs)
    hp_10Hz.start_time += reduced_triggers_time
    template_aligned_10Hz = hp_10Hz
    psd = from_txt(psd_path, length=len(template_aligned_10Hz.data)//2+1, 
                   delta_f=template_aligned_10Hz.delta_f, low_freq_cutoff=flow_1, is_asd_file=False)
    optimal_snr_10Hz = sigma(template_aligned_10Hz, psd=psd, low_frequency_cutoff=flow_1)

    hp_5Hz, _ = get_td_waveform(approximant="IMRPhenomD", 
                                mass1=reduced_triggers_mass1, mass2=reduced_triggers_mass2,
                                distance=440, f_lower=flow_2, f_ref=flow_1, delta_t=1.0/fs, f_final=f_final)
    hp_5Hz.start_time += reduced_triggers_time
    template_aligned_5Hz = hp_5Hz/optimal_snr_10Hz
    template_aligned_5Hz_rescaled = (template_aligned_5Hz.to_frequencyseries()*reduced_triggers_snr).to_timeseries()

    # Combining all the triggers' rescaled templates.
    triggers_waveform = triggers_waveform.add_into(template_aligned_5Hz_rescaled)

    reduced_triggers_file.write("%.18f\t %.18f\t %.18f\t %.8f\t %.8f\t %.8f\t %.8f\t %.8f\t %.8f\n" % 
                                (reduced_triggers_time, reduced_triggers_mass1, reduced_triggers_mass2,
                                np.real(triggers_tot_snr[l][i]), np.imag(triggers_tot_snr[l][i]),
                                np.real(triggers_op_snr[l][i]), np.imag(triggers_op_snr[l][i]),
                                np.real(triggers_noise_snr[l][i]), np.imag(triggers_noise_snr[l][i])))

reduced_triggers_file.close()

print("Saving data...")
# Save triggers waveform as .png file.
plt.figure(dpi=150)
plt.title("triggers_waveform_%s_%ds" % ("BNS", confusion_noise.duration))
plt.xlabel("Time(s)")
plt.ylabel("Strain")
plt.plot(confusion_noise.sample_times, triggers_waveform)
plt.savefig("%s/triggers_waveform_%s_%ds.png" % (output_path, "BNS", confusion_noise.duration))
plt.close()

# Save triggers waveform as .gwf file.
data_triggers_waveform = TimeSeries(triggers_waveform, delta_t=1.0/fs, epoch=epoch)
frame.write_frame("%s/triggers_waveform_%s_%ds.gwf" % (output_path, "BNS", confusion_noise.duration), 
                  "H1:LDAS-STRAIN", data_triggers_waveform)

# Save residual noise as .gwf file.
data_residual_noise = TimeSeries(confusion_noise.data - triggers_waveform, delta_t=1.0/fs, epoch=epoch)
frame.write_frame("%s/residual_noise_%ds.gwf" % (output_path, confusion_noise.duration), 
                  "H1:LDAS-STRAIN", data_residual_noise)

# Save residual noise as .png file.
plt.figure(dpi=150)
plt.title("residual_noise_%ds" % confusion_noise.duration)
plt.xlabel("Time(s)")
plt.ylabel("Strain")
plt.plot(confusion_noise.sample_times, data_residual_noise)
plt.savefig("%s/residual_noise_%ds.png" % (output_path, confusion_noise.duration))
plt.close()
print("Done.")
