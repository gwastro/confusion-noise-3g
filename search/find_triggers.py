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
import shutil
import argparse
import h5py
from math import log
import numpy as np
import matplotlib.pyplot as plt
from pycbc import frame
from pycbc.filter import matched_filter
from pycbc.filter.matchedfilter import make_frequency_series, sigmasq
from pycbc.types.timeseries import TimeSeries
from pycbc.types.frequencyseries import load_frequencyseries
from pycbc.psd import interpolate, inverse_spectrum_truncation
from pycbc.waveform import get_fd_waveform
from pycbc.conversions import eta_from_mass1_mass2
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
parser.add_argument('--det-noise-path', required=True,
                    help="The detector noise to be loaded.")
parser.add_argument('--bank-path', required=True,
                    help="The template bank to be used.")
parser.add_argument('--template-id-start', required=True,
                    help="The id of the first template to be used.")
parser.add_argument('--template-id-end', required=True,
                    help="The id of the last template to be used.")
parser.add_argument('--psd-path', required=False,
                    help="The pre-generated PSD.")
parser.add_argument('--snr-threshold', required=True,
                    help="The SNR threshold to be used.")
parser.add_argument('--fmin', required=True,
                    help="The fmin to be used.")
parser.add_argument('--output-folder', required=True,
                    help="The output folder to save triggers' files.")
parser.add_argument('--save-fig', required=True,
                    help="Whether to store triggers' SNR time series pictures.")
parser.add_argument("--force", action="store_true", default=False,
                    help="If the output-file already exists, overwrite it. "
                         "Otherwise, an OSError is raised.")
parser.add_argument("--verbose", action="store_true", default=False,
                    help="Print logging messages.")

opts = parser.parse_args()
init_logging(opts.verbose)

confusion_noise_path = opts.confusion_noise_path
det_noise_path = opts.det_noise_path
bank_path = opts.bank_path
template_id_start = int(opts.template_id_start)
template_id_end = int(opts.template_id_end)
psd_path = opts.psd_path
snr_threshold = np.float64(opts.snr_threshold)
fmin = np.float64(opts.fmin)
output_folder = opts.output_folder
save_fig = opts.save_fig


def chirp_time(m_total, m_eta, f_min):
    const_Msun = 1.989e30                      # kg
    const_G = 6.67259e-11                      # m^3/kg/s^2
    const_c = 2.99792458e8                     # m/s
    const_Euler = 0.577216

    mu = (np.pi*const_G*m_total*const_Msun*f_min/const_c**3)**(1./3)

    # 3.5 PN
    chirp_t = (5/256/m_eta)*(const_G*m_total*const_Msun/const_c**3)*( \
                (mu**(-8)) +(743/252+11/3*m_eta)*mu**(-6)-32*np.pi/5*mu**(-5)+ \
                    (3058673/508032+5429/504*m_eta+617/72*m_eta**2)*mu**(-4)+ \
                        (13*m_eta/3-7729/252)*np.pi*mu**(-3)+ \
                (6848*const_Euler/105-10052469856691/23471078400+128*np.pi**2/3+ \
                    (3147553127/3048192-451*np.pi**2/12)*m_eta-15211/1728*m_eta**2+ \
                        25565/1296*m_eta**3+6848/105*log(4*mu))*mu**(-2)+ \
                            (14809/378*m_eta**2-75703/756*m_eta-15419335/127008)*np.pi*mu**(-1))
    return chirp_t


confusion_noise = frame.read_frame(location=confusion_noise_path, channels='H1:LDAS-STRAIN')
det_noise = frame.read_frame(location=det_noise_path, channels='H1:LDAS-STRAIN')

fs = int(det_noise.sample_rate)
segment_step = 1800 # s
segment_start = 0
segment_end = 3600
f = h5py.File(bank_path, 'r')
found_triggers = 0
psd_duration = 16

if template_id_end > len(np.array(f['mass1'])) - 1:
    template_id_end = len(np.array(f['mass1'])) - 1
output_path = '%s/template_id_%d_%d_triggers' % (output_folder, template_id_start, template_id_end)
os.makedirs(output_path)
triggers_file = open("%s/triggers_info.txt" % output_path, 'w')
triggers_file.write("# time\t segment_start\t segment_end\t mass1\t mass2\t snr_optimal\t snr_tot_real\t snr_tot_imag\t snr_op_real\t snr_op_imag\t snr_noise_real\t snr_noise_imag\n")

data_tot_segments_list = []
data_op_segments_list = []
data_noise_segments_list = []
segment_path = '%s/data_segments' % output_folder


if not os.path.exists(segment_path):
    print("Saving the frequency-domain version of each data segment into the disk and the list.")
    os.makedirs(segment_path)

    for segment_index in range(int(det_noise.duration/segment_step)-1):
        data_tot_segment_td = TimeSeries(confusion_noise.time_slice(segment_start, segment_end).data + \
                                         det_noise.time_slice(segment_start, segment_end).data,
                                         delta_t=1.0/det_noise.sample_rate, epoch=segment_start)
        if segment_index == 0:
            data_512s = data_tot_segment_td.time_slice(segment_start, segment_start+512)
            psd = data_512s.psd(psd_duration)
            psd = interpolate(psd, data_tot_segment_td.delta_f) # use same df as waveform generation
            psd = inverse_spectrum_truncation(psd, int(psd_duration * data_512s.sample_rate),
                                              low_frequency_cutoff=fmin)
            psd.save("%s/psd_%.2fs_%.2fs.txt" % (segment_path, segment_start, segment_start+512))
        data_tot_segment_fd = data_tot_segment_td.to_frequencyseries()
        data_tot_segment_fd.save("%s/data_tot_segment_%.2fs_%.2fs.hdf" % (segment_path, segment_start, segment_end))
        data_op_segment_fd = confusion_noise.time_slice(segment_start, segment_end).to_frequencyseries()
        data_op_segment_fd.save("%s/data_op_segment_%.2fs_%.2fs.hdf" % (segment_path, segment_start, segment_end))
        data_noise_segment_fd = det_noise.time_slice(segment_start, segment_end).to_frequencyseries()
        data_noise_segment_fd.save("%s/data_noise_segment_%.2fs_%.2fs.hdf" % (segment_path, segment_start, segment_end))
        data_tot_segments_list.append(data_tot_segment_fd)
        data_op_segments_list.append(data_op_segment_fd)
        data_noise_segments_list.append(data_noise_segment_fd)
        segment_start += segment_step
        segment_end += segment_step

else:
    print("Loading the frequency-domain version of each data segment into the list.")
    psd = load_frequencyseries("%s/psd_%.2fs_%.2fs.txt" % (segment_path, 0, 512))
    for segment_index in range(int(det_noise.duration/segment_step)-1):
        data_tot_segment_fd = load_frequencyseries("%s/data_tot_segment_%.2fs_%.2fs.hdf" % (segment_path, segment_start, segment_end))
        data_tot_segments_list.append(data_tot_segment_fd)
        data_op_segment_fd = load_frequencyseries("%s/data_op_segment_%.2fs_%.2fs.hdf" % (segment_path, segment_start, segment_end))
        data_op_segments_list.append(data_op_segment_fd)
        data_noise_segment_fd = load_frequencyseries("%s/data_noise_segment_%.2fs_%.2fs.hdf" % (segment_path, segment_start, segment_end))
        data_noise_segments_list.append(data_noise_segment_fd)
        segment_start += segment_step
        segment_end += segment_step


for index in tqdm(range(template_id_start, template_id_end+1)):
    print("Doing the matched filtering with the template %d." % index)
    m1 = np.array(f['mass1'])[index]
    m2 = np.array(f['mass2'])[index]
    segment_start = 0
    segment_end = 3600

    hp, hc = get_fd_waveform(approximant="IMRPhenomD",
                             mass1=m1, mass2=m2, f_lower=fmin, f_ref=fmin,
                             delta_f=1./data_tot_segments_list[0].duration)
    hp.resize(len(data_tot_segments_list[0]))
    template = hp
    htilde = make_frequency_series(template)
    h_norm = sigmasq(htilde, psd, fmin)
    template_duration = chirp_time(m1+m2, eta_from_mass1_mass2(m1,m2), fmin)

    for segment_index in range(int(det_noise.duration/segment_step)-1):

        snr_tot = matched_filter(template, data_tot_segments_list[segment_index], psd=psd, low_frequency_cutoff=fmin, sigmasq=h_norm)
        snr_op = matched_filter(template, data_op_segments_list[segment_index], psd=psd, low_frequency_cutoff=fmin, sigmasq=h_norm)
        snr_noise = matched_filter(template, data_noise_segments_list[segment_index], psd=psd, low_frequency_cutoff=fmin, sigmasq=h_norm)
        # Avoid the SNR time series wrapper-around effect.
        snr_tot = snr_tot.crop(psd_duration+template_duration, psd_duration)
        snr_op = snr_op.crop(psd_duration+template_duration, psd_duration)
        snr_noise = snr_noise.crop(psd_duration+template_duration, psd_duration)
        if segment_index == 0:
            print("We have dropped first %.2fs in this data segment." % template_duration)

        # Find triggers above the SNR threshold.
        triggers_index = np.where(abs(snr_tot)>snr_threshold)[0]

        if len(triggers_index) > 0:
            found_triggers += 1
            triggers_time = snr_tot.sample_times[triggers_index]
            triggers_tot_snr = snr_tot[triggers_index]
            triggers_op_snr = snr_op[triggers_index]
            triggers_noise_snr = snr_noise[triggers_index]
            print("Template %d has found triggers:\n" % index, triggers_time)

            # Save triggers' SNR time series.
            if save_fig == "True":
                plt.figure(figsize=[10, 4])
                plt.title("template id %d" % index)
                plt.axhline(y=snr_threshold, color='red', label=r"$\rho_{th}$=%.2f" % snr_threshold)
                plt.plot(snr_tot.sample_times, abs(snr_tot))
                plt.ylabel('SNR')
                plt.xlabel('Time (s)')
                plt.legend(loc='lower right')
                plt.savefig("%s/snr_timeseries_%.3fs_%.3fs_id_%d.png" % (output_path, segment_start, segment_end, index))
                plt.close()

            # Save triggers' info.
            for i in range(len(triggers_time)):
                triggers_file.write("%.18f\t %d\t %d\t %.18f\t %.18f\t %.8f\t %.8f\t %.8f\t %.8f\t %.8f\t %.8f\t %.8f\n" % 
                                    (triggers_time[i], segment_start, segment_end, m1, m2, np.sqrt(h_norm),
                                     np.real(triggers_tot_snr[i]), np.imag(triggers_tot_snr[i]),
                                     np.real(triggers_op_snr[i]), np.imag(triggers_op_snr[i]),
                                     np.real(triggers_noise_snr[i]), np.imag(triggers_noise_snr[i])))

        segment_start += segment_step
        segment_end += segment_step

triggers_file.close()

# Delete the folder if there is no trigger.
if found_triggers == 0:
    shutil.rmtree(output_path)
