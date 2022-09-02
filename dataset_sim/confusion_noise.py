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
from threading import Thread
from time import sleep
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import lal
from pycbc.detector import Detector
from pycbc.waveform import get_td_waveform
from pycbc.frame import write_frame
from pycbc.io import record
from pycbc.types.config import InterpolatingConfigParser
from pycbc.types.timeseries import TimeSeries
from pycbc.distributions import read_params_from_config
from pycbc.distributions.utils import draw_samples_from_config
from tqdm import tqdm


parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('--signal-type', required=True,
                    help="The type of the signal to be injected, "
                    "choose from 'BBH', 'BNS' or 'NSBH'.")
parser.add_argument('--detector-name', required=True,
                    help="The name of the detector to be injected, "
                    "choose from 'ET-D', 'CE'.")
parser.add_argument('--rotation', type=str, required=True, default=False,
                    help="The rotation effect of the earth when project GW onto the detector, "
                    "default is 'False'.")
parser.add_argument('--duration', type=float, required=True,
                    help="The duration of simulated data. In the unit of hour.")
parser.add_argument('--average-time-interval', type=float, required=True,
                    help="The average time interval between two nearby signals.")
parser.add_argument('--seed', type=int, default=0,
                    help="The seed to be used for the random number generator. Default is 0.")
parser.add_argument('--input-config', required=True,
                    help="The config file to be loaded and draw samples from.")
parser.add_argument('--input-params-file', default=False,
                    help="The params file to be loaded and simulate signals from.")
parser.add_argument('--output-folder', required=True,
                    help="The output folder to save simulated noise and parameters file.")
parser.add_argument('--thread-number', type=int, required=True,
                    help="The number of threads to be used.")
parser.add_argument("--force", action="store_true", default=False,
                    help="If the output-file already exists, overwrite it. "
                         "Otherwise, an OSError is raised.")
parser.add_argument("--verbose", action="store_true", default=False,
                    help="Print logging messages.")

opts = parser.parse_args()

if os.path.exists(opts.output_folder) and not opts.force:
    raise OSError("output-file already exists; use --force if you wish to "
                  "overwrite it.")
os.makedirs(opts.output_folder)

gw_type = opts.signal_type
detector = opts.detector_name
rotation = opts.rotation
duration = 3600*opts.duration #  In seconds.
tau = opts.average_time_interval
random_seed = opts.seed
config_path = opts.input_config
params_path = opts.input_params_file
output_path = opts.output_folder
thread_number = opts.thread_number

# Initialise InterpolatingConfigParser class.
cp = InterpolatingConfigParser()
# Read the file.
file = open(config_path, 'r')
cp.read_file(file)
file.close()

# Get the static arguments from the config file.
_, static_params = read_params_from_config(
                        cp, sargs_section='static_params',
                        vargs_section='variable_params')
fs = int(static_params['fs'])
f_low = static_params['f_lower']
f_ref = static_params['f_ref']
approximant = static_params['approximant']


def ht_generator(det='ET-D', signal_type='BBH', model='IMRPhenomTPHM', fs=4096,
                 flow=2, fref=2, tc=0, params=None, rotation=False):

    if signal_type == "BBH":
        # Generate a waveform at the detector-frame.
        hp, hc = get_td_waveform(approximant=model, 
                    mass1=params['mass1'], mass2=params['mass2'],
                    spin1x=params['spin1x'], spin1y=params['spin1y'],
                    spin1z=params['spin1z'], spin2x=params['spin2x'],
                    spin2y=params['spin2y'], spin2z=params['spin2z'],
                    distance=params['distance'], coa_phase=params['coa_phase'],
                    inclination=params['inclination'], f_lower=flow,
                    f_ref=fref, delta_t=1.0/fs)

    elif signal_type in ["NSBH", "BNS"]:
        # Generate a waveform at the detector-frame.
        hp, hc = get_td_waveform(approximant=model, 
                    mass1=params['mass1'], mass2=params['mass2'],
                    spin1z=params['spin1z'],
                    spin2z=params['spin2z'],
                    distance=params['distance'], coa_phase=params['coa_phase'],
                    inclination=params['inclination'], f_lower=flow,
                    f_ref=fref, delta_t=1.0/fs)

    else:
        raise ValueError("`signal_type` must be chosen from ['BBH', 'BNS', 'NSBH'].")

    # Set merger time to 'tc'.
    hp.start_time += tc
    hc.start_time += tc

    # Project GW waveform onto GW detectors.
    ra = params['ra']
    dec = params['dec']
    psi = params['polarization']
    time = hp.start_time

    # TODO: Implement `load_detector_config`.
    if det == 'ET-D':
        det_1 = Detector("E1")
        det_2 = Detector("E2")
        det_3 = Detector("E3")
        num_detectors = 3
    elif det == 'CE':  #  Assume same position and orientation as aLIGO.
        det_1 = Detector("H1")
        det_2 = Detector("L1")
        num_detectors = 2
    else:
        raise NotImplementedError("No such detector.")

    fp_1, fc_1 = det_1.antenna_pattern(
                    right_ascension=ra, declination=dec, polarization=psi, t_gps=tc)
    fp_2, fc_2 = det_2.antenna_pattern(
                    right_ascension=ra, declination=dec, polarization=psi, t_gps=tc)
    if num_detectors == 3:
        fp_3, fc_3 = det_3.antenna_pattern(
                        right_ascension=ra, declination=dec, polarization=psi, t_gps=tc)

    if rotation == "True":
        # Take the rotation of the earth into account by using the "project_wave" function.
        ht_1 = det_1.project_wave(hp=hp, hc=hc, ra=ra, dec=dec, polarization=psi)
        ht_2 = det_2.project_wave(hp=hp, hc=hc, ra=ra, dec=dec, polarization=psi)
        ht_list = [ht_1, ht_2]
        if num_detectors == 3:
            ht_3 = det_3.project_wave(hp=hp, hc=hc, ra=ra, dec=dec, polarization=psi)
            ht_list.append(ht_3)
    elif rotation == "False":
        # Not take the rotation of the earth into account.
        ht_1 = fp_1*hp + fc_1*hc
        ht_2 = fp_2*hp + fc_2*hc
        ht_list = [ht_1, ht_2]
        if num_detectors == 3:
            ht_3 = fp_3*hp + fc_3*hc
            ht_list.append(ht_3)
    else:
        raise NotImplementedError("Must choose from 'True' or 'False'.")

    return ht_list


# Total simulated time.
total_time = int(duration) #s

if os.path.exists(params_path) and os.path.getsize(params_path) > 0:
    params = np.loadtxt(params_path, unpack=True, skiprows=1)
    if params.ndim == 1:
        num_signals = 1
        params = np.expand_dims(params, axis=1)
    else:
        num_signals = len(params[0])
    print("Parameters loaded, numbers of signals: ", num_signals)

else:
    num_signals = int(total_time/tau)
    # Draw samples according to the .ini file.
    samples = draw_samples_from_config(path=config_path, num=num_signals, seed=random_seed)
    params_file = open("%s/injection_params.txt" % output_path, 'w')
    params_file.write(
        """# coa_time mass1 mass2 spin1x spin1y spin1z spin2x spin2y spin2z distance coa_phase inclination ra dec polarization\n""")

    print("Parameters have been drawn, numbers of signals: ", num_signals)
    print("Saving injections' parameters into .txt file...")

    for i in tqdm(range(num_signals)):
        coa_time = np.random.uniform(0, total_time)
        # Create a dict which contains all parameters.
        samples_dict = {}
        for name in samples[:][i].dtype.names:
            samples_dict[name] = samples[:][i][name]
        params_dict = {**samples_dict, **static_params}

        if gw_type in ["BNS", "NSBH"]:
            params_file.write(
                """%.18f\t%.18f\t%.18f\t%.18f\t%.18f\t%.18f\t%.18f\t%.18f\t%.18f\t%.18f\t%.18f\t%.18f\t%.18f\t%.18f\t%.18f\n""" % 
                    (coa_time, params_dict['mass1'], params_dict['mass2'], 
                        0, 0, params_dict['spin1z'], 0, 0, params_dict['spin2z'],
                        params_dict['distance'], params_dict['coa_phase'], params_dict['inclination'],
                        params_dict['ra'], params_dict['dec'], params_dict['polarization']))
        else:
            params_file.write(
                """%.18f\t%.18f\t%.18f\t%.18f\t%.18f\t%.18f\t%.18f\t%.18f\t%.18f\t%.18f\t%.18f\t%.18f\t%.18f\t%.18f\t%.18f\n""" % 
                    (coa_time, params_dict['mass1'], params_dict['mass2'], 
                        params_dict['spin1x'], params_dict['spin1y'], params_dict['spin1z'],
                        params_dict['spin2x'], params_dict['spin2y'], params_dict['spin2z'],
                        params_dict['distance'], params_dict['coa_phase'], params_dict['inclination'],
                        params_dict['ra'], params_dict['dec'], params_dict['polarization']))        
    params_file.close()
    params = np.loadtxt("%s/injection_params.txt" % output_path, unpack=True, skiprows=1)
    if params.ndim == 1:
        params = np.expand_dims(params, axis=1)


class MyThread(Thread):
    def __init__(self, func, args):
        Thread.__init__(self)
        self.func = func
        self.args = args
        self.result = None
    def run(self):
        self.result = self.func(*self.args)
    def get_result(self):
        return self.result

def simulation_on_a_thread(index_start, index_end):
    print("[%d, %d)" % (index_start, index_end))
    # Create empty strain time series.
    time_series = np.linspace(start=0, stop=total_time, num=total_time*fs, endpoint=True)
    if detector == 'CE':
        strain_series_list = [np.zeros(total_time*fs), np.zeros(total_time*fs)]
    elif detector == 'ET-D':
        strain_series_list = [np.zeros(total_time*fs), np.zeros(total_time*fs), np.zeros(total_time*fs)]
    else:
        raise NotImplementedError("No such detector.")

    for i in tqdm(range(index_start, index_end)):
        samples_dict = {'mass1': params[1][i], 'mass2': params[2][i],
                        'spin1x': params[3][i], 'spin1y': params[4][i], 'spin1z': params[5][i],
                        'spin2x': params[6][i], 'spin2y': params[7][i], 'spin2z': params[8][i],
                        'distance': params[9][i], 'coa_phase': params[10][i], 'inclination': params[11][i],
                        'ra': params[12][i], 'dec': params[13][i], 'polarization': params[14][i]}
        coa_time = params[0][i]
        params_dict = {**samples_dict, **static_params}

        strain = ht_generator(det=detector, signal_type=gw_type, model=params_dict['approximant'],
                              tc=coa_time, params=params_dict, rotation=rotation)

        for det_index in range(len(strain_series_list)):
            signal_time = np.array(strain[det_index].sample_times)
            signal_start_time = signal_time[0]
            signal_end_time = signal_time[-1]

            # Inject signal.  # TODO: Try to not use `int()`.
            if signal_start_time >= time_series[0] and signal_end_time <= time_series[-1]:
                time_index = int(signal_start_time*fs+0.5)
                strain_series_list[det_index][time_index:time_index+len(strain[det_index])] += strain[det_index]
            elif signal_start_time <= time_series[0] and signal_end_time <= time_series[-1]:
                time_index = int(signal_end_time*fs+0.5)
                strain_series_list[det_index][:time_index] += strain[det_index][len(strain[det_index])-time_index:]
            elif signal_start_time >= time_series[0] and signal_end_time >= time_series[-1]:
                time_index = int(signal_start_time*fs+0.5)
                strain_series_list[det_index][time_index:] += strain[det_index][:len(strain_series_list[det_index])-time_index]
            else:
                time_index = int((time_series[0]-signal_start_time)*fs+0.5)
                strain_series_list[det_index][:] += strain[det_index][time_index:time_index+len(strain_series_list[det_index])]
    return strain_series_list


print("Simulating signals...")
strain = []
thread_pool = []

index_start = 0
stride = int(np.ceil(num_signals/thread_number))
index_end = index_start + stride

for i in range(thread_number):
    thread_pool.append(MyThread(simulation_on_a_thread, (index_start, index_end)))
    index_start += stride
    if index_end + stride > num_signals - 1:
        index_end = num_signals
    else:
        index_end += stride
for t in thread_pool:
    t.start()
    sleep(5)
for t in thread_pool:
    t.join()
    strain.append(t.get_result())

print("Combining strain data on all threads for each detector...")
strain = np.array(strain)
strain_series_list = np.sum(strain, axis=0)

print("Simulation finished, start to write results...")
time_series = np.linspace(start=0, stop=total_time, num=total_time*fs, endpoint=True)
for det_index in tqdm(range(len(strain_series_list))):
    # Plot the confusion noise made by overlapping GW signals.
    duration = len(time_series)/fs
    plt.figure(dpi=150)
    plt.title("confusion_noise_%s_%d_%s_%ds" % (detector, det_index+1, gw_type, duration))
    plt.xlabel("Time(s)")
    plt.ylabel("Strain")
    plt.plot(time_series, strain_series_list[det_index])
    plt.savefig("%s/confusion_noise_%s_%d_%s_%ds.png" % (output_path, detector, det_index+1, gw_type, duration))
    plt.close()

    # Save confusion noise as .gwf files.
    epoch = lal.LIGOTimeGPS(time_series[0])
    data_confution_noise = TimeSeries(strain_series_list[det_index], delta_t=1.0/fs, epoch=epoch)
    write_frame("%s/confusion_noise_%s_%d_%s_%ds.gwf" % (output_path, detector, det_index+1, gw_type, duration), 
                "H1:LDAS-STRAIN", data_confution_noise)

print("Finished.")
