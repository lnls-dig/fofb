"""Classes for FOFB experiment plots: Power Spectral Density and RMS Beam Size."""

import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy import integrate
from scipy import signal


class FOFBPlots:
    """Class to obtain RMS and PSD plots."""

    def __init__(self):

        # Paths to files
        self.fname_FOFB_off = '' # acquisition file, FOFB rate, FOFB OFF
        self.fname_FAcq_off = '' # acquisition file, FAcq rate, FOFB OFF
        self.fname_FOFB_on = '' # acquisition file, FOFB rate, FOFB ON
        self.fname_FAcq_on = '' # acquisition file, FAcq rate, FOFB ON
        self.fname_FOFB_psd1 = '' # first acquisition file for PSD plot
        self.fname_FOFB_psd2 = '' # second acquisition file for PSD plot
        self.beam_sizes_file = '' # .txt format, optional

        # standard
        self.bpm_type_offsets = { 'M1':-1, 'M2':0, 'C1-1':1, 'C1-2':2, 'C2':3, 'C3-1':4, 'C3-2':5, 'C4':6 }
        self.all_bpms = np.arange(0,160)

        # RMS plot params
        self.xtick_label_size = 19
        self.rms_xlabel_font_size = 18
        self.rms_ylabel_font_size = 18
        self.rms_subtitle_fontsize = 20
        self.rms_ylim_vertical = [0.8, 30]
        self.rms_ylim_horizontal = [0.08, 68]
        self.legend_font_size = 20
        self.box_to_anchor = [1.0028, 1.005]
        self.rms_figsize = [12,8]
        # PSD plot params
        self.label_FOFB_psd1 = ''
        self.label_FOFB_psd2 = ''
        self.bpm_psd = 1
        self.psd_ylim_horizontal = [1e-7, 1]
        self.psd_ylim_vertical = [1e-7, 1]
        self.psd_xlabel_font_size = 18
        self.psd_ylabel_font_size = 18
        self.psd_subtitle_fontsize = 15
        self.psd_figsize = [11,8]

    def _read_beam_sizes_from_file(self):
        """Read beam sizes from .txt file, not necessary if it has not changed.
        """
        beam_sizes_x = np.loadtxt(self.beam_sizes_file, usecols=0)
        beam_sizes_y = np.loadtxt(self.beam_sizes_file, usecols=1)
        return beam_sizes_x, beam_sizes_y


    def _get_beam_sizes(self):  
        """Get beam sizes values.
        """           

        beam_sizes_x = np.array([66.13, 81.16, 50.45, 24.37, 61.37, 30.14, 80.95, 44.17, 43.73, 80.93, 50.39, 24.41,
                                60.79, 29.89, 81.75, 44.33, 43.76, 80.89, 50.29, 24.32, 61.26, 30.02, 81.07, 44.53,
                                43.77, 80.93, 50.42, 24.27, 60.74, 29.83, 81.83, 67.51, 67.08, 81.98, 51.00, 24.64,
                                61.14, 30.24, 80.77, 43.37, 44.2, 81.12, 50.44, 24.81, 60.95, 30.03, 81.17, 43.62,
                                43.84, 81.56, 50.67, 24.47, 61.49, 30.27, 81.02, 43.82, 43.55, 80.72, 50.18, 24.52,
                                61.00, 29.95, 81.42, 66.37, 66.75, 81.21, 50.6, 24.47, 60.57, 29.86, 81.37, 44.01,
                                44.09, 81.05, 50.43, 24.43, 60.97, 29.97, 80.82, 44.09, 44.05, 81.24, 50.61, 24.38,
                                60.76, 29.92, 81.46, 44.27, 44.06, 80.88, 50.33, 24.45, 60.93, 29.92, 81.08, 66.99,
                                66.57, 81.15, 50.47, 24.73, 60.94, 30.03, 81.24, 43.67, 43.89, 81.46, 50.62, 24.5,
                                61.47, 30.24, 81.02, 43.91, 43.68, 80.9,  50.31, 24.55, 61.12, 30.01, 81.69, 44.22,
                                43.57, 81.05, 50.37, 24.34, 61.54, 30.16, 81.29, 67.08, 67.33, 81.25, 50.59, 24.68,
                                60.74, 29.95, 80.8,  43.79, 44.47, 81.73, 50.9,  24.57, 60.87, 30.09, 81.09, 43.63,
                                44.24, 81.08, 50.43, 24.71, 60.96, 30.03, 80.85, 43.66, 44.12, 81.69, 50.83, 24.49,
                                61.11, 30.15, 81.08, 66.08])
        

        beam_sizes_y = np.array([ 8.75,  9.3,  13.25,  5.5,   4.32, 10.39,  7.72,  8.41,  6.82,  7.23, 10.94,  6.64,
                                5.96, 13.48,  9.12,  8.15,  9.46,  9.5,  13.6, 5.69,  4.23,  9.97,  7.34,  8.04,
                                6.82,  6.95, 10.43,  6.45,  6.04, 13.89,  9.48,  8.11,  7.17,  8.66, 13.16,  7.16,
                                5.47, 11.69,  7.65,  6.75,  8.88,  8.28, 11.59,  4.93,  4.44, 10.95,  8.01,  8.29,
                                6.95,  7.62, 11.47,  6.61,  5.7,  12.79,  8.6,   7.73,  9.18,  9.02, 12.86,  5.74,
                                4.84, 11.57,  8.43,  8.15,  7.95,  8.22, 11.73,  5.32,  4.72, 11.39,  8.29,  8.61,
                                7.42,  7.93, 11.89,  6.79, 5.85, 13.13,  8.89,  8.08,  9.33,  9.25, 13.23,  5.52,
                                4.27, 10.2,   7.55,  8.24,  6.8,   7.1,  10.69, 6.09,  5.46, 12.49,  8.56,  7.42,
                                6.96,  8.14, 12.22,  6.45,  5.03, 10.98,  7.35,  6.77,  8.12,  7.76, 10.99,  5.,
                                4.49, 10.86,  7.91,  8.15,  6.95,  7.56, 11.33,  5.96,  4.77, 10.56,  7.16,  6.69,
                                7.65,  7.45, 10.64,  4.92,  4.35, 10.43,  7.56,  7.18,  6.42,  6.81,  9.88,  5.4,
                                5.11, 12.05,  8.49,  8.13,  7.95,  8.43, 12.38,  5.84,  4.63, 10.61,  7.61,  7.83,
                                7.23,  7.42, 10.97,  6.16,  5.71, 13.33,  9.4,   9.05,  8.79,  9.33, 13.73,  6.68,
                                4.82, 10.39,  7.,  6.37])
        
        return beam_sizes_x, beam_sizes_y
    
    def _get_fofb_bpms(self):
        sec_nums = list(range(1,21))
        sec_types = ['M1','M2','C2','C3-1']

        idx_bases = 8*(np.array(sec_nums) - 1)
        offsets = np.array([self.bpm_type_offsets[key] for key in sec_types])
        idx_bpms = np.sort(np.array([idx_bases + offset for offset in offsets]).ravel())

        fofb_bpms = idx_bpms
        fofb_bpms = np.append(fofb_bpms, 159)
        fofb_bpms = np.delete(fofb_bpms, 0)
        
        return fofb_bpms, idx_bpms
    
    def _get_data_FOFB_off(self):

        if self.fname_FOFB_off == '':
            raise ValueError('Class attribute "fname_FOFB_off" not defined. It must be the acquisition filename.')
        
        with open(self.fname_FOFB_off, "rb") as f:
            data = pickle.load(f)

        orbx = data['data']['orbx']
        orby = data['data']['orby']
        fs = data['data']['sampling_frequency']
        num_samples = data['data']['nrsamples_post']

        # Upper limit of integration for the slower acquisition rate and upper for the faster one
        intermediate_frequency = fs/num_samples
        rms_FOFB_x_off, rms_FOFB_y_off, frequency = self._get_spectral_analysis(orbx, orby, fs,
                                                                        numperseg=num_samples,
                                                                        freq_low=intermediate_frequency,
                                                                        freq_high=1000,
                                                                        window='boxcar',
                                                                        method='rms')
        
        if self.fname_FAcq_off == '':
            raise ValueError('Class attribute "fname_FAcq_off not defined". It must be the acquisition filename.')
        
        with open(self.fname_FAcq_off, "rb") as f:
            data = pickle.load(f)

        orbx = data['data']['orbx']
        orby = data['data']['orby']
        fs = data['data']['sampling_frequency']

        # orbx, orby, time = metrics.get_time_orbits(fname_Monit1, 1000, subtract_mean=True)
        rms_Monit1_x, rms_Monit1_y, frequency = self._get_spectral_analysis(orbx, orby, fs,
                                                                            numperseg=10*fs,
                                                                            freq_low=0.1,
                                                                            freq_high=intermediate_frequency,
                                                                            window='boxcar',
                                                                            method='rms')

        rms_FOFB_off_x = np.sqrt(rms_FOFB_x_off[-1]**2 + rms_Monit1_x[-1]**2)
        rms_FOFB_off_y = np.sqrt(rms_FOFB_y_off[-1]**2 + rms_Monit1_y[-1]**2)

        return rms_FOFB_off_x, rms_FOFB_off_y

    
    def _get_data_FOFB_on(self):

        if self.fname_FOFB_on == '':
            raise ValueError('Class attribute "fname_FOFB_on" not defined. It must be the acquisition filename.')

        with open(self.fname_FOFB_on, "rb") as f:
            data = pickle.load(f)

        orbx = data['data']['orbx']
        orby = data['data']['orby']
        fs = data['data']['sampling_frequency']
        num_samples = data['data']['nrsamples_post']


        intermediate_frequency = fs/num_samples


        rms_FOFB_x_100hz, rms_FOFB_y_100hz, frequency = self._get_spectral_analysis(orbx, orby, fs,
                                                                        numperseg=num_samples,
                                                                        freq_low=intermediate_frequency,
                                                                        freq_high=100,
                                                                        window='boxcar',
                                                                        method='rms')

        rms_FOFB_x_1k, rms_FOFB_y_1k, frequency = self._get_spectral_analysis(orbx, orby, fs,
                                                                        numperseg=num_samples,
                                                                        freq_low=intermediate_frequency,
                                                                        freq_high=1000,
                                                                        window='boxcar',
                                                                        method='rms')
        
        if self.fname_FAcq_on == '':
            raise ValueError('Class attribute "fname_FAcq_on not defined". It must be the acquisition filename.')
        with open(self.fname_FAcq_on, "rb") as f:
            data = pickle.load(f)

        orbx = data['data']['orbx']
        orby = data['data']['orby']
        fs = data['data']['sampling_frequency']

        rms_Monit1_x, rms_Monit1_y, frequency = self._get_spectral_analysis(orbx, orby, 1000,
                                                                            freq_low=0.1,
                                                                            freq_high=intermediate_frequency,
                                                                            numperseg=10*fs,
                                                                            window='boxcar',
                                                                            method='rms')



        rms_FOFB_on_x_100hz = np.sqrt(rms_FOFB_x_100hz[-1]**2 + rms_Monit1_x[-1]**2)
        rms_FOFB_on_y_100hz = np.sqrt(rms_FOFB_y_100hz[-1]**2 + rms_Monit1_y[-1]**2)

        rms_FOFB_on_x_1k = np.sqrt(rms_FOFB_x_1k[-1]**2 + rms_Monit1_x[-1]**2)
        rms_FOFB_on_y_1k = np.sqrt(rms_FOFB_y_1k[-1]**2 + rms_Monit1_y[-1]**2)

        return rms_FOFB_on_x_100hz, rms_FOFB_on_y_100hz, rms_FOFB_on_x_1k, rms_FOFB_on_y_1k



    def _get_spectral_analysis(self, orbx, orby, fs, freq_low=0, freq_high=np.inf, window='hann', numperseg=None, numoverlap=None,
                        method='psd', columns=None):

        num_samples = orbx.shape[0]

        orby -= np.mean(orby, axis=0)
        orbx -= np.mean(orbx, axis=0)

        # Default behavior uses a window that corresponds to 1/10 of the size of the sampled signal
        if numperseg is None:
            numperseg = num_samples//10

        # Given a matrix of signals, determines the spectrum of only the required columns
        if columns is None:
            freq_array, aux_x = signal.welch(orbx, fs, window=window, nperseg=numperseg, noverlap=numoverlap, axis=0)
            _, aux_y = signal.welch(orby, fs, window=window, nperseg=numperseg, noverlap=numoverlap, axis=0)
        else:
            freq_array, aux_x = signal.welch(orbx[0:,columns], fs, window=window, nperseg=numperseg, noverlap=numoverlap, axis=0)
            _, aux_y = signal.welch(orby[0:,columns], fs, window=window, nperseg=numperseg, noverlap=numoverlap, axis=0)

        index = [i for i, f in enumerate(freq_array) if (f >= freq_low and f <= freq_high)]
        if len(index) >= 2:
            freq_array = freq_array[index]
            aux_x = aux_x[index]
            aux_y = aux_y[index]
        else:
            raise ValueError('There is no data specified in the frequency interval delimited by freq_low and freq_high\
                            arguments')

        if method == 'psd':
            spectrum_x = aux_x
            spectrum_y = aux_y
        elif method == 'rms':
            spectrum_x = np.sqrt(integrate.cumtrapz(aux_x, x=freq_array, axis=0, initial=0))
            spectrum_y = np.sqrt(integrate.cumtrapz(aux_y, x=freq_array, axis=0, initial=0))
        elif method == 'sqrtpsd':
            spectrum_x = np.sqrt(aux_x)
            spectrum_y = np.sqrt(aux_y)

        return spectrum_x, spectrum_y, freq_array
    
    
    def get_rms_plot(self):

        beam_sizes_x, beam_sizes_y = self._get_beam_sizes()
        fofb_bpms, _ = self._get_fofb_bpms()
        all_bpms = np.arange(0,160)

        rms_FOFB_on_x_100hz, rms_FOFB_on_y_100hz, rms_FOFB_on_x_1k, rms_FOFB_on_y_1k = self._get_data_FOFB_on()
        rms_FOFB_off_x, rms_FOFB_off_y = self._get_data_FOFB_off()

    
        marker_size = 2.5
        fig, axs = plt.subplots(2, 1, constrained_layout=True, figsize=(self.rms_figsize[0],self.rms_figsize[1]))
        bpm_index_array = np.arange(rms_FOFB_on_x_1k.size)

    
        color_plot = 'indigo'
        axs[0].semilogy(bpm_index_array, 100*(rms_FOFB_off_x[all_bpms]/beam_sizes_x[all_bpms]),
                    label='FOFB Off', marker='D', 
                    fillstyle='none', markersize=marker_size, color=color_plot)
        # --- BPMS on the correction loop ---
        axs[0].semilogy(fofb_bpms, 100*(rms_FOFB_off_x[fofb_bpms]/beam_sizes_x[fofb_bpms]), 
                    'D', label='FOFB BPMs',
                    markerfacecolor='k', markeredgecolor='k', markersize=marker_size)

        axs[1].semilogy(bpm_index_array, 100*(rms_FOFB_off_y[all_bpms]/beam_sizes_y[all_bpms]), marker='D', 
                    fillstyle='none', markersize=marker_size, color=color_plot)
        # --- BPMs on the correction loop ---
        axs[1].semilogy(fofb_bpms, 100*(rms_FOFB_off_y[fofb_bpms]/beam_sizes_y[fofb_bpms]), 'D', 
                    markerfacecolor='k', markeredgecolor='k', markersize=marker_size)



        color_plot = 'dodgerblue'
        axs[0].semilogy(bpm_index_array, 100*(rms_FOFB_on_x_1k[all_bpms]/beam_sizes_x[all_bpms]),
                    label='FOFB On (0.1 Hz - 1 kHz)', marker='D', 
                    fillstyle='none', markersize=marker_size, color=color_plot)
        # --- BPMs on the correction loop ---
        axs[0].semilogy(fofb_bpms, 100*(rms_FOFB_on_x_1k[fofb_bpms]/beam_sizes_x[fofb_bpms]), 'D', 
                    markerfacecolor='k', markeredgecolor='k', markersize=marker_size)

        axs[1].semilogy(bpm_index_array, 100*(rms_FOFB_on_y_1k[all_bpms]/beam_sizes_y[all_bpms]), marker='D', 
                    fillstyle='none', markersize=marker_size, color=color_plot)
        # --- BPMs on the correction loop ---
        axs[1].semilogy(fofb_bpms, 100*(rms_FOFB_on_y_1k[fofb_bpms]/beam_sizes_y[fofb_bpms]),'D', 
                    markerfacecolor='k', markeredgecolor='k', markersize=marker_size)



        color_plot = 'darkorange'
        axs[0].semilogy(bpm_index_array, 100*(rms_FOFB_on_x_100hz[all_bpms]/beam_sizes_x[all_bpms]), 
                    label='FOFB On (0.1Hz - 100Hz)', marker='D', 
                    fillstyle='none', markersize=marker_size, color=color_plot)
        # --- BPMS on the correction loop ---
        axs[0].semilogy(fofb_bpms, 100*(rms_FOFB_on_x_100hz[fofb_bpms]/beam_sizes_x[fofb_bpms]), 'D', 
                    markerfacecolor='k', markeredgecolor='k', markersize=marker_size)

        axs[1].semilogy(bpm_index_array, 100*(rms_FOFB_on_y_100hz[all_bpms]/beam_sizes_y[all_bpms]), marker='D', 
                    fillstyle='none', markersize=marker_size, color=color_plot)
        # --- BPMS on the correction loop ---
        axs[1].semilogy(fofb_bpms, 100*(rms_FOFB_on_y_100hz[fofb_bpms]/beam_sizes_y[fofb_bpms]), 'D', 
                    markerfacecolor='k', markeredgecolor='k', markersize=marker_size)


        axs[0].set_xticks(np.arange(0, 160, 8))
        axs[0].tick_params(which='major', width=1.5, length=5, direction='in', axis='both',
                        labelsize=self.xtick_label_size)
        axs[0].tick_params(which='minor', width=1, length=2.5, direction='in', axis='y')
        axs[0].legend(bbox_to_anchor=(self.box_to_anchor[0], self.box_to_anchor[1]), fontsize=self.legend_font_size, ncol=4, loc='upper right')


        axs[0].margins(x=0.005, y = 1.)
        axs[0].set_ylabel('RMS / beam size (%)', fontsize=self.rms_ylabel_font_size)
        axs[0].set_title('Horizontal', fontsize=self.rms_subtitle_fontsize, pad=20)
        axs[0].set_ylim(self.rms_ylim_horizontal[0], self.rms_ylim_horizontal[1])
        axs[0].grid(which='both', color=(0.9,0.9,0.9))


        axs[1].set_xticks(np.arange(0, 160, 8))
        axs[1].tick_params(which='major', width=1.5, length=5, direction='in', axis='both', 
                        labelsize=self.xtick_label_size)
        axs[1].tick_params(which='minor', width=1, length=2.5, direction='in', axis='y')
        axs[1].margins(x=0.005, y = 5.)
        axs[1].set_ylabel('RMS / beam size (%)', fontsize=self.rms_ylabel_font_size)
        axs[1].set_xlabel('BPM Index', fontsize=self.rms_xlabel_font_size)
        axs[1].set_title('Vertical', fontsize=self.rms_subtitle_fontsize, pad = 20)
        axs[1].grid(which='both', color=(0.9,0.9,0.9))
        axs[1].set_ylim(self.rms_ylim_vertical[0], self.rms_ylim_vertical[1])

        for ax in axs.flat:
            ax.label_outer()

        plt.show()
        plt.savefig('RMS_plot.png')


    def get_psd_plot(self):

        _, idx_bpms = self._get_fofb_bpms()

        with open(self.fname_FOFB_psd1, "rb") as f:
            data = pickle.load(f)

        orbx = data['data']['orbx']
        orby = data['data']['orby']
        fs = data['data']['sampling_frequency']

        psd_FOFB_on_x1, psd_FOFB_on_y1, frequency = self._get_spectral_analysis(orbx, orby, fs,
                                                                            window='boxcar',
                                                                            method='psd')
        
        with open(self.fname_FOFB_psd2, "rb") as f:
            data = pickle.load(f)

        orbx = data['data']['orbx']
        orby = data['data']['orby']
        fs = data['data']['sampling_frequency']

        psd_FOFB_on_x2, psd_FOFB_on_y2, frequency = self._get_spectral_analysis(orbx, orby, fs,
                                                                            window='boxcar',
                                                                            method='psd')
        
        fig, axs = plt.subplots(2, 1, constrained_layout=True, figsize=(self.psd_figsize[0], self.psd_figsize[1]))

        fig.suptitle(f'BPM Index = {self.bpm_psd}')

        axs[0].loglog(frequency, psd_FOFB_on_x1[:, idx_bpms[self.bpm_psd]], label=self.label_FOFB_psd1)
        axs[0].loglog(frequency, psd_FOFB_on_x2[:, idx_bpms[self.bpm_psd]], label=self.label_FOFB_psd2)
        axs[0].set_xlabel('Frequency [Hz]',  fontsize=self.psd_xlabel_font_size)
        axs[0].set_ylabel('PSD [$\\mu m/  \\sqrt{Hz}$]',  fontsize=self.psd_ylabel_font_size)
        axs[0].set_title('Horizontal', fontsize=self.psd_subtitle_fontsize, pad = 10)
        axs[0].legend()
        axs[0].grid()

        axs[1].loglog(frequency, psd_FOFB_on_y1[:, idx_bpms[self.bpm_psd]])
        axs[1].loglog(frequency, psd_FOFB_on_y2[:, idx_bpms[self.bpm_psd]])
        axs[1].set_xlabel('Frequency [Hz]',  fontsize=self.psd_xlabel_font_size)
        axs[1].set_ylabel('PSD [$\\mu m/  \\sqrt{Hz}$]',  fontsize=self.psd_ylabel_font_size)
        axs[1].set_title('Vertical', fontsize=self.psd_subtitle_fontsize, pad = 10)
        axs[1].grid()

        axs[0].set_ylim([self.psd_ylim_horizontal[0], self.psd_ylim_horizontal[1]])
        axs[1].set_ylim([self.psd_ylim_vertical[0], self.psd_ylim_vertical[1]])

        plt.show()
        plt.savefig('PSD_plot.png')





