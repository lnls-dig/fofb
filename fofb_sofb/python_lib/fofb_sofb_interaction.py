"""."""

import numpy as np
import matplotlib.pyplot as mplt


class FOFB_SOFB:
    """."""

    def __init__(self, mat_sofb, mat_fofb):
        """."""
        self.mat_sofb = mat_sofb
        self.mat_fofb = mat_fofb
        self.imat_sofb = np.zeros(mat_sofb.shape[::-1], dtype=float)
        self.imat_fofb = np.zeros(mat_fofb.shape[::-1], dtype=float)
        self._bpmxenbl_sofb = np.ones(160, dtype=bool)
        self._bpmyenbl_sofb = np.ones(160, dtype=bool)
        self._bpmxenbl_fofb = np.ones(160, dtype=bool)
        self._bpmyenbl_fofb = np.ones(160, dtype=bool)
        self._chenbl_sofb = np.ones(120, dtype=bool)
        self._cvenbl_sofb = np.ones(160, dtype=bool)
        self._rfenbl_sofb = True
        self._chenbl_fofb = np.ones(80, dtype=bool)
        self._cvenbl_fofb = np.ones(80, dtype=bool)
        self._rfenbl_fofb = True
        self.minsv_sofb = 0.0
        self.regc_sofb = 0.0
        self.minsv_fofb = 0.0
        self.regc_fofb = 0.0

        self.gaini_sofb = 0.5
        self.gaini_fofb = 2500
        self.gainp_sofb = 0
        self.gainp_fofb = 0
        self.gaini_down = 0.02*25
        self.gaini_orbu = 0.02*25
        self.f0_fil_sofb = 25
        self.f0_fil_fofb = 25e3
        self.f0_plt_sofb = 100
        self.f0_plt_fofb = 10e3
        self.f0_fil_down = 25
        self.f0_ac_fofb = 1
        self.dly_fil_sofb = 30e-3
        self.dly_plt_sofb = 10e-3
        self.dly_fil_fofb = 200e-6
        self.dly_plt_fofb = 50e-6
        self.proj_null_space = False
        self.zero_distortion = False
        self.ac_fofb = False
        self.imat_sofb = self.calc_inverse_matrix('sofb')
        self.imat_fofb = self.calc_inverse_matrix('fofb')

    def __str__(self):
        """."""
        tmpf = '{:12s}: {:.2f} {:s}\n'
        tmps = '{:12s}: {:s}\n'
        stg = ''

        stg1 = ''
        stg1 += tmpf.format('minsv_sofb', self.minsv_sofb, '')
        stg1 += tmpf.format('regc_sofb', self.regc_sofb, '')
        stg1 += tmpf.format('gaini_sofb', self.gaini_sofb, '')
        stg1 += tmpf.format('gainp_sofb', self.gainp_sofb, '')
        stg1 += tmpf.format('f0_fil_sofb', self.f0_fil_sofb, '[Hz]')
        stg1 += tmpf.format('f0_plt_sofb', self.f0_plt_sofb, '[Hz]')
        stg1 += tmpf.format('dly_fil_sofb', self.dly_fil_sofb*1e3, '[ms]')
        stg1 += tmpf.format('dly_plt_sofb', self.dly_plt_sofb*1e3, '[ms]')
        stg1 = stg1.splitlines()
        size1 = max([len(st) for st in stg1]) + 8

        stg2 = ''
        stg2 += tmpf.format('minsv_fofb', self.minsv_fofb, '')
        stg2 += tmpf.format('regc_fofb', self.regc_fofb, '')
        stg2 += tmpf.format('gaini_fofb', self.gaini_fofb, '')
        stg2 += tmpf.format('gainp_fofb', self.gainp_fofb, '')
        stg2 += tmpf.format('f0_fil_fofb', self.f0_fil_fofb/1e3, '[kHz]')
        stg2 += tmpf.format('f0_plt_fofb', self.f0_plt_fofb/1e3, '[kHz]')
        stg2 += tmpf.format('dly_fil_fofb', self.dly_fil_fofb*1e6, '[us]')
        stg2 += tmpf.format('dly_plt_fofb', self.dly_plt_fofb*1e6, '[us]')
        stg2 = stg2.splitlines()
        size2 = max([len(st) for st in stg2]) + 8

        stg3 = ''
        stg3 += tmpf.format('gaini_down', self.gaini_down, '')
        stg3 += tmpf.format('gaini_orbu', self.gaini_orbu, '')
        stg3 += tmpf.format('f0_fil_down', self.f0_fil_down, '[Hz]')
        stg3 += tmpf.format('f0_ac_fofb', self.f0_ac_fofb, '[Hz]')
        stg3 += tmps.format('ac_fofb', str(self.ac_fofb))
        stg3 += tmps.format('proj_null_space', str(self.proj_null_space))
        stg3 += tmps.format('zero_distortion', str(self.zero_distortion))
        stg3 = stg3.splitlines()

        stg += '{:s}{:s}{:s}\n'.format(
            'SOFB Parameters:'.ljust(size1), 'FOFB Parameters:'.ljust(size2),
            'Interaction Parameters:')
        for i, (st1, st2, st3) in enumerate(zip(stg1, stg2, stg3)):
            stg += st1.ljust(size1) + st2.ljust(size2) + st3 + '\n'
        i += 1
        for st1, st2 in zip(stg1[i:], stg2[i:]):
            stg += st1.ljust(size1) + st2.ljust(size2) + '\n'

        stg1 = self._print_enable_lists(
            self._bpmxenbl_sofb, self._bpmyenbl_sofb,
            self._chenbl_sofb, self._cvenbl_sofb, self._rfenbl_sofb)
        stg2 = self._print_enable_lists(
            self._bpmxenbl_fofb, self._bpmyenbl_fofb,
            self._chenbl_fofb, self._cvenbl_fofb, self._rfenbl_fofb)

        stg1 = stg1.splitlines()
        stg2 = stg2.splitlines()
        size = max([len(st) for st in stg1]) + 8

        stg += '\n{:s}{:s}\n'.format(
            'SOFB Enable Lists:'.ljust(size), 'FOFB Enable Lists:')
        for st1, st2 in zip(stg1, stg2):
            stg += st1.ljust(size) + st2 + '\n'
        return stg

    @staticmethod
    def _print_enable_lists(bpmx, bpmy, ch, cv, rf):
        bpmx = np.where(np.roll(bpmx, 1).reshape(20, -1), 'o', 'x')
        bpmy = np.where(np.roll(bpmy, 1).reshape(20, -1), 'o', 'x')
        ch = np.where(np.roll(ch, 1).reshape(20, -1), 'o', 'x')
        cv = np.where(np.roll(cv, 1).reshape(20, -1), 'o', 'x')

        bx = 'BPMx'
        by = 'BPMy'
        cx = 'CH'
        cy = 'CV'
        stg = f'RF = {str(rf):s}\n'
        stg += '      '
        stg += f'{bx.center(bpmx.shape[-1])}  '
        stg += f'{by.center(bpmy.shape[-1])}  '
        stg += f'{cx.center(ch.shape[-1])}  '
        stg += f'{cy.center(cv.shape[-1])}  '
        stg += '\n'
        for i in range(20):
            stg += f'{i+1:02d} -> '
            stg += ''.join(bpmx[i]) + '  '
            stg += ''.join(bpmy[i]) + '  '
            stg += ''.join(ch[i]) + '  '
            stg += ''.join(cv[i]) + '  '
            stg += '\n'
        return stg

    @property
    def bpmxenbl_sofb(self):
        """."""
        return self._bpmxenbl_sofb

    @bpmxenbl_sofb.setter
    def bpmxenbl_sofb(self, value):
        self._bpmxenbl_sofb = np.array(value, dtype=bool)
        self.imat_sofb = self.calc_inverse_matrix('sofb')

    @property
    def bpmyenbl_sofb(self):
        """."""
        return self._bpmyenbl_sofb

    @bpmyenbl_sofb.setter
    def bpmyenbl_sofb(self, value):
        self._bpmyenbl_sofb = np.array(value, dtype=bool)
        self.imat_sofb = self.calc_inverse_matrix('sofb')

    @property
    def bpmxenbl_fofb(self):
        """."""
        return self._bpmxenbl_fofb

    @bpmxenbl_fofb.setter
    def bpmxenbl_fofb(self, value):
        self._bpmxenbl_fofb = np.array(value, dtype=bool)
        self.imat_fofb = self.calc_inverse_matrix('fofb')

    @property
    def bpmyenbl_fofb(self):
        """."""
        return self._bpmyenbl_fofb

    @bpmyenbl_fofb.setter
    def bpmyenbl_fofb(self, value):
        self._bpmyenbl_fofb = np.array(value, dtype=bool)
        self.imat_fofb = self.calc_inverse_matrix('fofb')

    @property
    def chenbl_sofb(self):
        """."""
        return self._chenbl_sofb

    @chenbl_sofb.setter
    def chenbl_sofb(self, value):
        self._chenbl_sofb = np.array(value, dtype=bool)
        self.imat_sofb = self.calc_inverse_matrix('sofb')

    @property
    def cvenbl_sofb(self):
        """."""
        return self._cvenbl_sofb

    @cvenbl_sofb.setter
    def cvenbl_sofb(self, value):
        self._cvenbl_sofb = np.array(value, dtype=bool)
        self.imat_sofb = self.calc_inverse_matrix('sofb')

    @property
    def rfenbl_sofb(self):
        """."""
        return self._rfenbl_sofb

    @rfenbl_sofb.setter
    def rfenbl_sofb(self, value):
        self._rfenbl_sofb = np.array(value, dtype=bool)
        self.imat_sofb = self.calc_inverse_matrix('sofb')

    @property
    def chenbl_fofb(self):
        """."""
        return self._chenbl_fofb

    @chenbl_fofb.setter
    def chenbl_fofb(self, value):
        self._chenbl_fofb = np.array(value, dtype=bool)
        self.imat_fofb = self.calc_inverse_matrix('fofb')

    @property
    def cvenbl_fofb(self):
        """."""
        return self._cvenbl_fofb

    @cvenbl_fofb.setter
    def cvenbl_fofb(self, value):
        self._cvenbl_fofb = np.array(value, dtype=bool)
        self.imat_fofb = self.calc_inverse_matrix('fofb')

    @property
    def rfenbl_fofb(self):
        """."""
        return self._rfenbl_fofb

    @rfenbl_fofb.setter
    def rfenbl_fofb(self, value):
        self._rfenbl_fofb = value
        self.imat_fofb = self.calc_inverse_matrix('fofb')

    def calc_inverse_matrix(self, system='sofb', return_svd=False):
        """."""
        if system.lower().startswith('fofb'):
            bpm_enbl = np.r_[self._bpmxenbl_fofb, self._bpmyenbl_fofb]
            cor_enbl = np.r_[
                self._chenbl_fofb, self._cvenbl_fofb, self._rfenbl_fofb]
            mat = self.mat_fofb
            min_sv = self.minsv_fofb
            regc = self.regc_fofb
        else:
            bpm_enbl = np.r_[self._bpmxenbl_sofb, self._bpmyenbl_sofb]
            cor_enbl = np.r_[
                self._chenbl_sofb, self._cvenbl_sofb, self._rfenbl_sofb]
            mat = self.mat_sofb
            min_sv = self.minsv_sofb
            regc = self.regc_sofb

        sel_mat = bpm_enbl[:, None] * cor_enbl[None, :]
        mats = mat[sel_mat].reshape(np.sum(bpm_enbl), np.sum(cor_enbl))
        uuu, sing, vvv = np.linalg.svd(mats, full_matrices=False)
        idcs = sing > min_sv
        singr = sing[idcs]

        regc = regc * regc
        inv_s = np.zeros(sing.size, dtype=float)
        inv_s[idcs] = singr/(singr*singr + regc)

        singp = np.zeros(sing.size, dtype=float)
        singp[idcs] = 1/inv_s[idcs]
        imat = np.dot(vvv.T*inv_s, uuu.T)

        inv_mat = np.zeros(mat.shape, dtype=float).T
        inv_mat[sel_mat.T] = imat.ravel()

        if system.lower().startswith('fofb'):
            # FOFB will never correct with RF:
            inv_mat[-1] *= 0

        if return_svd:
            return inv_mat, (uuu, sing, vvv)
        return inv_mat

    @staticmethod
    def get_filter(s, f0, delay=0, mode='lowpass'):
        """."""
        w0 = 2*np.pi*f0
        fil = 1/(1+s/w0)
        if mode.lower().startswith('high'):
            fil = 1 - fil
        fil *= np.exp(-delay*s)
        return fil

    @staticmethod
    def get_controller(s, mat, gaini=1, gainp=0):
        """."""
        return (gaini/s + gainp) * mat

    def get_system_matrices(self, freq):
        """."""
        s = 2j*np.pi*freq

        Hs = self.get_filter(
            s, f0=self.f0_fil_sofb, delay=self.dly_fil_sofb)
        Hf = self.get_filter(
            s, f0=self.f0_fil_fofb, delay=self.dly_fil_fofb)

        Cs = self.get_controller(
            s, mat=self.imat_sofb, gaini=self.gaini_sofb,
            gainp=self.gainp_sofb)
        if self.proj_null_space:
            Cs -= Cs @ self.mat_fofb @ self.imat_fofb
        if self.zero_distortion:
            nsel = ~np.r_[self.bpmxenbl_fofb, self.bpmyenbl_fofb]
            Cs[:, nsel] = 0

        Cf = self.get_controller(
            s, mat=self.imat_fofb, gaini=self.gaini_fofb,
            gainp=self.gainp_fofb)
        if self.ac_fofb:
            Cf *= self.get_filter(s, f0=self.f0_ac_fofb, mode='highpass')

        Hk = self.get_filter(s, f0=self.f0_fil_down)
        Dk = self.get_controller(
            s, mat=self.imat_sofb @ self.mat_fofb, gaini=self.gaini_down)
        Dk *= Hk

        Ou = self.get_controller(
            s, mat=self.mat_sofb, gaini=self.gaini_orbu)

        Gs = self.mat_sofb * self.get_filter(
            s, f0=self.f0_plt_sofb, delay=self.dly_plt_sofb)
        Gf = self.mat_fofb * self.get_filter(
            s, f0=self.f0_plt_fofb, delay=self.dly_plt_fofb)

        CsHs = Cs * Hs
        CfHf = Cf * Hf
        CfOuCsHs = Cf @ Ou @ CsHs
        CfHf_CfOuCsHs = CfHf + CfOuCsHs
        T = Gs @ CsHs + (Gs @ Dk + Gf) @ CfHf_CfOuCsHs

        Td2y = np.linalg.pinv(np.eye(320) + T)
        Td2kf = - CfHf_CfOuCsHs @ Td2y
        Td2ks = - CsHs @ Td2y - Dk @ Td2kf

        return Td2y, Td2kf, Td2ks

    def calc_singular_values(self, freqs, only_fofb=False, print_status=True):
        """."""
        svals_d2y, svals_d2kf, svals_d2ks = [], [], []
        for i, freq in enumerate(freqs):
            if print_status and not i % 10:
                print(f'{i:03d}/{freqs.size:03d} -> {freq:.3f} Hz')
            try:
                Td2y, Td2kf, Td2ks = self.get_system_matrices(freq)
                if only_fofb:
                    bpmenbl = np.r_[self.bpmxenbl_fofb, self.bpmyenbl_fofb]
                    Td2y = Td2y[bpmenbl, :]
                sd2y = np.linalg.svd(
                    Td2y, full_matrices=False, compute_uv=False)
                sd2kf = np.linalg.svd(
                    Td2kf, full_matrices=False, compute_uv=False)
                sd2ks = np.linalg.svd(
                    Td2ks, full_matrices=False, compute_uv=False)
                svals_d2y.append(sd2y)
                svals_d2kf.append(sd2kf)
                svals_d2ks.append(sd2ks)
            except np.linalg.LinAlgError:
                svals_d2y.append(np.full(320, np.nan))
                svals_d2kf.append(np.full(161, np.nan))
                svals_d2ks.append(np.full(281, np.nan))
        svals_d2y = np.array(svals_d2y)
        svals_d2kf = np.array(svals_d2kf)
        svals_d2ks = np.array(svals_d2ks)
        if print_status:
            print('finished!')
        return svals_d2y, svals_d2kf, svals_d2ks

    def sweep_parameter(
            self, freqs, param, values, names, only_fofb=False,
            print_status=True):
        """."""
        svals = []
        init = getattr(self, param)
        for i, val in enumerate(values):
            setattr(self, param, val)
            if print_status:
                print('\nCalulating configuration: '+names[i])
                print('\n'.join(str(self).splitlines()[:9]))
            svals.append(self.calc_singular_values(
                freqs, only_fofb=only_fofb, print_status=print_status))

        setattr(self, param, init)
        return svals

    def compare_loop_states(
            self, freqs, variables, states, names, only_fofb=False,
            print_status=True):
        """."""
        svals = []
        init = [getattr(self, var) for var in variables]

        for i, stt in enumerate(states):
            for j, var in enumerate(variables):
                val = type(init[j])(0)
                if stt & (1 << j):
                    val = init[j]
                setattr(self, var, val)
            if print_status:
                print('\nCalculating configuration: '+names[i])
                print('\n'.join(str(self).splitlines()[:9]))
            svals.append(self.calc_singular_values(
                freqs, only_fofb=only_fofb, print_status=print_status))

        for i, var in enumerate(variables):
            setattr(self, var, init[i])
        return svals

    def plot_singular_values(
            self, freqs, svals, names, svs=0, hide_propties=None):
        """."""
        fig = mplt.figure(figsize=(11, 10))
        gs = mplt.GridSpec(
            3, 2, hspace=0.2, left=0.115, right=0.98, top=0.89, bottom=0.08,
            wspace=0.04, width_ratios=[4, 2])
        ady = fig.add_subplot(gs[0, 0])
        adf = fig.add_subplot(gs[1, 0], sharex=ady)
        ads = fig.add_subplot(gs[2, 0], sharex=ady)
        atx = fig.add_subplot(gs[:, 1])
        atx.axis('off')

        txt = self.__str__()
        stg = '\n'.join([t[:25].rstrip() for t in txt.splitlines()[:9]])
        stg += '\n\n'
        stg += '\n'.join([t[33:58].rstrip() for t in txt.splitlines()[:9]])
        stg += '\n\n'
        stg += '\n'.join([t[66:].rstrip() for t in txt.splitlines()[:9]])
        stg += '\n\n'
        stg += '\n'.join([t[52:].rstrip() for t in txt.splitlines()[10:]])

        if hide_propties is not None:
            stg = '\n'.join([
                t for t in stg.splitlines()
                if not any(x in t for x in hide_propties)])

        atx.text(
            0, 1, stg, fontsize='xx-small', fontfamily='monospace',
            horizontalalignment='left', verticalalignment='top',
            transform=atx.transAxes)

        cycle = mplt.rcParams['axes.prop_cycle'].by_key()['color']

        fig.suptitle('Largest Singular Value')

        idx = 0  # transfer from d to y
        for i, name in enumerate(names):
            lines = ady.plot(
                freqs, 10*np.log10(svals[i][idx][:, svs]), color=cycle[i])
            lines[0].set_label(name)
        ady.set_xscale('log')
        ady.legend(
            loc='lower center', fontsize='xx-small', ncol=4,
            bbox_to_anchor=(0.5, 1))
        ady.set_ylabel('D --> Orbit\n[dB]', fontsize='small')

        idx = 1  # transfer from d to Kf
        for i, name in enumerate(names):
            lines = adf.plot(freqs, svals[i][idx][:, svs], color=cycle[i])
            lines[0].set_label(name)
        adf.set_xscale('log')
        adf.set_yscale('log')
        adf.set_ylabel('D --> Kicks FOFB\n[urad/um]', fontsize='small')

        idx = 2  # transfer from d to Ks
        for i, name in enumerate(names):
            lines = ads.plot(freqs, svals[i][idx][:, svs], color=cycle[i])
            lines[0].set_label(name)
        ads.set_xscale('log')
        ads.set_yscale('log')
        ads.set_ylabel('D --> Kicks SOFB\n[urad/um]', fontsize='small')
        ads.set_xlabel('Frequency [Hz]')

        fig.show()
        return fig, (ady, adf, ads)
