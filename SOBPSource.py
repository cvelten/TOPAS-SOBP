from __future__ import division

from argparse import ArgumentParser
import os

import numpy as np
import pandas as pd


class SOBPSource(object):
    RecommendedValues = pd.DataFrame(
        {
            50.: [1.72, 1.71, 1.71, 1.71, 1.71, 1.73],
            100.: [1.63, 1.68, 1.68, 1.67, 1.68, 1.68],
            150.: [1.69, 1.68, 1.66, 1.66, 1.64, 1.62],
            200.: [1.62, 1.60, 1.58, 1.56, 1.54, 1.53],
            250.: [1.57, 1.56, 1.53, 1.52, 1.51, 1.48],
        },
        index=[0.15, 0.20, 0.25, 0.30, 0.35, 0.40])

    @staticmethod
    def get_p(energy, chi):
        """Interpolate from recommended power law values.

        Args:
            energy (float): Maximum proton energy in MeV.
            chi (float): SOBP width (fraction of the maximum range).

        Returns:
            float: Interpolated power law parameter.
        """
        # add a new energy/width column/row to the DataFrame
        tmp = SOBPSource.RecommendedValues.copy()
        tmp.loc[chi, energy] = np.nan
        tmp = tmp.sort_index(0).sort_index(1)
        # interpolate width values for a given energy first
        tmp = tmp.interpolate(method='index', axis='index')
        # interpolate between different energies
        tmp = tmp.interpolate(method='index', axis='columns')
        return tmp.loc[chi, energy]

    @staticmethod
    def _rk(k, n, chi, R0):
        """Range of the k-th beamlet.

        Args:
            k (int): Number of the beamlet, within [0, n].
            n (int): Maximum of beamlets.
            chi (float): SOBP width (fraction of the maximum range).
            R0 (float): Range of the highest energy proton.

        Returns:
            float: Proton range of the k-th beamlet
        """
        return (1 - (1 - k / n) * chi) * R0

    @staticmethod
    def _ek(r_k, alpha, p0):
        """Energy of the k-th beamlet.

        Args:
            r_k (float): Proton range of the k-th beamlet.
            alpha (float): Power law scaling parameter .
            p0 (float): Conventional power law parameter for energy-to-range conversion.

        Returns:
            float: Energy in MeV of the k-th beamlet."""
        return np.power(r_k / alpha, 1 / p0)

    @staticmethod
    def _wk(k, n, p):
        """Weight of the k-th beamlet.

        Args:
            k (int): Number of the beamlet, within [0, n].
            n (int): Maximum of beamlets.
            p (float): Power law parameter p.

        Returns:
            float: Weight of the k-th beamlet.
        """
        a = 1 / (2 * n)
        b = 1 - 1 / p
        if k == 0:
            # k = 0: lowest energy beamlet
            return 1 - np.power(1 - a, b)
        elif k < n:
            return np.power(1 -
                            (k - 0.5) / n, b) - np.power(1 - (k + 0.5) / n, b)
        else:
            # k = n: highest energy beamlet
            return np.power(a, b)

    def __init__(self, alpha=2.2E-3, p0=1.77):
        self.alpha = alpha
        self.p0 = p0
        self.energy = 0
        self.chi = 0
        self.E, self.T = [], []
        pass

    def generate_weights(self, energy, chi, p, n=19, delta=0.05):
        """Generate spread-out Bragg peak beamlet energies and weights.

        Args:
            energy (float): Maximum proton energy.
            chi (float): SOBP width (fraction of the maximum range).
            p (float): Power law parameter p.
            n (int): Minimum number of beamlets.
            delta (float): Beamlet spacing, used if it leads to more than n beamlets.

        Returns:
            tuple(list, list): List of beamlet energies and cumulative weights (TOPAS time points).
        """
        R0 = self.alpha * np.power(energy, self.p0)

        # set number of beamlets according to the minimum spacing (in centimeters)
        # or use minimum number of beamlets
        n_ = int(np.ceil(chi * R0 / delta))
        n = n if n_ < n else n_

        E, T = [], [0]
        for k in range(n + 1):
            r = SOBPSource._rk(k, n, chi, R0)
            e = SOBPSource._ek(r, self.alpha, self.p0)
            w = SOBPSource._wk(k, n, p)
            E.append(e)
            T.append(T[-1] + w)

        self.energy, self.chi = energy, chi
        self.E, self.T = E, T[1:]
        return E, T

    def print(self):
        n = len(self.E)
        print('# SOBP Time Feature (0 to 1000 ms)')
        print('# Energy = {:3.1f} MeV'.format(self.E[-1]))
        print('# Width  = {:1.2f}'.format(self.chi))

        print('s:Tf/BeamEnergy/Function = "Step"')
        print('dv:Tf/BeamEnergy/Times = {}'.format(n))
        for t in self.T:
            print('{:.1f} '.format(t * 1000), end='')
        print('ms')
        print('dv:Tf/BeamEnergy/Values = {}'.format(n))
        for e in self.E:
            print('{:.1f} '.format(e), end='')
        print('MeV')


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Generate TOPAS time features to generate spread-out Bragg peaks for protons in water.'
    )
    parser.add_argument('-e',
                        '--energy',
                        type=float,
                        required=True,
                        help='maximum energy in MeV')
    parser.add_argument('-c',
                        '--chi',
                        type=float,
                        required=True,
                        help='SOBP width as fraction of max. range')
    parser.add_argument(
        '-p',
        '--powerp',
        type=float,
        default=None,
        help='power law parameter p; if empty will interpolate from recommended values'
    )
    parser.add_argument('-n',
                        '--nbeams',
                        type=int,
                        default=19,
                        help='minimum number of beamlets')
    parser.add_argument('--delta',
                        type=float,
                        default=0.05,
                        help='beamlet spacing in cm')
    parser.add_argument('--recommended',
                        action='store_true',
                        help='print recommended p-value table')
    args = parser.parse_args()

    if args.recommended:
        print(SOBPSource.RecommendedValues)

    if args.powerp is None:
        args.powerp = SOBPSource.get_p(args.energy, args.chi)

    # Create object with default parameters, calculate beam weights and print them to stdout
    source = SOBPSource()
    source.generate_weights(args.energy, args.chi, args.powerp, args.nbeams,
                            args.delta)
    source.print()
