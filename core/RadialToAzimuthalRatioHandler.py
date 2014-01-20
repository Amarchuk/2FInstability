__author__ = 'amarch'
# -*- coding: utf-8 -*-

from utils import strutils as infoutils
import itertools

from scipy.integrate import *

from RotationCurve import *
from Galaxy import *
from utils import strutils as infoutils
import itertools
import copy


class RadialToAzimuthalRatioHandler():

    def __init__(self, galaxy):
        self.galaxy = galaxy

    def eval_sigPhi_to_sigR(self):

        self.sigPhi_to_sigR_expscale = 0
        self.sigPhi_to_sigR_nullp = 0

        self.poly_star = self.galaxy.star_rc.poly_fit
        (Rmin, Rmax) = self.galaxy.rc_handler.rcs_radii_intersection(self.galaxy.star_rc, self.galaxy.gas_rc)
        step = (Rmax - Rmin)/1000
        self.xx = arange(Rmin, Rmax, step)
        self.minxx = min(self.xx)
        self.maxxx = max(self.xx)
        xxx = filter(lambda x: x < (self.minxx+(self.maxxx-self.minxx) / 3) and x > 1, self.xx)
        yy = [sigPhi_to_sigR_real(self.poly_star, x) for x in xxx]
        maxyy = max(filter(lambda x: x < 1, yy))
        self.maxyyy = (maxyy-0.5)/math.exp(1) + 0.5
        maxyy_x = yy.index(maxyy)
        maxyy_x = xxx[maxyy_x]
        yy = zip(xxx, yy)
        self.intersect_list = []
        inters = 0
        for y in enumerate(yy):
            if inters == 2:
                break
            else:
                if y[0] == yy.__len__() - 1:
                    break
                if y[1][1] <= self.maxyyy and yy[y[0] + 1][1] > self.maxyyy:
                    self.intersect_list.append(y[1][0])
                    self.sigPhi_to_sigR_expscale += y[1][0]
                    inters += 1
                if y[1][1] > self.maxyyy and yy[y[0] + 1][1] <= self.maxyyy:
                    self.intersect_list.append(y[1][0])
                    self.sigPhi_to_sigR_expscale += y[1][0]
                    inters += 1
        if inters > 0:
            self.sigPhi_to_sigR_expscale = self.sigPhi_to_sigR_expscale / inters - maxyy_x
            self.sigPhi_to_sigR_nullp = (maxyy-0.5)*math.exp(maxyy_x/self.sigPhi_to_sigR_expscale)
        else:
            #TODO: НЕ ПРОВЕРЕНО!
            expfit = poly1d(polyfit(xxx, map(math.log, [po[1] for po in yy]), deg=1))
            self.sigPhi_to_sigR_expscale = (-1 / expfit.coeffs[0])
            self.sigPhi_to_sigR_nullp = math.exp(expfit.coeffs[1])




    def plot_sigPhi2_to_sigR2(self, color='red'):
        for y in self.intersect_list:
            plt.axvline(x=y, ls='--')
        plt.plot(self.xx, [self.sigPhi2_to_sigR2(x) for x in self.xx], '.-', color=color)
        plt.plot(self.xx, [self.sigPhi_to_sigR_real(x) for x in self.xx], '.-')
        plt.axvline(x=(self.minxx+(self.maxxx-self.minxx) / 3), ls='-')
        plt.axhline(y=(self.maxyyy), ls='-.')
        plt.axhline(y=0)
        plt.axhline(y=0.5)
        plt.axhline(y=1)
        plt.xlabel("$R,\ arcsec$")
        plt.ylabel(r"$\sigma_{\varphi}^2/\sigma_{R}^2$")

    def sigPhi_to_sigR_real(self, R):
        return 0.5 * (1 + R * self.poly_star.deriv()(R) / self.poly_star(R))


    def sigPhi2_to_sigR2(self, R):
        return 0.5 + self.sigPhi_to_sigR_nullp * math.exp(-R/self.sigPhi_to_sigR_expscale)