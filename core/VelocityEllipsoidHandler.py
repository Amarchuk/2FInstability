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
from RadialToAzimuthalRatioHandler import *


class VelocityEllipsoidHandler():

    def __init__(self, galaxyname, galaxy):
        self.name = galaxyname + " Velocity Ellipsoid Handler"
        self.galaxy = galaxy
        self.sve_dict = {}
        self.sigPhi_to_sigR_ratio_handler = None


    def get_sve(self, name):
        return self.sve_dict[name]

    def sigPhi2_to_sigR2(self, R):
        if self.sigPhi_to_sigR_ratio_handler is None:
            self.sigPhi_to_sigR_ratio_handler = RadialToAzimuthalRatioHandler(self.galaxy)
            self.sigPhi_to_sigR_ratio_handler.eval_sigPhi_to_sigR()
        return self.sigPhi_to_sigR_ratio_handler.sigPhi2_to_sigR2(R)

    def plot_sigPhi2_to_sigR2(self, color='red'):
        if self.sigPhi_to_sigR_ratio_handler is None:
            self.sigPhi_to_sigR_ratio_handler = RadialToAzimuthalRatioHandler(self.galaxy)
            self.sigPhi_to_sigR_ratio_handler.eval_sigPhi_to_sigR()
        self.sigPhi_to_sigR_ratio_handler.plot_sigPhi2_to_sigR2(color)

    def sigPhi_to_sigR(self, R):
        return math.sqrt(self.sigPhi2_to_sigR2(R))

    def sigZ2_to_sigR2(self, R):
        sig_los_mi2 = 0
        if R in self.galaxy.sig_los_mi.radii():
            index = self.galaxy.sig_los_mi.radii().index(R)
            sig_los_mi2 = self.galaxy.sig_los_mi.dispersions()[index] **2
        sig_los_ma2 = self.galaxy.sig_los_ma.poly_fit(abs(R)) **2
        return (sig_los_mi2*self.sigPhi2_to_sigR2(R) - sig_los_ma2)*(math.sin(self.galaxy.incl*math.pi/180)**2)/\
               ((sig_los_ma2 - sig_los_mi2)*(math.cos(self.galaxy.incl*math.pi/180)**2))

    def plot_sigZ2_to_sigR2(self):
        xx = self.galaxy.sig_los_mi.radii()
        yy = [self.sigZ2_to_sigR2(R) for R in xx]
        plt.plot(map(abs, xx), yy, 's')
        plt.ylim(-5, 5)
        plt.axhline(y=1, color='black')
        plt.axhline(y=-1, color='black')
        plt.axhline(y=0, color='black')
        for R in map(abs, xx):
            plt.axvline(x=R)
        plt.xlabel("$R,\ arcsec$")
        plt.ylabel(r"$\sigma_{z}^2/\sigma_{R}^2$")