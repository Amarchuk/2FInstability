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

    def sigPhi_to_sigR(self, R):
        if self.sigPhi_to_sigR_ratio_handler is None:
            self.sigPhi_to_sigR_ratio_handler = RadialToAzimuthalRatioHandler(self.galaxy)
            self.sigPhi_to_sigR_ratio_handler.eval_sigPhi_to_sigR()
        return self.sigPhi_to_sigR_ratio_handler.sigPhi_to_sigR(R)

    def plot_sigPhi_to_sigR(self, color='red'):
        if self.sigPhi_to_sigR_ratio_handler is None:
            self.sigPhi_to_sigR_ratio_handler = RadialToAzimuthalRatioHandler(self.galaxy)
            self.sigPhi_to_sigR_ratio_handler.eval_sigPhi_to_sigR()
        self.sigPhi_to_sigR_ratio_handler.plot_sigPhi_to_sigR(color)