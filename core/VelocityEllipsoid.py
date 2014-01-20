__author__ = 'amarch'
# -*- coding: utf-8 -*-

from bisect import bisect
import sys
import numpy
from scipy import interpolate
from core.velocityEllipsoidReconstr import *
from utils import strutils as infoutils
from utils.bcolors import *
import matplotlib.pyplot as plt
from VelocityDispersion import *

class VelocityEllipsoid():

    '''SVE for short (stellar velocity ellipsoid)'''

    def __init__(self, name, description=None):
        self.name = name
        self.description = description
        self.sigR = VelocityDispersion(None, "Radial dispersions for " + name)
        self.sigPhi = VelocityDispersion(None, "Azimuthal dispersions for " + name)
        self.sigZ = VelocityDispersion(None, "Vertical dispersions for " + name)

    def print_info(self, indent):
        infoutils.print_header("Velocity Ellipsoid", self.name, self.description, indent)
        self.sigR.print_info(indent+1)
        self.sigPhi.print_info(indent+1)
        self.sigZ.print_info(indent+1)

    def plot(self, label, color='red'):
        self.sigR.plot_on_one_side(label='$\sigma_{R}$', color='blue')
        self.sigPhi.plot_on_one_side(label=r'$\sigma_{\varphi}$', color='green')
        self.sigZ.plot_on_one_side(label='$\sigma_{z}$', color='red')
        plt.legend(loc='upper right').draw_frame(False)