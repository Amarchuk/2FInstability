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
import scipy.optimize


class RadialToVerticalRatioHandler():

    def __init__(self, galaxy):
        self.galaxy = galaxy
        self.sigZ_to_sigR = 0.0
        self.sig_R_0 = 0.0

    def residuals(self, params, xdata, ydata):
            return (ydata - numpy.dot(xdata, params))

    def experimental_alpha_evaluation(self, normalize=False):

        x0 = [0.25, 0.25]
        sig_max = self.galaxy.sig_los_mi.bezier(0.0)**2
        points = map(lambda p : [abs(p[0]), p[1]], self.galaxy.sig_los_ma.data_points)
        points.sort()
        radii = [p[0] for p in points]

        if normalize:
            ydata = numpy.concatenate(([sig_max],[(p[1]**2)/(self.norm_sig_los_mi(p[0])**2) for p in points]))
            xdata = numpy.transpose(numpy.array([numpy.concatenate(([1.0],[self.galaxy.sve_handler.sigPhi2_to_sigR2(x) for x in radii])),
                                                 numpy.concatenate(([1.0],[1.0 for x in radii]))]))
        else:
            ydata = numpy.concatenate(([sig_max],[p[1]**2 for p in points]))
            xdata = numpy.transpose(numpy.array([numpy.concatenate(([1.0],[(self.norm_sig_los_mi(x)**2)*self.galaxy.sve_handler.sigPhi2_to_sigR2(x) for x in radii])),
                                                 numpy.concatenate(([1.0],[(self.norm_sig_los_mi(x)**2) for x in radii]))]))
        solution = scipy.optimize.leastsq(self.residuals, x0, args=(xdata, ydata))[0]
        print 'Solution: <',solution[0],' : ',solution[1],'>'
        if solution[0] > 0 and solution[1] > 0:
            tan = math.tan(self.galaxy.incl*math.pi/180.0)
            sin = math.sin(self.galaxy.incl*math.pi/180.0)
            self.sig_R_0 = math.sqrt(solution[0])/sin
            self.sigZ_to_sigR = math.sqrt(solution[1]/solution[0])*tan
            print 'sig_R_0: ', self.sig_R_0
            print 'sigZ/sigR: ', self.sigZ_to_sigR
            # self.set_sigZ_to_sigR(0.19)

    def set_sigZ_to_sigR(self, alpha):
        self.sigZ_to_sigR = 0.2
        self.sig_R_0 = self.galaxy.sig_los_mi.bezier(0.0)/math.sqrt(math.sin(self.galaxy.incl*math.pi/180.0)**2 +
                                     (self.sigZ_to_sigR*math.cos(self.galaxy.incl*math.pi/180.0))**2)


    def norm_sig_los_mi(self, x):
        sig_max = self.galaxy.sig_los_mi.bezier(0.0)
        return self.galaxy.sig_los_mi.bezier(x)/sig_max

    def plot_residuals(self):
        # plt.plot([0.0] + radii, map(abs, self.residuals((solution[0], solution[1]), xdata, ydata)), 'x-')
        pass

    def plot_sig_R(self):
        points = map(lambda p : [abs(p[0]), p[1]], self.galaxy.sig_los_ma.data_points)
        points.sort()
        radii = [p[0] for p in points]
        plt.plot(radii, [self.sig_R(x) for x in radii], 'x-', label=(r'$\sigma_{R}^{\alpha=%s}$' % self.sigZ_to_sigR))

    def sig_R(self, x):
        return self.sig_R_0*self.norm_sig_los_mi(x)

    def sig_Z(self, x):
        return self.sigZ_to_sigR*self.sig_R(x)

    def plot_sig_Z(self):
        points = map(lambda p : [abs(p[0]), p[1]], self.galaxy.sig_los_ma.data_points)
        points.sort()
        radii = [p[0] for p in points]
        plt.plot(radii, [self.sig_Z(x) for x in radii], 'x-', label=(r'$\sigma_{Z}^{\alpha=%s}$' % self.sigZ_to_sigR))

    def plot_reconstructed_sig_los_mi(self):
        points = map(lambda p : [abs(p[0]), p[1]], self.galaxy.sig_los_ma.data_points)
        points.sort()
        radii = [p[0] for p in points]
        def zero_or_positive(x):
            return 0 if x < 0 else x
        def new_sig_mi_2(x):
            return self.sig_R(x)**2*(math.sin(self.galaxy.incl*math.pi/180.0)**2 +
                                     (self.sigZ_to_sigR*math.cos(self.galaxy.incl*math.pi/180.0))**2)
        new_sig_los_mi = [math.sqrt(zero_or_positive(new_sig_mi_2(x))) for x in radii]
        plt.plot(radii, new_sig_los_mi, 'v-', label=(r'$\sigma_{mi}^{\alpha=%s}$' % self.sigZ_to_sigR))

