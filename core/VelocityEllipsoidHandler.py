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

    def build_sigR_bezier(self, R):
        sig_los_ma2 = self.galaxy.sig_los_ma.bezier(R)**2
        sig_los_mi2 = self.galaxy.sig_los_mi.bezier(R)**2
        result = (sig_los_mi2 - sig_los_ma2)/((1-self.sigPhi2_to_sigR2(R))*math.sin(self.galaxy.incl*math.pi/180)**2)
        return 0 if result < 0 else math.sqrt(result)

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

    def plot_bezier_sigR(self):
        xx = map(abs, self.galaxy.sig_los_mi.radii())
        xx.sort()
        yy = [self.build_sigR_bezier(R) for R in xx]
        plt.plot(xx, yy, 's-', label = r"$\sigma_{R}$")
        if max(yy) > 250:
            plt.ylim(0, 250)
        elif min(yy) < 0:
            plt.ylim(0)

    def plot_bezier_sigPhi(self):
        xx = map(abs, self.galaxy.sig_los_mi.radii())
        xx.sort()
        yy = [self.build_sigR_bezier(R)*self.sigPhi_to_sigR(R) for R in xx]
        plt.plot(xx, yy, 'p-', label = r"$\sigma_{\phi}$")

    def zero_or_positive(self, x):
        return 0 if x < 0 else x

    def plot_sigZ_mi(self, sig_R):
        xx = map(abs, self.galaxy.sig_los_mi.radii())
        xx.sort()
        def sig_z2_mi(R, sig_R):
            return (self.galaxy.sig_los_mi.bezier(R)**2 - (sig_R(R)*math.sin(self.galaxy.incl*math.pi/180))**2)/\
                   (math.cos(self.galaxy.incl*math.pi/180))
        yy = [math.sqrt(self.zero_or_positive(sig_z2_mi(R, sig_R))) for R in xx]
        plt.plot(xx, yy, 'd-', label = r"$\sigma_{z}^{mi}$")

    def plot_sigZ_ma(self, sig_Phi):
        xx = map(abs, self.galaxy.sig_los_ma.radii())
        xx.sort()
        def sig_z2_ma(R, sig_Phi):
            return (self.galaxy.sig_los_ma.bezier(R)**2 - (sig_Phi(R)*math.sin(self.galaxy.incl*math.pi/180))**2)/\
                   (math.cos(self.galaxy.incl*math.pi/180))
        yy = [math.sqrt(self.zero_or_positive(sig_z2_ma(R, sig_Phi))) for R in xx]
        plt.plot(xx, yy, 'D-', label = r"$\sigma_{z}^{ma}$")
        
    def plot_sve_from_bezier(self):
        self.plot_bezier_sigR()
        self.plot_bezier_sigPhi()
        self.plot_sigZ_mi(self.build_sigR_bezier)
        self.plot_sigZ_ma(lambda x: self.build_sigR_bezier(x)*self.sigPhi_to_sigR(x))
        self.galaxy.sig_los_ma.plot_on_one_side('$\sigma_{los}^{maj}$', all_data=True)
        self.galaxy.sig_los_mi.plot_on_one_side('$\sigma_{los}^{min}$', color='blue', all_data=True)
        plt.legend()

    def experimental_alpha_evaluation(self, normalize=False):

        def residuals(params, xdata, ydata):
            return (ydata - numpy.dot(xdata, params))

        x0 = [0.25, 0.25]
        sig_max = self.galaxy.sig_los_mi.bezier(0.0)**2
        points = map(lambda p : [abs(p[0]), p[1]], self.galaxy.sig_los_ma.data_points)
        points.sort()
        radii = [p[0] for p in points]

        def F(x):
            return (self.galaxy.sig_los_mi.bezier(x)**2)/sig_max

        if normalize:
            ydata = numpy.concatenate(([sig_max],[(p[1]**2)/F(p[0]) for p in points]))
            xdata = numpy.transpose(numpy.array([numpy.concatenate(([1.0],[self.sigPhi2_to_sigR2(x) for x in radii])),
                                                 numpy.concatenate(([1.0],[1.0 for x in radii]))]))
        else:
            ydata = numpy.concatenate(([sig_max],[p[1]**2 for p in points]))
            xdata = numpy.transpose(numpy.array([numpy.concatenate(([1.0],[F(x)*self.sigPhi2_to_sigR2(x) for x in radii])),
                                                 numpy.concatenate(([1.0],[F(x) for x in radii]))]))
        solution = scipy.optimize.leastsq(residuals, x0, args=(xdata, ydata))[0]
        print '<',solution[0],' : ',solution[1],'>'
        if solution[0] > 0 and solution[1] > 0:
            tan = math.tan(self.galaxy.incl*math.pi/180.0)
            sin = math.sin(self.galaxy.incl*math.pi/180.0)
            sig_R_0 = math.sqrt(solution[0])/sin
            print 'sig_R_0: ', sig_R_0
            print 'sig_Z/sigR: ', math.sqrt(solution[1]/solution[0])*tan
        plt.plot([0.0] + radii, map(abs, residuals((solution[0], solution[1]), xdata, ydata)), 'x-')