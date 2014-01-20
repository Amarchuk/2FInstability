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



class RotationCurveHandler():
    def __init__(self, galaxyname, galaxy):
        self.name = galaxyname + " Rotation Curve Handler"
        self.galaxy = galaxy
        self.bended_star_ma_rc = None
        self.bended_gas_ma_rc = None
        self.star_ma_rc = None
        self.star_mi_rc = None
        self.gas_rcs = {}

    def set_stellar_ma_rc(self, star_ma_rc):
        self.star_ma_rc = star_ma_rc

    def set_stellar_mi_rc(self, star_mi_rc):
        self.star_mi_rc = star_mi_rc

    def add_gas_rc(self, gas_rc, rc_name):
        if not self.gas_rcs.__contains__(rc_name):
            self.gas_rcs[rc_name] = gas_rc

    def bend_rc_integral(self, zero_point, rc):
        '''
        Simulate bending of stellar ma RC if zero is in
        point (zero_point[0],zero_point[1])
        and calculate integral as value to minimize.
        '''
        new_r = map(lambda x: x - zero_point[1], rc.radii())
        new_v = map(lambda x: x - zero_point[0], rc.velocities())
        negative = filter(lambda x: x[0] < 0 and x[1] < 0, zip(new_r, new_v))
        positive = filter(lambda x: x[0] >= 0 and x[1] >= 0, zip(new_r, new_v))
        integral = 0
        if negative.__len__() > 0:
            integral += numpy.abs(simps(zip(*negative)[1], zip(*negative)[0]))
        if positive.__len__() > 0:
            integral += numpy.abs(simps(zip(*positive)[1], zip(*positive)[0]))
        return integral

    def find_zero_point(self, rc):
        zero_point = minimize(lambda x: self.bend_rc_integral(x, rc),
                              [numpy.mean(rc.velocities()), numpy.mean(rc.radii())],
                              method='nelder-mead', options={'xtol': 1e-8, 'disp': True})
        return zero_point.x


    def incline_velocity(self, v, angle):
        return v / math.sin(angle * math.pi / 180)

    def correct_rc(self, rc, zero_point, angle):
        bended_star_ma_rc = RotationCurve(os.path.dirname(rc.path) + "/bended_" + os.path.basename(rc.path),
                                          "Bended and inclined " + rc.name)
        if self.is_corrected(rc):
            bended_star_ma_rc.data_points = rc.data_points
        else:
            incl_lambda = lambda x: self.incline_velocity(x, angle)
            bended_star_ma_rc.data_points = zip(map(lambda x: abs(x - zero_point[1]), rc.radii()),
                                                map(incl_lambda,
                                                    map(lambda x: abs(x - zero_point[0]), rc.velocities())),
                                                map(incl_lambda, rc.delta_velocities()))
            bended_star_ma_rc.data_points.sort()
        return bended_star_ma_rc


    def get_corrected_star_ma_rc(self, zero_point, angle):
        self.bended_star_ma_rc = self.correct_rc(self.star_ma_rc, zero_point, angle)
        return self.bended_star_ma_rc

    def get_corrected_gas_ma_rc(self, rc_name, zero_point, angle):
        self.bended_gas_ma_rc = self.correct_rc(self.gas_rcs[rc_name], zero_point, angle)
        return self.bended_gas_ma_rc

    def interpolate_poly_rc(self, rc, deg):
        rc.poly_fit = poly1d(polyfit(zip(*(rc.data_points + rc.fake_data_points))[0],
                                     zip(*(rc.data_points + rc.fake_data_points))[1], deg=deg,
                                     w=map(lambda x: 1 / (x + 0.1) ** 2,
                                           zip(*(rc.data_points + rc.fake_data_points))[2])))

    def get_all_rcs(self):
        all_rcs = []
        if self.star_ma_rc is not None:
            all_rcs.append(self.star_ma_rc)
        if self.star_mi_rc is not None:
            all_rcs.append(self.star_mi_rc)
        if self.bended_star_ma_rc is not None:
            all_rcs.append(self.bended_star_ma_rc)
        if self.bended_gas_ma_rc is not None:
            all_rcs.append(self.bended_gas_ma_rc)
        for gas_rc_name in self.gas_rcs:
            if self.gas_rcs[gas_rc_name] is not None:
                all_rcs.append(self.gas_rcs[gas_rc_name])
        return all_rcs

    def all(self, iterable):
        for element in iterable:
            if not element:
                return False
        return True

    def is_corrected(self, rc):
        return self.all(i >= 0 for i in rc.radii()) and \
               self.all(i < 500 for i in rc.velocities())


    def print_info(self, indent):
        infoutils.print_header("Rotation Curve Handler", self.name, indent=indent)
        gas_corrected_and_fitted = True if (
        self.bended_gas_ma_rc is not None and self.bended_gas_ma_rc.poly_fit is not poly1d([0])) else False
        star_corrected_and_fitted = True if (
        self.bended_star_ma_rc is not None and self.bended_star_ma_rc.poly_fit is not poly1d([0])) else False
        infoutils.print_simple_param(indent + 1, "Star RC corected and fitted", str(star_corrected_and_fitted))
        infoutils.print_simple_param(indent + 1, "Gas RC corected and fitted", str(gas_corrected_and_fitted))
        all_rcs = self.get_all_rcs()
        for rc in all_rcs:
            rc.print_info(indent + 1)

    def centralize_rc(self, rc):
        '''
        subtract mean value from rc velocities
        '''
        mean = numpy.mean(rc.velocities())
        centr_rc = copy.deepcopy(rc)
        centr_rc.data_points = zip(rc.radii(), map(lambda x: x-mean, rc.velocities()), rc.delta_velocities())
        return centr_rc

    def plot_all_uncorrected_centralize(self):
        plt.title("All Not corrected rotation curves for " + self.name)
        colors = ['r', 'b', 'g', 'm', 'k', 'y', 'c']
        all_rcs = self.get_all_rcs()
        if self.bended_gas_ma_rc is not None:
            all_rcs.remove(self.bended_gas_ma_rc)
        if self.bended_star_ma_rc is not None:
            all_rcs.remove(self.bended_star_ma_rc)
        for rc_and_col in zip(all_rcs, itertools.cycle(colors)):
            rc = rc_and_col[0]
            if not self.is_corrected(rc):
                rc = self.centralize_rc(rc_and_col[0])
            color = rc_and_col[1]
            rc.plot(rc.name, color=color)

    def rcs_radii_intersection(self, rc1, rc2):
        Rmin = max(min(rc1.radii()), min(rc2.radii()))
        Rmax = min(max(rc1.radii()), max(rc2.radii()))
        if Rmin < Rmax:
            return [Rmin, Rmax]
        else:
            return [0, 0]