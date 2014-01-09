__author__ = 'amarch'
# -*- coding: utf-8 -*-

from scipy.integrate import *

from RotationCurve import *
from Galaxy import *


class RotationCurveHandler():
    star_ma_rc = None
    star_mi_rc = None
    gas_rcs = {}

    def __init__(self, galaxyname, galaxy):
        self.name = galaxyname + " Rotation Curve Handler"
        self.galaxy = galaxy

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
        if reduce(lambda x,y: (x > 0) and (y>0), rc.radii()):
            bended_star_ma_rc.data_points  = rc.data_points
        else:
            incl_lambda = lambda x: self.incline_velocity(x, angle)
            bended_star_ma_rc.data_points = zip(map(lambda x: abs(x - zero_point[1]), rc.radii()),
                                                map(incl_lambda, map(lambda x: abs(x - zero_point[0]), rc.velocities())),
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



    # r = handler.star_ma_rc.radii()
    # dv_ma = handler.star_ma_rc.delta_velocities()
    # s1 = scipy.interpolate.splrep(handler.star_ma_rc.radii(), handler.star_ma_rc.velocities(),
    #                                                     w=map(lambda x: 1 / (x + 0.1) ** 2, dv_ma), k=1)
    # plt.plot(numpy.arange(min(r), max(r), 0.1), interpolate.splev(numpy.arange(min(r), max(r), 0.1),s1,der=0), '-', label="k=1")
    # s1 = scipy.interpolate.splrep(handler.star_ma_rc.radii(), handler.star_ma_rc.velocities(),
    #                                                     w=map(lambda x: 1 / (x + 0.1) ** 2, dv_ma), k=2)
    # plt.plot(numpy.arange(min(r), max(r), 0.1), interpolate.splev(numpy.arange(min(r), max(r), 0.1),s1,der=0), '-', label="k=2")
    # s1 = scipy.interpolate.splrep(handler.star_ma_rc.radii(), handler.star_ma_rc.velocities(),
    #                                                     w=map(lambda x: 1 / (x + 0.1) ** 2, dv_ma), k=3)
    # plt.plot(numpy.arange(min(r), max(r), 0.1), interpolate.splev(numpy.arange(min(r), max(r), 0.1),s1,der=0), '-', label="k=3")
    # s1 = scipy.interpolate.splrep(handler.star_ma_rc.radii(), handler.star_ma_rc.velocities(),
    #                                                     w=map(lambda x: 1 / (x + 0.1) ** 2, dv_ma), k=4)
    # plt.plot(numpy.arange(min(r), max(r), 0.1), interpolate.splev(numpy.arange(min(r), max(r), 0.1),s1,der=0), '-', label="k=4")
    # s1 = scipy.interpolate.splrep(handler.star_ma_rc.radii(), handler.star_ma_rc.velocities(),
    #                                                     w=map(lambda x: 1 / (x + 0.1) ** 2, dv_ma), k=5)
    # plt.plot(numpy.arange(min(r), max(r), 0.1), interpolate.splev(numpy.arange(min(r), max(r), 0.1),s1,der=0), '-', label="k=5")
    # plt.plot(handler.star_ma_rc.radii(), handler.star_ma_rc.velocities(), 'o')
    # plt.legend()
    # plt.show()
    # p_star = poly1d(polyfit(r_ma, map(abs, v_ma), deg=pol_degree_star, w=map(lambda x: 1 / (x + 0.1) ** 2, dv_ma)))

