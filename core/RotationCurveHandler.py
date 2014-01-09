__author__ = 'amarch'
# -*- coding: utf-8 -*-

from scipy.integrate import *

from RotationCurve import *
from Galaxy import *


class RotationCurveHandler():
    def __init__(self, galaxyname, galaxy):
        self.star_ma_rc = None
        self.star_mi_rc = None
        self.gas_rcs = {}
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
        zero_point = minimize(self.bend_rc_integral,
                              [numpy.mean(rc.velocities()), numpy.mean(rc.radii())],
                              method='nelder-mead', options={'xtol': 1e-8, 'disp': True})
        return zero_point.x

    def incline_velocity(self, v, angle):
        return v / math.sin(angle * math.pi / 180)

    def bend_star_ma_rc(self):
        rc = self.star_ma_rc
        zero_point = self.find_zero_point(rc)
        self.bended_star_ma_rc = RotationCurve(os.path.basename(self.star_ma_rc.path) + "/v_star_ma_bended.dat",
                                               "Bended RC: " + self.star_ma_rc.name)
        incl_lambda = lambda x: self.incline_velocity(self.galaxy.incl, x)
        self.bended_star_ma_rc.data_points = zip(abs(self.star_ma_rc.radii() - zero_point[1]),
                                                 map(incl_lambda, abs(self.star_ma_rc.velocities() - zero_point[0])),
                                                 map(incl_lambda, self.star_ma_rc.delta_velocities()))
        self.bended_star_ma_rc.data_points.sort()
        return self.bended_star_ma_rc

    def bend_gas_rc(self, rc_name):
    # self.gas_rc = RotationCurve(os.path.basename(self.gas_rcs[rc_name]) + "/v_gas.dat",
    #                                        "Main gas RC: " + self.gas_rcs[rc_name].name)
    # incl_lambda = lambda x: self.incline_velocity(self.galaxy.incl, x)
    # self.gas_rc.data_points = zip(self.star_ma_rc.radii(),
    #                                          map(incl_lambda, abs(self.star_ma_rc.velocities() - zero_point[0])),
    #                                          map(incl_lambda, self.star_ma_rc.delta_velocities()))


if __name__ == "__main__":
    ngc338 = Galaxy(name="NGC 338 (UGC 624)", path="../data/ngc338", incl="64.0", delta_incl="7.5",
                    description="Inclination according to Zasov 2012, photometry in I band.",
                    resolution="311", image="../data/ngc338/ngc338_SDSS.jpeg")
    handler = RotationCurveHandler("ngc 338", ngc338)
    handler.set_stellar_ma_rc(RotationCurve("../data/ngc338/v_stars_ma.dat", "Stars MA RC NGC338"))
    zero = handler.find_zero_point()
    print zero


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

