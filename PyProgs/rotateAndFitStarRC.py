__author__ = 'amarch'
# -*- coding: utf-8 -*-

# Отражаем звездную кривую вращения, перегибаем ее, используя данные по большой и малой оси,
# затем приблизить полиномом вручную; полученный полином и картинку сохранить если потребуется.

import matplotlib.pyplot as plt
from bisect import bisect
import math
import numpy
from scipy import interpolate
from velocityEllipsoidReconstr import *
import sys

# Класс для логирования в файл
class Tee(object):
    def __init__(self, name, mode):
        self.file = open(name, mode)
        self.stdout = sys.stdout
        sys.stdout = self

    def __del__(self):
        sys.stdout = self.stdout
        self.file.close()

    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)


def run_once(f):
    def wrapper(*args, **kwargs):
        if not globals().has_key(str(f.__name__) + '_plot'):
            globals()[str(f.__name__) + '_plot'] = True
            return f(*args, **kwargs)

    return wrapper


def meanMinorVelocity(path):
    '''Вычисляем среднюю скорость по малой оси для корректного перегибания кривой вращения по главной оси.'''
    print '#!!!!!!!!!!!!# Calculate mean velocity for the minor axis. File headers:'
    v_star_minor_file = open(path + "/v_stars_mi.dat")
    r_mi = []
    v_mi = []
    dv_mi = []
    v_mi_mean = 0
    for line in v_star_minor_file:
        if line[0] == '#':
            print line[:-1]
        else:
            line = filter(lambda x: x != '', line.split("  "))
            r_mi.append(float(line[0]))
            v_mi.append(float(line[1]))
            dv_mi.append(float(line[2]))
    v_mi_mean = array(v_mi).mean()
    print '#!!!!!!!!!!!!# v_mi min ', array(v_mi).min(), ' v_mi max ', array(v_mi).max(), ' v_mi mean ', array(
        v_mi).mean()

    @run_once
    def plot0():
        plt.figure(0)
        plt.axhline(y=v_mi_mean, color='red', label='mean')
        plt.plot(r_mi, v_mi, 'x', label='velocity minor')
        plt.errorbar(r_mi, v_mi, yerr=dv_mi, fmt=None, marker=None, mew=0)
        plt.xlabel("$r,\ ''$")
        plt.ylabel("$V_{mi},\ km/s$")
        plt.legend()
        v_star_minor_file.close()
        plt.savefig(path + "/mean_mi_v.png")

    plot0()

    return v_mi_mean


def lstqBend(r_in, v_in):
    '''Поиск лучшего порога перегибания кривой вращения с минимизацией суммы квадратов невязок.'''
    r = r_in[:]
    v = v_in[:]
    bestLSQ = inf
    bestBorder = 0
    h = 0.4 * (max(v) - min(v)) / 100
    for border in arange((min(v) * 0.7 + max(v) * 0.3), (min(v) * 0.3 + max(v) * 0.7), h):
        LSQ = 0
        greater = filter(lambda x: x[1] >= border, zip(r, v))
        greater.sort()
        smaller = filter(lambda x: x[1] < border, zip(r, v))
        for vel in smaller:
            ind = bisect(zip(*greater)[0], vel[0])
            p1 = greater[ind - 1]
            LSQ += numpy.square(p1[1] - abs(vel[1] - border) - border)
        for res in filter(lambda x: x[0] < min(zip(*smaller)[0]), greater):
            LSQ += numpy.square(res[1])
        if bestLSQ > LSQ:
            bestLSQ = LSQ
            bestBorder = border
    print "#!!!!!!!!!!!!# Best LSQ is", math.sqrt(bestLSQ)
    print "#!!!!!!!!!!!!# Best border is", bestBorder
    return bestBorder


def bendStarRC(correctGasData, correctStarData, path, incl, zero_level, manual, pol_degree_star, pol_degree_gas, name,
               scale, gascorr, monte_carlo):
    '''Перегибаем и отражаем кривую вращения вдоль большой оси, используя как середину скорость zero_level.
    После приведения к одному четверти кривая и ошибки исправляются за наклон и приближаются полиномом фиксированной
    степени. Также приближается кривая для газа. Есть возможность вручную откорректировать кривую или посчитать
    по минимальной квадратичной ошибке.'''

    print '#!!!!!!!!!!!!# Star RC rotating and fitting program for ' + name + ' start!'

    v_star_major_file = open(path + "/v_stars_ma.dat")
    r_ma1 = []
    v_ma1 = []
    dv_ma1 = []
    print '\n#!!!!!!!!!!!!# Read velocities for the major axis. File headers:'
    for line in v_star_major_file:
        if line[0] == '#':
            print line[:-1]
        else:
            line = filter(lambda x: x != '', line.split("  "))
            r_ma1.append(abs(float(line[0])) * scale)
            v_ma1.append((float(line[1])))
            dv_ma1.append(float(line[2]) / math.sin(incl * math.pi / 180))

    if monte_carlo:
        v_ma1 = [numpy.random.normal(loc=mc[0], scale=mc[1]) for mc in zip(v_ma1, dv_ma1)]

    if manual:
        correction = -60 + zero_level
    else:
        correction = lstqBend(r_ma1, v_ma1)

    print '#!!!!!!!!!!!!# Bending star on velocity level', correction

    v_ma1 = map(lambda x: (x - correction) / math.sin(incl * math.pi / 180), v_ma1)

    v_star_major_file.close()

    #Если необходимо выпрямить апроксимацию на краю - можно добавить несколько последних точек,
    #это должно помочь сгладить. Или обрезать по upperBord.

    #    upperBord = 3000
    #    r_ma, v_ma = zip(*(filter(lambda x: x[0] < upperBord, zip(r_ma, v_ma))))
    #    r_ma = list(r_ma)
    #    v_ma = list(v_ma)
    #
    #    multiplate = 5
    #    addition_points = 3
    #    r_points = heapq.nlargest(addition_points, r_ma)
    #    v_points = []
    #    dv_points = []
    #    for po in r_points:
    #        v_points.append(v_ma[r_ma.index(po)])
    #        dv_points.append(dv_ma[r_ma.index(po)])
    #    r_ma = r_ma + [i[0] + scale * i[1] for i in zip(r_points * multiplate, range(1, multiplate * addition_points + 1))]
    #    v_ma = v_ma + v_points * multiplate
    #    dv_ma = dv_ma + dv_points * multiplate

    # Также можно добавить заранее сделанные точки в нужном количестве
    #   'U11914_N7217'
    #    add_points = 50
    #    r_points = [75]
    #    v_points = [221]
    #    dv_points = [0]

    #    UGC0624_N0338
    #    add_points = 50
    #    r_points = [49]
    #    v_points = [238]
    #    dv_points = [0]

    #    U5253_N2985
    #    add_points = 70
    #    r_points = [44]
    #    v_points = [226]
    #    dv_points = [0]
    #
    #    r_ma = r_ma + [i[0] + scale * i[1] for i in zip(r_points * add_points, range(1, add_points + 1))]
    #    v_ma = v_ma + v_points * add_points
    #    dv_ma = dv_ma + dv_points * add_points

    r_ma, v_ma, dv_ma = correctStarData(r_ma1, v_ma1, dv_ma1)

    #    p_star = poly1d(polyfit(r_ma, map(abs, v_ma), deg=pol_degree_star))
    p_star = poly1d(polyfit(r_ma, map(abs, v_ma), deg=pol_degree_star, w=map(lambda x: 1 / (x + 0.1) ** 2, dv_ma)))

    #    r_ma = r_ma[:-(multiplate * addition_points)]
    #    v_ma = v_ma[:-(multiplate * addition_points)]
    #    dv_ma = dv_ma[:-(multiplate * addition_points)]

    #    r_ma = r_ma[:-add_points]
    #    v_ma = v_ma[:-add_points]
    #    dv_ma = dv_ma[:-add_points]



    print '#!!!!!!!!!!!!# Fit star RC with polinomial with degree =', pol_degree_star
    print '#!!!!!!!!!!!!# Result is '
    print p_star
    print '#!!!!!!!!!!!!# Polinomial least numpy.square error', math.sqrt(
        sum(numpy.square(p_star(r_ma) - map(abs, v_ma))))

    gas = open(path + "/v_gas_ma.dat")
    r_g1 = []
    v_g1 = []
    dv_g1 = []
    if gascorr:
        corr_by_incl = 1
    else:
        corr_by_incl = math.sin(incl * math.pi / 180)
    for line in gas:
        if line[0] == '#':
            print line[:-1]
        else:
            line = filter(lambda x: x != '', line.split(" "))
            if line.__len__() >= 4:
                r_g1.append(float(line[0]) * scale)
                v_g1.append(float(line[2]) / corr_by_incl)
                dv_g1.append(float(line[3]) / corr_by_incl)
            else:
                r_g1.append(abs(float(line[0])) * scale)
                v_g1.append(abs(float(line[1])) / corr_by_incl)
                dv_g1.append(float(line[2]) / corr_by_incl)
    gas.close()

    if monte_carlo:
        v_g1 = [numpy.random.normal(loc=mc[0], scale=mc[1]) for mc in zip(v_g1, dv_g1)]

    #    correction = lstqBend(r_g1, v_g1)
    #
    #    print '#!!!!!!!!!!!!# Bending gas on velocity level', correction
    #
    #    v_g1 = map(lambda x: abs(x - correction), v_g1)
    #    r_g1, v_g1, dv_g1 = map(list, zip(*sorted(zip(r_g1, v_g1, dv_g1))))

    #Если необходимо выпрямить апроксимацию на краю - можно добавить несколько последних точек,
    #это должно помочь сгладить. Или обрезать по upperBord.

    #    upperBord = 200
    #    r_g, v_g, dv_g = zip(*(filter(lambda x: x[0] < upperBord, zip(r_g, v_g, dv_g))))
    #    r_g = list(r_g)
    #    v_g = list(v_g)
    #    dv_g = list(dv_g)
    #
    #    multiplate = 5
    #    addition_points = 2
    #    r_points = heapq.nlargest(addition_points, r_g)
    #    v_points = []
    #    dv_points = []
    #    for po in r_points:
    #        v_points.append(v_g[r_g.index(po)])
    #        dv_points.append(dv_g[r_g.index(po)])
    #    r_g = r_g + [i[0] + scale * i[1] for i in zip(r_points * multiplate, range(1, multiplate * addition_points + 1))]
    #    v_g = v_g + v_points * multiplate
    #    dv_g = dv_g + dv_points * multiplate

    # Также можно добавить заранее сделанные точки в нужном количестве
    #   'U11914_N7217'
    #    add_points = 50
    #    r_points = [74]
    #    v_points = [302]
    #    dv_points = [0]

    #    UGC0624_N0338
    #    add_points = 30
    #    r_points = [101]
    #    v_points = [296]
    #    dv_points = [0]

    #    U5253_N2985
    #    add_points = 30
    #    r_points = [70]
    #    v_points = [245]
    #    dv_points = [0]
    #
    #
    #    r_g = r_g + [i[0] + scale * i[1] for i in zip(r_points * add_points, range(1, add_points + 1))]
    #    v_g = v_g + v_points * add_points
    #    dv_g = dv_g + dv_points * add_points

    gaslen = r_g1.__len__()

    r_g, v_g, dv_g = correctGasData(r_g1, v_g1, dv_g1)

    #    p_gas = poly1d(polyfit(r_g, v_g, deg=pol_degree_gas))
    p_gas = poly1d(polyfit(r_g, v_g, deg=pol_degree_gas, w=map(lambda x: 1 / (x + 0.1) ** 2, dv_g)))


    #    r_g = r_g[:-(multiplate * addition_points)]
    #    v_g = v_g[:-(multiplate * addition_points)]
    #    dv_g = dv_g[:-(multiplate * addition_points)]

    #    r_g = r_g[:-add_points]
    #    v_g = v_g[:-add_points]
    #    dv_g = dv_g[:-add_points]

    print '#!!!!!!!!!!!!# Fit gas with polinomial with degree =', pol_degree_gas
    print '#!!!!!!!!!!!!# Result is '
    print p_gas

    positive = filter(lambda x: x[1] >= 0, zip(r_ma, v_ma, dv_ma))
    negative = filter(lambda x: x[1] > 0, zip(r_ma, map(lambda d: (-1) * d, v_ma), dv_ma))
    positive.sort()

    fn = interpolate.interp1d(zip(*negative)[0], zip(*negative)[1])
    fp = interpolate.interp1d(zip(*positive)[0], zip(*positive)[1])
    xx = arange(max(min(zip(*negative)[0]), min(zip(*positive)[0])),
        min(max(zip(*negative)[0]), max(zip(*positive)[0])), 0.1)

    SQE = 0
    for entry in negative:
        ind = bisect(zip(*positive)[0], entry[0])
        p1 = positive[ind - 1]
        SQE += numpy.square(p1[1] - entry[1])
    for resid in filter(lambda x: x[0] < min(zip(*negative)[0]), positive):
        SQE += numpy.square(resid[1])
    print '#!!!!!!!!!!!!# Least numpy.square error after bending', math.sqrt(SQE)

    a0, b0, c0 = zip(*positive)
    a1, b1, c1 = zip(*negative)

    @run_once
    def plot1():
        plt.figure(1)
        plt.plot(a0, b0, '.', color='blue', label='upper star RC')
        plt.plot(a1, b1, '.', color='red', label='bottom star RC')
        plt.plot(arange(min(r_ma), max(r_ma) + 25, 0.1), p_star(arange(min(r_ma), max(r_ma) + 25, 0.1)),
            label='polinom approx')
        #plt.plot(xx, [abs(ii - jj) for ii, jj in zip(fn(xx), fp(xx))], '-', label='residuals bottom vs upper')
        plt.errorbar(a0, b0, yerr=c0, fmt=None, marker=None, mew=0)
        plt.errorbar(a1, b1, yerr=c1, fmt=None, marker=None, mew=0)
        plt.plot(r_g, v_g, 'x', label='gas')
        plt.plot(arange(min(r_ma), max(r_g) + 25, 0.1), p_gas(arange(min(r_ma), max(r_g) + 25, 0.1)))
        plt.errorbar(r_g, v_g, yerr=dv_g, fmt=None, marker=None, mew=0)
        plt.xlabel("$r,\ ''$")
        plt.xlim(-10)
        plt.ylim(-10, 600)
        plt.ylabel("$V_{ma},\ km/s$")
        plt.legend(loc='upper left')
        plt.savefig(path + "/bended_fitted_star_gas_RC.png")

    plot1()

    out_RC_data = open(path + "/bended_RC.dat", 'w+')
    RCdata = zip(a0 + a1, b0 + b1, c0 + c1)
    RCdata.sort()
    for entry in RCdata:
        out_RC_data.write(str(entry[0]) + ' ' + str(entry[1]) + ' ' + str(entry[2]) + '\n')
    out_RC_data.close()
    return p_star, p_gas, zip(r_ma1, map(abs, v_ma1), dv_ma1), zip(r_g[:gaslen], map(abs, v_g[:gaslen]), dv_g[:gaslen])


def bendPolinom(polinom, Rmin, Rmax, incl_from, incl_to):
    '''Для корректного изменения наклона по углу для данных.'''
    degree = polinom.order
    xx = arange(Rmin, Rmax, 0.1)
    yy = [polinom(x)*math.sin(incl_from*math.pi/180.0)/math.sin(incl_to*math.pi/180.0) for x in xx]
    return poly1d(polyfit(xx, yy, deg=degree))