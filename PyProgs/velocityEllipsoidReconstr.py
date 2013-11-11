__author__ = 'amarch'
# -*- coding: utf-8 -*-

from sys import *
import os
from numpy import *
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import scipy
import math
import heapq
from scipy.interpolate import spline

def run_once(f):
    def wrapper(*args, **kwargs):
        if not globals().has_key(str(f.__name__) + '_plot'):
            globals()[str(f.__name__) + '_plot'] = True
            return f(*args, **kwargs)

    return wrapper


# Восстановление эллипсоида скоростей с использованием ассиметричного дрифта.
# Используются формулы (19)-(21) методички или (2)-(4) из Сильченко.

h_disc = 1 # экспоненциальный размер диска, исправляется ниже
# Экспоненциальный профиль SigPhi_to_sigR
sigH = 1.0
nullp = 0.5
# Экспоненциальный профиль SigLosMaj
sigLMaj = 1.0
# Экспоненциальный профиль sigR2 и sigZ2
h_kin = 0.0
sigR20 = 0.0
sigZ20 = 0.0
# Граница начала фитирования экспонентой
rexp_sigR2 = 600
rexp_sigZ2 = 600

def setStartExp(dist, name):
    '''Устанавливается граница начала фитирования экспонентой для name (варианты sigR2, sigZ2)'''
    global rexp_sigR2
    global rexp_sigZ2
    if name == 'sigR2':
        rexp_sigR2 = dist
    else:
        if name == 'sigZ2':
            rexp_sigZ2 = dist
        else:
            print '#!!!!!!!!!!!!# Wrong name parameter in setRExp\n'

def set5533Null():
    '''Для 5533 обнуляем масштабы, чтобы корректно посчитать AD.'''
    global sigZ20
    sigZ20 = 0.0


def evalStartExp(r_pc, sigma2, function):
    '''Вычисляем границу, от которой хотим начат фитировать следующим образом: находим все перегибы и смотрим, между
    какими перегибами среднее меньше.'''
    global h_disc
    global sigZ20
    # Случай 5533
    if h_disc == 34.4 and sigZ20 > 0:
        print "Is it 5533? "
        return 100.0
    else:
        approx = map(function, r_pc)
        diff = [sig - app for sig, app in zip(sigma2, approx)]
        signpoints = filter(lambda x: x[0] * x[1] < 0, zip(diff[:-1], diff[1:], r_pc[:-1]))
        signpoints = [x[2] for x in signpoints]
        signpoints.append(r_pc[0])
        signpoints.append(r_pc[-1])
        signpoints = sorted(list(set(signpoints)))
        absmed = 10000
        leftR = r_pc[0]
        data = zip(r_pc, map(abs, diff))
        for r in zip(signpoints[:-1], signpoints[1:]):
            currmed = scipy.mean([y[1] for y in filter(lambda x: x[0] > r[0] and x[0] < r[1], data)])
            if currmed < absmed:
                absmed = currmed
                leftR = r[0]
        print '#!!!!!!!!!!!!!!# leftR is ', leftR
        return leftR


def gauss(t, coeffs):
    return coeffs[0] + coeffs[1] * exp(- ((t - coeffs[2]) / coeffs[3]) ** 2)


def gauss_residuals(coeffs, y, t):
    return y - gauss(t, coeffs)


def fsech(tt, coeffs):
    return [coeffs[0] + coeffs[1] * (1 / math.cosh((t - coeffs[2]) / coeffs[3]) ** 2) for t in tt]


def fsech_residuals(coeffs, y, t):
    return [a - b for a, b in zip(y, fsech(t, coeffs))]


def evaluateSigLosWingsExpScale(path, r_eff_bulge):
    '''Обрезаем дисперсию и фитируем крылья экспонентой.'''
    global sigLMaj

    v_star_major_file = open(path + "/v_stars_ma.dat")
    r_ma = []
    sig_los_ma = []
    dsig_los_ma = []
    print '\n#!!!!!!!!!!!!# Read sig_los for the major axis.'
    for line in v_star_major_file:
        if line[0] == '#':
            pass
        else:
            line = filter(lambda x: x != '', line.split("  "))
            r_ma.append(abs(float(line[0])))
            sig_los_ma.append(float(line[3]))
            dsig_los_ma.append(float(line[4]))

    data = filter(lambda x: x[0] > r_eff_bulge, zip(r_ma, sig_los_ma, dsig_los_ma))
    r_ma1, sig_los_ma1, dsig_los_ma1 = map(list, zip(*data))

    sig_los_ma1 = map(math.log, sig_los_ma1)
    expfit = poly1d(polyfit(r_ma1, sig_los_ma1, deg=1, w=map(lambda x: 1 / (x + 0.1) ** 2, dsig_los_ma1)))
    sigLMaj = (-1 / expfit.coeffs[0])

    print '\n#!!!!!!!!!!!!# Fit sig_los_maj wing with exp, scale = ', sigLMaj, ' and const = ', math.exp(
        expfit.coeffs[1])

    @run_once
    def plot2():
        plt.figure(2)
        plt.plot(r_ma, sig_los_ma, 'o', color='blue')
        plt.errorbar(r_ma, sig_los_ma, yerr=dsig_los_ma, fmt=None, marker=None, mew=0)
        plt.plot(r_ma, map(lambda x: math.exp(-x / sigLMaj + expfit.coeffs[1]), r_ma))
        plt.xlabel("$r,\ ''$")
        plt.ylabel("$\sigma_{los}^{maj},\ km/s$")
        plt.savefig(path + "/sigma_expfit.png")

    plot2()

    v_star_major_file.close()


def evaluateAllSigLosWingsExpScale(path, r_eff_bulge):
    '''Обрезаем дисперсию и фитируем крылья экспонентой для обоих разрезов.'''
    global sigLMaj

    v_star_major_file = open(path + "/v_stars_ma.dat")
    r_ma = []
    sig_los_ma = []
    dsig_los_ma = []
    print '\n#!!!!!!!!!!!!# Read sig_los for the major axis.'
    for line in v_star_major_file:
        if line[0] == '#':
            pass
        else:
            line = filter(lambda x: x != '', line.split("  "))
            r_ma.append(abs(float(line[0])))
            sig_los_ma.append(float(line[3]))
            dsig_los_ma.append(float(line[4]))

    data = filter(lambda x: x[0] > r_eff_bulge, zip(r_ma, sig_los_ma, dsig_los_ma))
    r_ma1, sig_los_ma1, dsig_los_ma1 = map(list, zip(*data))

    sig_los_ma1 = map(math.log, sig_los_ma1)
    expfit = poly1d(polyfit(r_ma1, sig_los_ma1, deg=1, w=map(lambda x: 1 / (x + 0.1) ** 2, dsig_los_ma1)))
    sigLMaj = (-1 / expfit.coeffs[0])

    print '\n#!!!!!!!!!!!!# Fit sig_los_maj wing with exp, scale = ', sigLMaj, ' and const = ', math.exp(
        expfit.coeffs[1])

    v_star_min_file = open(path + "/v_stars_mi.dat")
    r_mi = []
    sig_los_mi = []
    dsig_los_mi = []
    for line in v_star_min_file:
        if line[0] == '#':
            pass
        else:
            line = filter(lambda x: x != '', line.split("  "))
            r_mi.append(abs(float(line[0])))
            sig_los_mi.append(float(line[3]))
            dsig_los_mi.append(float(line[4]))

    data = filter(lambda x: x[0] > r_eff_bulge, zip(r_mi, sig_los_mi, dsig_los_mi))
    r_mi1, sig_los_mi1, dsig_los_mi1 = map(list, zip(*data))

    sig_los_mi1 = map(math.log, sig_los_mi1)
    expfit_mi = poly1d(polyfit(r_mi1, sig_los_mi1, deg=1, w=map(lambda x: 1 / (x + 0.1) ** 2, dsig_los_mi1)))
    sigLMi = (-1 / expfit_mi.coeffs[0])

    print '\n#!!!!!!!!!!!!# Fit sig_los_mi wing with exp, scale = ', sigLMi, ' and const = ', math.exp(
        expfit_mi.coeffs[1])

    @run_once
    def plot2():
        plt.figure(2)
        plt.plot(r_ma, sig_los_ma, 'o', color='blue')
        plt.errorbar(r_ma, sig_los_ma, yerr=dsig_los_ma, fmt=None, marker=None, mew=0)
        plt.plot(r_ma, map(lambda x: math.exp(-x / sigLMaj + expfit.coeffs[1]), r_ma))
        plt.plot(r_mi, sig_los_mi, 'x', color='red')
        plt.errorbar(r_mi, sig_los_mi, yerr=dsig_los_mi, fmt=None, marker=None, mew=0)
        plt.plot(r_mi, map(lambda x: math.exp(-x / sigLMi + expfit_mi.coeffs[1]), r_mi))
        plt.xlabel("$r,\ ''$")
        plt.ylabel("$\sigma_{los}^{maj},\ km/s$")
        plt.savefig(path + "/sigma_expfit.png")

    plot2()

    v_star_major_file.close()
    v_star_min_file.close()


def setSigLosWingsExpScale(value):
    '''В случае если хорошо не сходится можно установить значение самому.'''
    global sigLMaj
    sigLMaj = value


def fitSechSigLosMaj(correctSigmaLosMaj, path, scale, incl):
    '''Приближение дисперсии скоростей вдоль луча зрения по большой оси sinh^2 и сохранение картинки.
    Функция возвращает параметры приближающей гауссианы.'''
    v_star_major_file = open(path + "/v_stars_ma.dat")
    r_ma1 = []
    sig_los_ma1 = []
    dsig_los_ma1 = []
    print '\n#!!!!!!!!!!!!# Read sig_los for the major axis.'
    for line in v_star_major_file:
        if line[0] == '#':
            pass
        else:
            line = filter(lambda x: x != '', line.split("  "))
            r_ma1.append(float(line[0]) * scale)
            sig_los_ma1.append(float(line[3]))
            dsig_los_ma1.append(float(line[4]))

    r_ma, sig_los_ma, dsig_los_ma, x0 = correctSigmaLosMaj(r_ma1, sig_los_ma1, dsig_los_ma1)

    # Если не сошлось - надо исправить начальное приближение ниже:
    x0 = array([0, 140, 0, 50])
    params, flag = leastsq(fsech_residuals, x0, args=(sig_los_ma, r_ma))
    print '\n#!!!!!!!!!!!!# Fit sig_los by sech^2. Fitting parameters: ', params

    @run_once
    def plot12(r_ma):
        plt.figure(12)
        plt.plot(r_ma, sig_los_ma, 'x', color='blue')
        plt.errorbar(r_ma, sig_los_ma, yerr=dsig_los_ma, fmt=None, marker=None, mew=0)
        r_ma = sorted(r_ma)
        plt.plot(r_ma, fsech(r_ma, params), label='sech^2 fit')
        plt.xlabel("$r,\ ''$")
        plt.ylabel("$\sigma_{los}^{maj},\ km/s$")
        plt.legend()
        plt.savefig(path + "/sigma_los_maj_sech2.png")

    plot12(r_ma)

    return params


def fitGaussSigLosMaj(correctSigmaLosMaj, path, scale, incl):
    '''Приближение дисперсии скоростей вдоль луча зрения по большой оси гауссианой и сохранение картинки.
    Функция возвращает параметры приближающей гауссианы.'''
    v_star_major_file = open(path + "/v_stars_ma.dat")
    r_ma1 = []
    sig_los_ma1 = []
    dsig_los_ma1 = []
    print '\n#!!!!!!!!!!!!# Read sig_los for the major axis.'
    for line in v_star_major_file:
        if line[0] == '#':
            pass
        else:
            line = filter(lambda x: x != '', line.split("  "))
            r_ma1.append(float(line[0]) * scale)
            sig_los_ma1.append(float(line[3]))
            dsig_los_ma1.append(float(line[4]))

    r_ma, sig_los_ma, dsig_los_ma, x0 = correctSigmaLosMaj(r_ma1, sig_los_ma1, dsig_los_ma1)

    params, flag = leastsq(gauss_residuals, x0, args=(sig_los_ma, r_ma))
    print '\n#!!!!!!!!!!!!# Fit sig_los by gaussian. Fitting parameters: ', params

    @run_once
    def plot3(r_ma):
        plt.figure(3)
        plt.plot(r_ma, sig_los_ma, 'x', color='blue')
        plt.errorbar(r_ma, sig_los_ma, yerr=dsig_los_ma, fmt=None, marker=None, mew=0)
        r_ma = sorted(r_ma)
        plt.plot(r_ma, gauss(r_ma, params), label='gaussian fit')
        plt.xlabel("$r,\ ''$")
        plt.ylabel("$\sigma_{los}^{maj},\ km/s$")
        plt.legend()
        plt.savefig(path + "/sigma_los_maj.png")

    plot3(r_ma)

    return params, zip(r_ma1, sig_los_ma1, dsig_los_ma1)


def fitPolySigLosMaj(correctSigmaLosMaj, path, scale, incl, poly_deg, monte_carlo, Rmax):
    '''Приближение дисперсии скоростей вдоль луча зрения по большой оси полиномом и сохранение картинки.
    Приближается перегнутая картинка, т.е. для R > 0, а рисуется для исходной отраженная положительная часть..'''
    v_star_major_file = open(path + "/v_stars_ma.dat")
    r_ma1 = []
    sig_los_ma1 = []
    dsig_los_ma1 = []
    print '\n#!!!!!!!!!!!!# Read sig_los for the major axis.'
    for line in v_star_major_file:
        if line[0] == '#':
            pass
        else:
            line = filter(lambda x: x != '', line.split("  "))
            r_ma1.append(float(line[0]) * scale)
            sig_los_ma1.append(float(line[3]))
            dsig_los_ma1.append(float(line[4]))

    if monte_carlo:
        sig_los_ma1 = [random.normal(loc=mc[0], scale=mc[1]) for mc in zip(sig_los_ma1, dsig_los_ma1)]

    r_ma, sig_los_ma, dsig_los_ma, x0 = correctSigmaLosMaj(r_ma1, sig_los_ma1, dsig_los_ma1)

    p_sig = poly1d(polyfit(map(abs, r_ma), map(abs, sig_los_ma), deg=poly_deg))
    #    p_sig = poly1d(polyfit(map(abs, r_ma), map(abs, sig_los_ma), deg=poly_deg, w=map(lambda x: 1/(x+0.1), dsig_los_ma)))
    print '\n#!!!!!!!!!!!!# Fit sig_los by polynom. Fitting parameters: '
    print p_sig

    #    r_ma = r_ma[:-(multiplate * addition_points)]
    #    sig_los_ma = sig_los_ma[:-(multiplate * addition_points)]
    #    dsig_los_ma = dsig_los_ma[:-(multiplate * addition_points)]

    @run_once
    def plot10():
        plt.figure(10)
        plt.plot(r_ma, sig_los_ma, 'x', color='blue')
        plt.errorbar(r_ma, sig_los_ma, yerr=dsig_los_ma, fmt=None, marker=None, mew=0)
        plt.plot(list(arange(-0.1, -Rmax, -0.1)[::-1]) + list(arange(0, Rmax, 0.1)),
            list(p_sig(arange(0.1, Rmax, 0.1))[::-1]) + list(p_sig(arange(0, Rmax, 0.1))), label='polinom approx')
        plt.xlabel("$r,\ ''$")
        plt.ylabel("$\sigma_{los}^{maj},\ km/s$")
        plt.ylim(0, 250)
        plt.legend()
        plt.savefig(path + "/sigma_los_maj_poly.png")

    plot10()

    return p_sig


def fitGaussSigLosMin(correctSigmaLosMin, path, scale, incl):
    '''Приближение дисперсии скоростей вдоль луча зрения по малой оси гауссианой и сохранение картинки.
    Функция возвращает параметры приближающей гауссианы. Еще исправляется за R*cos(i).'''

    radian = incl * math.pi / 180

    v_star_minor_file = open(path + "/v_stars_mi.dat")
    r_mi1 = []
    sig_los_mi1 = []
    dsig_los_mi1 = []
    print '\n#!!!!!!!!!!!!# Read sig_los for the minor axis.'
    for line in v_star_minor_file:
        if line[0] == '#':
            pass
        else:
            line = filter(lambda x: x != '', line.split("  "))
            r_mi1.append(float(line[0]))
#            r_mi1.append(float(line[0]) * scale / math.cos(radian))
            sig_los_mi1.append(float(line[3]))
            dsig_los_mi1.append(float(line[4]))

    r_mi, sig_los_mi, dsig_los_mi, x0 = correctSigmaLosMin(r_mi1, sig_los_mi1, dsig_los_mi1)

    #    # Можно обрезать в случае плохих краев
    #    r_mi = r_mi[1:-1]
    #    sig_los_mi = sig_los_mi[1:-1]
    #    dsig_los_mi  = dsig_los_mi[1:-1]

    #    # Если не сошлось - надо исправить начальное приближение гауссианы ниже:
    #    x0 = array([0, 10, 5, 10])

    params, flag = leastsq(gauss_residuals, x0, args=(sig_los_mi, r_mi))
    print '\n#!!!!!!!!!!!!# Fit sig_los_min by gaussian. Fitting parameters: ', params


#    @run_once
#    def plot8(r_mi):
#        plt.figure(8)
#        r_mi = sorted(r_mi)
#        plt.plot(r_mi, sig_los_mi, 'x', color='blue')
#        plt.errorbar(r_mi, sig_los_mi, yerr=dsig_los_mi, fmt=None, marker=None, mew=0)
#        plt.plot(r_mi, gauss(r_mi, params), label='gaussian fit')
#        plt.xlabel("$r,\ ''$")
#        plt.ylabel("$\sigma_{los}^{mi},\ km/s$")
#        plt.legend()
#        plt.savefig(path + "/sigma_los_mi.png")
#
#    plot8(r_mi)

    return params, zip(r_mi1, sig_los_mi1, dsig_los_mi1)


def fitPolySigLosMin(correctSigmaLosMin, path, scale, incl, poly_deg, monte_carlo, Rmax):
    '''Приближение дисперсии скоростей вдоль луча зрения по малой оси полиномом и сохранение картинки.
    Приближается перегнутая картинка, т.е. для R > 0, а рисуется для исходной отраженная положительная часть..'''

    radian = incl * math.pi / 180

    v_star_minor_file = open(path + "/v_stars_mi.dat")
    r_mi1 = []
    sig_los_mi1 = []
    dsig_los_mi1 = []
    print '\n#!!!!!!!!!!!!# Read sig_los for the major axis.'
    for line in v_star_minor_file:
        if line[0] == '#':
            pass
        else:
            line = filter(lambda x: x != '', line.split("  "))
            r_mi1.append(float(line[0]))
#            r_mi1.append(float(line[0]) * scale / math.cos(radian))
            sig_los_mi1.append(float(line[3]))
            dsig_los_mi1.append(float(line[4]))

            # Можно обрезать в случае плохихи краев
            #    r_mi = r_mi[1:-1]
            #    sig_los_mi = sig_los_mi[1:-1]
            #    dsig_los_mi  = dsig_los_mi[1:-1]

    if monte_carlo:
        sig_los_mi1 = [random.normal(loc=mc[0], scale=mc[1]) for mc in zip(sig_los_mi1, dsig_los_mi1)]

    r_mi, sig_los_mi, dsig_los_mi, x0 = correctSigmaLosMin(r_mi1, sig_los_mi1, dsig_los_mi1)

    p_sig = poly1d(polyfit(map(abs, r_mi), map(abs, sig_los_mi), deg=poly_deg))
    #    p_sig = poly1d(polyfit(map(abs, r_mi), map(abs, sig_los_mi), deg=poly_deg, w=map(lambda x: 1/(x+0.1), dsig_los_mi)))
    print '\n#!!!!!!!!!!!!# Fit sig_los_min by polynom. Fitting parameters: '
    print p_sig

    @run_once
    def plot11():
        plt.figure(11)
        plt.plot(r_mi, sig_los_mi, 'x', color='blue')
        plt.errorbar(r_mi, sig_los_mi, yerr=dsig_los_mi, fmt=None, marker=None, mew=0)
        plt.plot(list(arange(-0.1, -Rmax, -0.1)[::-1]) + list(arange(0, Rmax, 0.1)),
            list(p_sig(arange(0.1, Rmax, 0.1))[::-1]) + list(p_sig(arange(0, Rmax, 0.1))), label='polinom approx')
        plt.xlabel("$r,\ ''$")
        plt.ylabel("$\sigma_{los}^{min},\ km/s$")
        plt.legend()
        plt.savefig(path + "/sigma_los_min_poly.png")

    plot11()
    return p_sig


def sigPhi_to_sigR_real(poly_star, x):
    '''По формуле (20) считаем отношение sigmaPhi^2 / sigmaR^2 в точке R=x'''
    return 0.5 * (1 + x * poly_star.deriv()(x) / poly_star(x))


def sigPhi_to_sigR(poly_star, x):
    '''Сглаженное особым образом приближение (20), чтобы убрать волнистость и сделать правильную ассимптотику.
    Сделано немного странно, потому что мне было лень убирать везде зависимости.'''
    global sigH
    global nullp
    return 0.5 + nullp * math.exp(-x / sigH)


def eval_SigPhi_to_sigR(poly_star, Rmin, Rmax, step, path):
    '''Вычисление масштаба приближения SigPhi_to_sigR и отрисовка картинки.'''
    global sigH
    global nullp
    xx = arange(Rmin, Rmax, step)
    minxx = min(xx)
    maxxx = max(xx)
    xxx = filter(lambda x: x < (minxx+(maxxx-minxx) / 3) and x > 1, xx)
    yy = [sigPhi_to_sigR_real(poly_star, x) for x in xxx]
    maxyy = max(filter(lambda x: x < 1, yy))
    maxyyy = (maxyy-0.5)/math.exp(1) + 0.5
#    maxyy = 2 * maxyy / math.exp(1)
    maxyy_x = yy.index(maxyy)
    maxyy_x = xxx[maxyy_x]
    yy = zip(xxx, yy)
    intersect_list = []
    inters = 0
    for y in enumerate(yy):
        if inters == 2:
            break
        else:
            if y[0] == yy.__len__() - 1:
                break
            if y[1][1] <= maxyyy and yy[y[0] + 1][1] > maxyyy:
                intersect_list.append(y[1][0])
                sigH += y[1][0]
                inters += 1
            if y[1][1] > maxyyy and yy[y[0] + 1][1] <= maxyyy:
                intersect_list.append(y[1][0])
                sigH += y[1][0]
                inters += 1
    if inters > 0:
        sigH = sigH / inters - maxyy_x
        nullp = (maxyy-0.5)*math.exp(maxyy_x/sigH)
        print '\n#!!!!!!!!!!!!# SigPhi_to_sigR nullp = ', nullp
    else:
        #TODO: НЕ ПРОВЕРЕНО!
        expfit = poly1d(polyfit(xxx, map(math.log, [po[1] for po in yy]), deg=1))
        sigH = (-1 / expfit.coeffs[0])
        nullp = math.exp(expfit.coeffs[1])


    print '\n#!!!!!!!!!!!!# SigPhi_to_sigR exp scale = ', sigH

    @run_once
    def plot4():
        plt.figure(4)
        for y in intersect_list:
            plt.axvline(x=y, ls='--')
        plt.plot(xx, [sigPhi_to_sigR(poly_star, x) for x in xx], '.-')
        plt.plot(xx, [sigPhi_to_sigR_real(poly_star, x) for x in xx], '.-')
        plt.axvline(x=(minxx+(maxxx-minxx) / 3), ls='-')
        plt.axhline(y=(maxyyy), ls='-.')
        plt.axhline(y=0)
        plt.axhline(y=0.5)
        plt.axhline(y=1)
        plt.xlabel("$R,\ ''$")
        plt.ylabel(r"$\sigma_{\varphi}^2/\sigma_{R}^2$")
        plt.savefig(path + "/sigPhi_to_sigR.png")

    plot4()


def correctDistanceInterval(path, scale):
    '''Подсчет интервала, где есть данные и по скорости звезд и по скорости газа.'''
    gas = open(path + "/v_gas_ma.dat")
    r_g = []
    for line in gas:
        if line[0] == '#':
            pass
        else:
            line = filter(lambda x: x != '', line.split(" "))
            r_g.append(abs(float(line[0])) * scale)
    gas.close()

    bended_RC = open(path + "/bended_RC.dat")
    r_s = []
    for line in bended_RC:
        if line[0] == '#':
            pass
        else:
            line = filter(lambda x: x != '', line.split(" "))
            r_s.append(float(line[0]))
    bended_RC.close()
    Rmin = max(min(r_g), min(r_s))
    Rmax = min(max(r_g), max(r_s))
    print '#!!!!!!!!!!!!# Intersection between gas and star data = [' + str(Rmin) + ":" + str(Rmax) + "]"
    if Rmin < Rmax:
        return [Rmin, Rmax]
    else:
        return [0, 0]


def simpleSVEfromSigma(r_pc, h, path, poly_star, sigmaLosMajPoints, poly_sigma_maj, ratio, upperLimit, incl):
    '''В случе если не получается ничего вытащить из ассиметричного дрифта - можно попробовать простой алгоритм с
    фиксированным sigma_z/sigma_R = ratio.'''

    global h_kin
    global h_disc
    global sigR20
    global sigZ20
    global rexp_sigR2
    global rexp_sigZ2

    h_disc = h

    radian = incl * math.pi / 180
    upper = upperLimit
    r_gt_1h = filter(lambda x: x > h and x < upper, r_pc)

    sigPhi2_to_sigR2 = [sigPhi_to_sigR(poly_star, x) for x in r_gt_1h]
    sigLosMa = map(poly_sigma_maj, r_gt_1h)
    sig_R2 = [(x[0] ** 2) / (x[1] * math.sin(radian) ** 2 + (ratio * math.cos(radian) ** 2)) for x in
              zip(sigLosMa, sigPhi2_to_sigR2)]
    sig_Phi2 = [(x[0] * x[1]) for x in zip(sigPhi2_to_sigR2, sig_R2)]
    sig_Z2 = [x * (ratio ** 2) for x in sig_R2]

    sig_R2 = map(absSigma, sig_R2)
    expfit = poly1d(polyfit(r_gt_1h, map(math.log, sig_R2), deg=1))
    print '#!!!!!!!!!!!!# sig_R2 exp coeff ', math.exp(expfit.coeffs[1]), ' and exp scale = ', (-1 / expfit.coeffs[0])
    sigR20 = math.exp(expfit.coeffs[1])
    h_kin = (-1 / expfit.coeffs[0])
    sigZ20 = sigR20 * (ratio ** 2)

    rexp_sigR2 = evalStartExp(r_gt_1h, sig_R2, lambda x: sigR20 * math.exp(-x / h_kin))
    sig_R2_real = sig_R2[:]

    #Пересчет значений на более длинный промежуток c учетом знаний об эксп. масштабе
    sigPhi2_to_sigR2 = [sigPhi_to_sigR(poly_star, x) for x in r_pc]
    sigLosMa = map(poly_sigma_maj, r_pc)
    sig_R2 = []
    sig_Z2 = []
    sig_Phi2 = []
    for R in r_pc:
        ind = r_pc.index(R)
        if R < rexp_sigR2:
            sig_R2.append((sigLosMa[ind] ** 2) / (
            sigPhi2_to_sigR2[ind] * math.sin(radian) ** 2 + (ratio * math.cos(radian) ** 2)))
        else:
            sig_R2.append(sigR20 * math.exp(-R / h_kin))
        sig_Phi2.append(sigPhi2_to_sigR2[ind] * sig_R2[ind])
        sig_Z2.append(sig_R2[ind] * (ratio ** 2))

    rexp_sigZ2 = evalStartExp(r_pc, sig_Z2, lambda x: sigZ20 * math.exp(-x / h_kin))

    r_ma, sig_los_ma, dsig_los_ma = zip(*sigmaLosMajPoints)

    @run_once
    def plot7():
        plt.figure(7)
        #        plt.axvline(x=r_gt_1h[-sig_R2.__len__()/2])
        #        plt.plot(r_gt_1h, sig_R2, '.', color='green', label='$\sigma_{R}^2$')
        #        plt.plot(r_pc, [sigR20*math.exp(-R/h_kin) for R in r_pc], '.-', color='red', label='$expfitR2$')
        #        plt.plot(r_gt_1h, sig_Z2, 'x', color='red', label='$\sigma_{Z}^2$')
        #        plt.plot(r_pc, [sigZ20*math.exp(-R/h_kin) for R in r_pc], '.-', color='green', label='$expfitZ2$')
        #        plt.plot(r_gt_1h, sig_Phi2, 'o', color='blue', label=r'$\sigma_{\varphi}^2$')
        #        plt.plot(map(abs, r_ma), map(lambda x: x ** 2, sig_los_ma), 'o')
        #        plt.plot(r_gt_1h, [x ** 2 for x in sigLosMa], '.', color='black', label=r'$\sigma_{los,maj}^2$')
        plt.plot(r_gt_1h, sig_R2_real, '.', color='red', label='$real \sigma_{R}^2$')
        plt.plot(r_pc, sig_R2, '.', color='green', label='expfit $\sigma_{R}^2$')
        plt.plot(r_pc, sig_Z2, 'x', color='red', label='$\sigma_{Z}^2$')
        plt.plot(r_pc, sig_Phi2, 'o', color='blue', label=r'$\sigma_{\varphi}^2$')
        plt.plot(map(abs, r_ma), map(lambda x: x ** 2, sig_los_ma), 'o')
        plt.plot(r_pc, [x ** 2 for x in sigLosMa], '.', color='black', label=r'$\sigma_{los,maj}^2$')
        plt.axhline(y=0)
        mamamax = max(sig_R2 + sig_Z2)
        if mamamax > 40000:
            plt.ylim(0, 40000)
        else:
            plt.ylim(0, mamamax)
        plt.xlabel("$R,\ ''$")
        plt.ylabel("$\sigma^2,\ km/s$")
        plt.legend()
        plt.savefig(path + "/all_sigmas.png")

    plot7()

    return sig_R2, sig_Phi2, sig_Z2


def simpleSVEwhenPhiEqualsZ(r_pc, h, path, poly_star, sigmaLosMajPoints, poly_sigma_maj, ratio, upperLimit, incl):
    '''Простой алгоритм с главным предположением что sigma_z = sigma_phi. Т.о. они просто получаются равны siglaLosMaj.'''

    global h_kin
    global sigR20
    global sigZ20
    global rexp_sigR2
    global rexp_sigZ2

    radian = incl * math.pi / 180
    upper = upperLimit
    r_gt_1h = filter(lambda x: x > h and x < upper, r_pc)

    sigPhi2_to_sigR2 = [sigPhi_to_sigR(poly_star, x) for x in r_gt_1h]
    sigLosMa = map(poly_sigma_maj, r_gt_1h)
    sig_Phi2 = [x ** 2 for x in sigLosMa]
    sig_Z2 = [x ** 2 for x in sigLosMa]
    sig_R2 = [(x[0] / x[1]) for x in zip(sig_Phi2, sigPhi2_to_sigR2)]

    sig_R2 = map(absSigma, sig_R2)
    expfit = poly1d(polyfit(r_gt_1h, map(math.log, sig_R2), deg=1))
    print '#!!!!!!!!!!!!# sig_R2 exp coeff ', math.exp(expfit.coeffs[1]), ' and exp scale = ', (-1 / expfit.coeffs[0])
    sigR20 = math.exp(expfit.coeffs[1])
    h_kin = (-1 / expfit.coeffs[0])
    sigZ20 = sigR20 * (ratio ** 2)

    rexp_sigR2 = evalStartExp(r_gt_1h, sig_R2, lambda x: sigR20 * math.exp(-x / h_kin))

    sig_R2_real = sig_R2[:]

    #Пересчет значений на более длинный промежуток c учетом знаний об эксп. масштабе
    sigPhi2_to_sigR2 = [sigPhi_to_sigR(poly_star, x) for x in r_pc]
    sigLosMa = map(poly_sigma_maj, r_pc)
    sig_R2 = []
    sig_Z2 = [x ** 2 for x in sigLosMa]
    sig_Phi2 = [x ** 2 for x in sigLosMa]
    for R in r_pc:
        ind = r_pc.index(R)
        if R < rexp_sigR2:
            sig_R2.append(sig_Phi2[ind] / sigPhi2_to_sigR2[ind])
        else:
            sig_R2.append(sigR20 * math.exp(-R / h_kin))
    sig_Z2 = []
    sig_Phi2 = []
    for R in r_pc:
        ind = r_pc.index(R)
        if R < rexp_sigR2:
            sig_Z2.append(sigLosMa[ind]**2)
            sig_Phi2.append(sigLosMa[ind]**2)
        else:
            sig_Z2.append(sig_R2[ind]*sigPhi2_to_sigR2[ind])
            sig_Phi2.append(sig_R2[ind]*sigPhi2_to_sigR2[ind])

    r_ma, sig_los_ma, dsig_los_ma = zip(*sigmaLosMajPoints)

    @run_once
    def plot7():
        plt.figure(7)
        #        plt.axvline(x=r_gt_1h[-sig_R2.__len__()/2])
        #        plt.plot(r_gt_1h, sig_R2, '.', color='green', label='$\sigma_{R}^2$')
        #        plt.plot(r_pc, [sigR20*math.exp(-R/h_kin) for R in r_pc], '.-', color='red', label='$expfitR2$')
        #        plt.plot(r_gt_1h, sig_Z2, 'x', color='red', label='$\sigma_{Z}^2$')
        #        plt.plot(r_pc, [sigZ20*math.exp(-R/h_kin) for R in r_pc], '.-', color='green', label='$expfitZ2$')
        #        plt.plot(r_gt_1h, sig_Phi2, 'o', color='blue', label=r'$\sigma_{\varphi}^2$')
        #        plt.plot(map(abs, r_ma), map(lambda x: x ** 2, sig_los_ma), 'o')
        #        plt.plot(r_gt_1h, [x ** 2 for x in sigLosMa], '.', color='black', label=r'$\sigma_{los,maj}^2$')
        plt.plot(r_gt_1h, sig_R2_real, '.', color='red', label='$real \sigma_{R}^2$')
        plt.plot(r_pc, sig_R2, '.', color='green', label='$expfit\ sigma_{R}^2$')
        plt.plot(r_pc, sig_Z2, 'x', color='red', label='$\sigma_{Z}^2$')
        plt.plot(r_pc, sig_Phi2, 'o', color='blue', label=r'$\sigma_{\varphi}^2$')
        plt.plot(map(abs, r_ma), map(lambda x: x ** 2, sig_los_ma), 'o')
        plt.plot(r_pc, [x ** 2 for x in sigLosMa], '.', color='black', label=r'$\sigma_{los,maj}^2$')
        plt.axhline(y=0)
        mamamax = max(sig_R2 + sig_Z2)
        if mamamax > 40000:
            plt.ylim(0, 40000)
        else:
            plt.ylim(0, mamamax)
        plt.xlabel("$R,\ ''$")
        plt.ylabel("$\sigma^2,\ km/s$")
        plt.legend()
        plt.savefig(path + "/all_sigmas.png")

    plot7()

    return sig_R2, sig_Phi2, sig_Z2


def ratioSVEfromSigma(r_pc, h, path, poly_star, poly_sigma_maj, poly_sigma_min, upperLimit, incl):
    '''Попытка решения системы для обоих разрезов восстановить чему равно ratio в sigma_z/sigma_R = ratio.'''

    radian = incl * math.pi / 180
    upper = upperLimit
    r_gt_1h = filter(lambda x: x > h and x < upper, r_pc)
    rations = []

    sigPhi2_to_sigR2 = [sigPhi_to_sigR(poly_star, x) for x in r_gt_1h]
    sigLosMa = map(poly_sigma_maj, r_gt_1h)
    sigLosMi = map(poly_sigma_min, r_gt_1h)
    for r in r_gt_1h:
        ind = r_gt_1h.index(r)
        a = array([[math.sin(radian) ** 2, math.cos(radian) ** 2],
                   [sigPhi2_to_sigR2[ind] * math.sin(radian) ** 2, math.cos(radian) ** 2]])
        b = array([sigLosMi[ind] ** 2, sigLosMa[ind] ** 2])
        sigR2, sigR2A = linalg.solve(a, b)
        rations.append(sigR2A / sigR2)
    rations = filter(lambda x: x > 0, rations)
    rations = map(sqrt, rations)
    print '#!!!!!!!!!!!!# Mean ratio A = ' + str(scipy.mean(rations)) + 'with std = ' + str(scipy.std(rations))


def ratioZtoR(r_pc, h, path, poly_star, poly_gas, poly_sigma_maj, poly_sigma_min, upperLimit, incl):
    '''Попытка подобрать константу sigZ=A*sigR через приближение по двум разрезам.'''
    global h_kin

    radian = incl * math.pi / 180
    upper = upperLimit
    r_gt_1h = filter(lambda x: x > h and x < upper, r_pc)
    sigPhi2_to_sigR2 = [sigPhi_to_sigR(poly_star, x) for x in r_gt_1h]
    sigLosMa = map(poly_sigma_maj, r_gt_1h)
    sigLosMi = map(poly_sigma_min, r_gt_1h)
    sigLos2 = [x[0]**2 + x[1]**2 for x in zip(sigLosMa,sigLosMi)]


    def fhkin(tt, coeffs):
#        return [ simpleSigR2Evaluation(t, h_kin, radian, poly_star, poly_sigma_maj, coeffs[0])*
#                 (sigPhi_to_sigR(poly_star, t)*math.sin(radian)**2 + coeffs[0]*math.cos(radian)**2) +
#                 simpleSigR2Evaluation(t/math.cos(radian), h_kin, radian, poly_star, poly_sigma_maj, coeffs[0])*
#                 (math.sin(radian)**2 + coeffs[0]*math.cos(radian)**2) for t in tt]

        return [ ((poly_sigma_maj(t) ** 2) / ( sigPhi_to_sigR(poly_star, t) * math.sin(radian) ** 2
                                                                + ((coeffs[0]* math.cos(radian)) ** 2)))*
                         (math.sin(radian)**2 + (coeffs[0]*math.cos(radian))**2) for t in tt]

    def fhkin_residuals(coeffs, y, t):
        return [a - b for a, b in zip(y, fhkin(t, coeffs))]

    x0 = array([2.0])
    params, flag = leastsq(fhkin_residuals, x0, args=([t**2 for t in sigLosMi], r_gt_1h))
    print '#!!!!!!!!!!!!# Params 0 is ', params[0]
#    print mean(fhkin_residuals([0],[t**2 for t in sigLosMi], r_gt_1h))
#    print mean(fhkin_residuals([0.5],[t**2 for t in sigLosMi], r_gt_1h))
#    print mean(fhkin_residuals([0.8],[t**2 for t in sigLosMi], r_gt_1h))
#    print mean(fhkin_residuals([3.0],[t**2 for t in sigLosMi], r_gt_1h))

    ratio = 0.5
#    ratio = math.sqrt(params[0])
#    print '#!!!!!!!!!!!!# Ratio sigZ/sigR = ' + str(ratio)

    @run_once
    def plot8():
        plt.figure(8)
        plt.plot(r_gt_1h, [t**2 for t in sigLosMa], 'x', color='blue')
        plt.plot(r_gt_1h, [t**2 for t in sigLosMi], 'x', color='red')
#        plt.plot(r_gt_1h, [simpleSigR2Evaluation(t, h_kin, radian, poly_star, poly_sigma_maj, 0.5) for t in sigLosMi], 'x', color='green')
        plt.plot(r_gt_1h, [ ((poly_sigma_maj(t) ** 2) / ( sigPhi_to_sigR(poly_star, t) * math.sin(radian) ** 2
                                                                           + ((ratio* math.cos(radian)) ** 2)))*
                            (math.sin(radian)**2 + (ratio*math.cos(radian))**2) for t in r_gt_1h], label='min fit')
#        plt.plot(r_gt_1h, [ ((poly_sigma_maj(t) ** 2) / ( sigPhi_to_sigR_real(poly_star, t))) for t in r_gt_1h], label='mmin fit')
##        plt.plot(r_gt_1h, [simpleSigR2Evaluation(t/math.cos(radian), h_kin, radian, poly_star, poly_sigma_maj, ratio)*
##                           (math.sin(radian)**2 + (ratio**2)*math.cos(radian)**2)
##                           for t in r_gt_1h], label='min fit')
        plt.xlabel("$r,\ ''$")
        plt.ylabel("$\sigma_{los}^{2},\ km2/s2$")
        plt.legend()
        plt.savefig(path + "/ratio_sigma_los.png")

    plot8()


def asymmetricDriftEvaluation(r_pc, h, path, p_star, p_gas, upperLimit):
    '''Вычисление ассиметричного сдвига на основе формулы (21) из методички. Логарифмическая производная от радиальной
     дисперсии скоростей считается как предложено в статье Silchenko et al. 2011, экспонентой фитируется для R > 1h.
     Сами значения считаются только для тех точек, есть данные и по газу и по звездам.'''

    global h_kin
    global sigR20
    global rexp_sigR2
    global h_disc
    eps = 0.1
    h_kin = 0
    h_kin_next = h
    sigR2 = []
    upper = upperLimit
    r_gt_1h = filter(lambda x: x > h and x < upper, r_pc)
    expfit = poly1d(1)

    h_disc = h

    print '#!!!!!!!!!!!!# Asymmetric drift evaluation procedure with eps = ' + str(eps) + ' starts.'
    while(abs(h_kin - h_kin_next) > eps):
        h_kin = h_kin_next
        sigR2[:] = []
        for R in r_gt_1h:
            sigR2.append((p_gas(R) ** 2 - p_star(R) ** 2 ) / ( sigPhi_to_sigR(p_star, R) - 1 + R / h + R / h_kin ))
        #        sigR2 = map(abs, sigR2)
        sigR2 = map(math.log, sigR2)
        expfit = poly1d(polyfit(r_gt_1h, sigR2, deg=1))
        h_kin_next = (-1 / expfit.coeffs[0])
        print '#!!!!!!!!!!!!# Next approx h_kin =', h_kin_next

    h_kin = h_kin_next
    sigR2[:] = []
    for R in r_pc:
        sigR2.append((p_gas(R) ** 2 - p_star(R) ** 2 ) / ( sigPhi_to_sigR(p_star, R) - 1 + R / h + R / h_kin ))

    sigR20 = math.exp(expfit.coeffs[1])
    rexp_sigR2 = evalStartExp(r_pc, sigR2, lambda x: sigR20 * math.exp(-x / h_kin))

    @run_once
    def plot6():
        plt.figure(6)
        plt.plot(r_pc, sigR2, '.', color='green', label='points from asym drift')
        plt.plot(r_pc, map(lambda x: math.exp(expfit.coeffs[1] - x / h_kin), r_pc), '-', color='red', label='exp fit')
        plt.xlabel("$R,\ ''$")
        plt.ylabel(r"$\sigma_{R}^2,\ km^2/s^2$")
        #        plt.xlim(20,100)
        #        plt.ylim(0,90000)
        plt.legend()
        plt.savefig(path + "/sigR2.png")

    plot6()

    #    return h_kin, map(absSigma, sigR2)
    return h_kin, [sigR2Evaluation(R, h, h_kin, p_star, p_gas) for R in r_pc]


def velosityEllipsoid(h, r_pc, sigR2, path, incl, sigLosParams, poly_star):
    '''Окончательное восстановление эллипсоида скоростей по формулам (19)-(20) на заданном промежутке r_pc.
    Изображение профилей всех дисперсий на одной картинке. В качестве sigLosParams может быть как полином,
    так и параметры гауссианы. '''
    global sigZ20
    global rexp_sigZ2
    radian = incl * math.pi / 180

    sigLosMa = []
    if sigLosParams.__class__.__name__ == 'poly1d':
        sigLosMa = map(sigLosParams, r_pc)
    else:
        sigLosMa = [gauss(x, sigLosParams) for x in r_pc]

    sigZ2 = []
    sigPhi2 = []
    for R in r_pc:
        sigPhi2.append(sigPhi_to_sigR(poly_star, R) * sigR2[r_pc.index(R)])
        sigZ2.append((sigLosMa[r_pc.index(R)] ** 2 - (math.sin(radian) ** 2) * sigPhi2[-1]) / (math.cos(radian) ** 2))

    def fhkin(tt, coeffs):
        return [coeffs[0] * math.exp(-t / h_kin) for t in tt]

    def fhkin_residuals(coeffs, y, t):
        return [a - b for a, b in zip(y, fhkin(t, coeffs))]

    x0 = array([28000])
    params, flag = leastsq(fhkin_residuals, x0, args=(sigZ2[-sigZ2.__len__() / 2:], r_pc[-sigZ2.__len__() / 2:]))
    print '#!!!!!!!!!!!!# sigZ2 exp coeff ', params[0]
    sigZ20 = params[0]

    rexp_sigZ2 = evalStartExp(r_pc, sigZ2, lambda x: sigZ20 * math.exp(-x / h_kin))

    @run_once
    def plot7():
        plt.figure(7)
        #    plt.plot(r_pc, map(math.sqrt, sigR2), '.', color='green', label='$\sigma_{R}$')
        #    plt.plot(r_pc, map(math.sqrt, sigZ2), 'x', color='red', label='$\sigma_{Z}$')
        #    plt.plot(r_pc, map(math.sqrt, sigPhi2), 'o', color='blue', label=r'$\sigma_{\varphi}$')
        #    plt.plot(r_pc, [x for x in sigLosMa], '.', color='black', label=r'$\sigma_{los,maj}$')
        plt.plot(r_pc, map(lambda x: sigR20 * math.exp(-x / h_kin), r_pc), '.-', color='magenta', label='$expfitR2$')
        plt.plot(r_pc, map(lambda x: params[0] * math.exp(-x / h_kin), r_pc), '.-', color='cyan', label='$expfitZ2$')
        plt.plot(r_pc, sigR2, '.', color='green', label='$\sigma_{R}^2$')
        plt.plot(r_pc[::25], sigR2[::25], 's', color='green')
        plt.plot(r_pc, sigZ2, 'x', color='red', label='$\sigma_{Z}^2$')
        plt.plot(r_pc, sigPhi2, 'o', color='blue', label=r'$\sigma_{\varphi}^2$')
        plt.plot(r_pc, [x ** 2 for x in sigLosMa], '.', color='black', label=r'$\sigma_{los,maj}^2$')
        plt.axhline(y=0)
        mamamax = max(sigR2 + sigZ2)
        if mamamax > 40000:
            plt.ylim(0, 40000)
        else:
            plt.ylim(0, mamamax)
        plt.xlabel("$R,\ ''$")
        plt.ylabel("$\sigma^2,\ km/s$")
        plt.axvline(x=rexp_sigZ2, ymin=0, ymax=0.05)
        #plt.xlim(6000, 14000)
        #plt.ylim(0)
        plt.legend()
        plt.savefig(path + "/all_sigmas.png")

    plot7()

    sigZ2 = []
    for R in r_pc:
        if R < rexp_sigZ2:
            sigZ2.append(
                (sigLosMa[r_pc.index(R)] ** 2 - (math.sin(radian) ** 2) * sigPhi2[r_pc.index(R)]) / (
                math.cos(radian) ** 2))
        else:
            sigZ2.append(sigZ20 * math.exp(-R / h_kin))
        #        print '+++++++ r=' + str(R) + ' sigLosMa^2=' + str(sigLosMa[r_pc.index(R)] ** 2) + ' sigR^2= ' + str(sigR2[r_pc.index(R)])

    return map(absSigma, sigZ2), map(absSigma, sigPhi2)


def sigLosMiCheck(r_pc, sigLosGaussParamsMi, sigR2, sigZ2, incl, path):
    #ИСПРАВИТЬ НА ПОЛИНОМЫ!
    '''Проверка правильности восстановления эллипсоида через использование первого уравнения (19).
    Необходимо обратить внимание, что зависимость не от R, а от R*cos(i), т.е. проверка выполняется не совсем для того
    же самого интервала расстояний.'''
    radian = incl * math.pi / 180

    @run_once
    def plot9():
        plt.figure(9)
        plt.plot(r_pc, [s[0] * (math.sin(radian) ** 2) + s[1] * (math.cos(radian) ** 2) for s in zip(sigR2, sigZ2)],
            '.',
            color='green', label='$\sigma_{R}^2\sin^2(i) + \sigma_{Z}^2\cos^2(i)$')
        plt.plot(r_pc, [gauss(R, sigLosGaussParamsMi) ** 2 for R in r_pc], '.', color='black',
            label=r'$\sigma_{los,mi}^2$')
        plt.xlabel("$R,\ ''$")
        plt.ylabel("$\sigma^2,\ km/s$")
        plt.legend()
        plt.savefig(path + "/sigma_mi_check.png")

    plot9()


def sigR2Evaluation(R, h, h_kin, p_star, p_gas):
    '''Вычисление sigmaR^2 в случае, если уже известен кинетический масштаб.'''
    if R < rexp_sigR2:
        return (p_gas(R) ** 2 - p_star(R) ** 2 ) / ( sigPhi_to_sigR(p_star, R) - 1 + R / h + R / h_kin )
    else:
        return sigR20 * math.exp(-R / h_kin)

#    return (p_gas(R) ** 2 - p_star(R) ** 2 ) / ( sigPhi_to_sigR(p_star, R) - 1 + R / h + R / h_kin )

def simpleSigR2Evaluation(R, h_kin, radian, poly_star, sigLosMa, ratio):
    '''Вычисление в случае не ассим. дрифта.'''
    if R < rexp_sigR2:
      return ((sigLosMa(R) ** 2) / ( sigPhi_to_sigR(poly_star, R) * math.sin(radian) ** 2 + (ratio * math.cos(radian) ** 2)))
    else:
        return sigR20 * math.exp(-R / h_kin)


def sigZ2Evaluation(R, h, h_kin, p_star, p_gas, incl, sigLosParams ):
    '''Вычисление sigmaZ^2 для заданного R используя формулы (19)-(20).'''
    radian = incl * math.pi / 180
    sigLosMa = 0
    if sigLosParams.__class__.__name__ == 'poly1d':
        sigLosMa = sigLosParams(R)
    else:
        sigLosMa = gauss(R, sigLosParams)

    sigR2 = sigR2Evaluation(R, h, h_kin, p_star, p_gas)
    sigPhi2 = sigPhi_to_sigR(p_star, R) * sigR2
    if R < rexp_sigZ2:
        return (sigLosMa ** 2 - (math.sin(radian) ** 2) * sigPhi2) / (math.cos(radian) ** 2)
    else:
        return sigZ20 * math.exp(-R / h_kin)
    #    sigZ2 = (sigLosMa ** 2 - (math.sin(radian) ** 2) * sigPhi2) / (math.cos(radian) ** 2)
    # ПОПРАВИТЬ!

#    if sigZ2 > 0:
#        return sigZ2
#    else:
#        return 0.001

#    return sigZ2


def absSigma(sig2):
    '''Что делать, чесли в одном из значений дисперсии получилось отрицательное число.'''
    if sig2 > 0:
        return sig2
    else:
        return 0.1


def renameFilesByMethod(path, methodname):
    '''Меняем имена файлам, чтобы было проще понимать, где какой'''
    for f in os.listdir(path):
        os.rename(path+f, path+os.path.splitext(f)[0]+'_'+methodname+os.path.splitext(f)[1])