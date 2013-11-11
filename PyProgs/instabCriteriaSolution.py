__author__ = 'amarch'
# -*- coding: utf-8 -*-

from numpy import *
import matplotlib.pyplot as plt
from scipy.optimize import *
import scipy
from scipy.special import *
import math

def run_once(f):
    def wrapper(*args, **kwargs):
        if not globals().has_key(str(f.__name__) + '_plot'):
            globals()[str(f.__name__) + '_plot'] = True
            return f(*args, **kwargs)

    return wrapper

# Нахождение критерия неустойчивости в случае двухжидкостной неустойчивости с учетом конечной толщины диска (14) и
# нет(13). Используются статьи Rafikov 2001 и Jog,Solomon 1984. Честно пишется аналитическая производная условия
# неустойчивости, соответсвующая нужному дисперсионному уравнению, после чего производная приравнивается к нулю и
# ищутся корни. Корень, дающий максимальное значение, подставляется в исходное выражение и получаем искомое Qeff.

sunR = 4.42 # Звездная величина Солнца в полосе R
sunI = 4.08 # Звездная величина Солнца в полосе I
G = 4.32 #гравитационная постоянная в нужных еденицах
sound_vel = 6  # скорость звука в газе, км/с
h_disc = 1 # экспоненциальный размер диска, исправляется ниже
I0_plot = False
# Полином приближения и граница для эпициклической частоты.
epiExpfit = poly1d([0])
epiBorder = 0

def mass_to_light(color):
    '''Отношение масса светимость вычисляется по калибровке из статьи Bell E. 2003 Table7.
    Коэффициенты берем для потока полосы R, цвет B-R.'''
    aR = -0.523
    bR = 0.683
    return power(10, aR + bR * color)


def mass_to_light_Iband(color):
    '''Отношение масса светимость вычисляется по калибровке из статьи Bell E. 2003 Table7.
    Коэффициенты берем для потока полосы I, цвет B-R.'''
    aR = -0.405
    bR = 0.518
    return power(10, aR + bR * color)


def surfaceDensityStarR(massToLight, h_disc, R, mu0_c):
    '''R полоса'''
    global I0_plot
    I0 = 4.255 * math.pow(10, 8 + 0.4 * (sunR - mu0_c))
    if not I0_plot:
        print "#!!!!!!!!!!!!# I0_R = ", massToLight * I0
        I0_plot = True
    return massToLight * I0 * math.exp(-R / h_disc)


def surfaceDensityStarI(massToLight, h_disc, R, mu0_c):
    '''I полоса'''
    global I0_plot
    I0 = 4.255 * math.pow(10, 8 + 0.4 * (sunI - mu0_c))
    if not I0_plot:
        print "#!!!!!!!!!!!!# I0_I = ", massToLight * I0
        I0_plot = True
    return massToLight * I0 * math.exp(-R / h_disc)


def surfaceDensityStarForTwoDiscs(massToLight, h_1, mu0_c_1, h_2, mu0_c_2, R):
    '''Звездная плотность в случае работы с двумя звездными дисками.'''
    global I0_plot
    I1 = 4.255 * math.pow(10, 8 + 0.4 * (sunI - mu0_c_1))
    I2 = 4.255 * math.pow(10, 8 + 0.4 * (sunI - mu0_c_2))
    if not I0_plot:
        print "#!!!!!!!!!!!!# I1 = ", massToLight * I1
        print "#!!!!!!!!!!!!# I2 = ", massToLight * I2
        I0_plot = True
    return massToLight * I1 * math.exp(-R / h_1) + massToLight * I2 * math.exp(-R / h_2)


def surfaceDensityGas(path):
    '''Возвращяет точки, в которых есть данные по газовой плотности.'''
    gas_dens = open(path + '/gas_density.dat', 'r')
    r_g = []
    gas_d = []
    for line in gas_dens:
        if line[0] == '#':
            print line[:-1]
        else:
            line = filter(lambda x: x != '', line.split("  "))
            r_g.append(float(line[0]))
            gas_d.append(1.44 * float(line[1])) # Учет молекулярного газа и гелия через фактор 1.44
    gas_dens.close()
    return r_g, gas_d


def plotSurfDens(massToLight, h_disc, mu0_c, Rmin, Rmax, step, path, surfDensFunct):
    r_g, gas_d = surfaceDensityGas(path)
    xx = arange(Rmin, Rmax, step)

    @run_once
    def plot5():
        plt.figure(5)
        plt.plot(xx, [surfDensFunct(massToLight, h_disc, x, mu0_c) for x in xx], '-')
        plt.plot(r_g, gas_d, 'x', label='gas density')
        plt.legend()
        plt.xlabel("$R,\ ''$")
        plt.ylabel(r"$\Sigma_{s}(R),\ M_{sun}/pc^2$")
        plt.savefig(path + "/surfDensity.png")

    plot5()


def plotSurfDensForTwoDiscs(massToLight, h_1, mu0_c_1, h_2, mu0_c_2, Rmin, Rmax, step, path):
    r_g, gas_d = surfaceDensityGas(path)
    xx = arange(Rmin, Rmax, step)

    @run_once
    def plot5():
        plt.figure(5)
        plt.plot(xx, [surfaceDensityStarForTwoDiscs(massToLight, h_1, mu0_c_1, h_2, mu0_c_2, x) for x in xx], '-')
        plt.plot(r_g, gas_d, 'x', label='gas density')
        plt.legend()
        plt.xlabel("$R,\ ''$")
        plt.ylabel(r"$\Sigma(R),\ M_{sun}/pc^2$")
        plt.savefig(path + "/surfDensity.png")

    plot5()


def epicyclicFreq_real(poly_gas, R, resolution):
    '''Честное вычисление эпициклической частоты на расстоянии R.'''
    return sqrt(2.0 * (poly_gas(R) ** 2) * (1 + R * poly_gas.deriv()(R) / poly_gas(R))) / (R * resolution / 1000)


def epicyclicFreq(poly_gas, R, resolution):
    '''Вычисление эпициклической частоты на расстоянии R - до 1h честно, дальше приближение.'''
    global h_disc
    global epiExpfit
    global epiBorder
    if R < h_disc:
        return sqrt(2.0 * (poly_gas(R) ** 2) * (1 + R * poly_gas.deriv()(R) / poly_gas(R))) / (R * resolution / 1000)
    else:
        if R < epiBorder:
            return sqrt(2) * poly_gas(R) / (R * resolution / 1000)
        else:
            return math.exp(epiExpfit(R))


def evalEpyciclicFreq(poly_gas, r_ma, path, resolution, h):
    '''Записываем в глобальную переменную размер диска и рисуем эпциклические частоты.'''
    global h_disc
    global epiExpfit
    global epiBorder
    h_disc = h
    kappa = [epicyclicFreq_real(poly_gas, R, resolution) for R in r_ma]
    approx = [math.sqrt(2) * poly_gas(R) / (R * resolution / 1000) for R in r_ma]
    expfit = poly1d(polyfit(r_ma, map(math.log, map(abs,approx)), deg=1))
    epiExpfit = expfit
    epiBorder = max(r_ma)
    used = [epicyclicFreq(poly_gas, R, resolution) for R in r_ma]
    expf = [math.exp(expfit(r)) for r in r_ma]


    @run_once
    def plot13():
        plt.figure(13)
        plt.plot(r_ma, kappa, 'x', label='real')
        plt.plot(r_ma, approx, '.', label='approx')
        plt.plot(r_ma, expf, 'x', label='expfit')
        plt.plot(r_ma, used, '-', label='used')
        plt.axvline(x=h_disc, ymin=0, ymax=0.05)
        plt.xlabel("$R,\ ''$")
        plt.ylabel(r"$\kappa(R),\ km/s/kpc$")
        plt.legend()
        plt.savefig(path + "/epicyclic_freq.png")

    plot13()


def Qstar(R, poly_gas, star_density, sigma, resolution):
    '''Вычисление безразмерного параметра Тумре для звездного диска. Зависит от плотности звезд, дисперсии скоростей и
    эпициклической частоты. Вычисляется по формулам на стр.4 в двухжидкостном приближении'''
    return epicyclicFreq(poly_gas, R, resolution) * sigma / (math.pi * G * star_density)


def Qgas(R, poly_gas, gas_density, resolution):
    '''Вычисление безразмерного параметра Тумре для газового диска. Зависит от плотности газа и
     эпициклической частоты. Вычисляется по формулам на стр.4 в двухжидкостном приближении'''
    return epicyclicFreq(poly_gas, R, resolution) * sound_vel / (math.pi * G * gas_density)


def dimlessWavenumber(k, R, sigma, poly_gas, resolution):
    '''Вычисление безразмерного волнового числа, где sigma - соответствующая расстоянию R дисперсия звезд.'''
    return k * sigma / epicyclicFreq(poly_gas, R, resolution)


def findTwoFluidQeffs(r_arcs, poly_gas, gas_density, star_density, sigma, path, resolution, kmax):
    '''Двухжидкостная неустойчивость. Для каждого R из r_arcs строим график 1/Qeff по формуле (13), используя
    безразмерное волновое число. Затем находим максимум через решение производной и получаем значение Qeff.'''
    Qeffs = []
    plt.figure(14)
    plt.xlabel(r"$\bar{k}$")
    plt.ylabel(r"$\frac{1}{Q_eff}$")
    plt.axhline(y=1)
    for R in r_arcs:
        print "#!!!!!!!!!!!!# ========================================="
        print "#!!!!!!!!!!!!# 2F kinetics: Compute 1/Qeff for R =", R
        ind = r_arcs.index(R)
        Qs = Qstar(R, poly_gas, star_density[ind], sigma[ind], resolution)
        Qg = Qgas(R, poly_gas, gas_density[ind], resolution)
        dimlessK = [dimlessWavenumber(k, R, sigma[ind], poly_gas, resolution) for k in arange(0.01, kmax, 0.01)]
        s = sound_vel / sigma[ind]
        print "#!!!!!!!!!!!!# Star density ", star_density[ind], " gas density ", gas_density[ind], " sigma ", sigma[
                                                                                                               ind]
        print "#!!!!!!!!!!!!# Qs ", Qs, " Qg ", Qg, " s ", s, "kappa", epicyclicFreq(poly_gas, R, resolution)
        TFcriteria = []
        root_for_max = solveDerivTwoFluidQeff(Qs, Qg, s, 0.00001, kmax)
        max_val = twoFluidQeff(Qs, Qg, s, root_for_max)
        print "#!!!!!!!!!!!!# Max 1/Qeff ", max_val, " and Qeff = ", 1 / max_val
        for barK in dimlessK:
            # Для точности используем i0e(x) = exp(-abs(x)) * i0(x)
            # Возможно использование асимптотики I_0(x) ~ exp(x)/sqrt(2pi*x)
            oneToQeff = (1 - scipy.special.i0e(barK ** 2)) * 2 / (Qs * barK)
            oneToQeff += 2 * s * barK / (Qg * (1 + (s * barK) ** 2))
            TFcriteria.append(oneToQeff)
        plt.plot(dimlessK, TFcriteria, '.', label=str(R))
        plt.plot(root_for_max, max_val, 'o')
        Qeffs.append([R, root_for_max, 1.0 / max_val])
    plt.legend()
    plt.xlim(0, 50)
    plt.savefig(path + "/qeff_2F_kinem.png")
    return Qeffs


def findTwoFluidHydroQeffs(r_arcs, poly_gas, gas_density, star_density, sigma, path, resolution, kmax):
    '''Двухжидкостная неустойчивость,гидродинамическое приближение. Для каждого R из r_arcs строим график 1/Qeff
    по формуле (9), используя безразмерное волновое число.
    Затем находим максимум через решение производной и получаем значение Qeff.'''
    Qeffs = []
    plt.figure(18)
    plt.xlabel(r"$\bar{k}$")
    plt.ylabel(r"$\frac{1}{Q_eff}$")
    plt.axhline(y=1)
    for R in r_arcs:
        print "#!!!!!!!!!!!!# ========================================="
        print "#!!!!!!!!!!!!# Hydro 2F: Compute 1/Qeff for R =", R
        ind = r_arcs.index(R)
        Qs = Qstar(R, poly_gas, star_density[ind], sigma[ind], resolution)
        Qg = Qgas(R, poly_gas, gas_density[ind], resolution)
        dimlessK = [dimlessWavenumber(k, R, sigma[ind], poly_gas, resolution) for k in arange(0.01, kmax, 0.01)]
        s = sound_vel / sigma[ind]
        print "#!!!!!!!!!!!!# Star density ", star_density[ind], " gas density ", gas_density[ind], " sigma ", sigma[
                                                                                                               ind]
        print "#!!!!!!!!!!!!# Qs ", Qs, " Qg ", Qg, " s ", s, "kappa", epicyclicFreq(poly_gas, R, resolution)
        TFcriteria = []
        root_for_max = solveDerivTwoFluidHydroQeff(Qs, Qg, s, 0.0000001, kmax)
        max_val = twoFluidHydroQeff(Qs, Qg, s, root_for_max)
        print "#!!!!!!!!!!!!# Max 1/Qeff ", max_val, " and Qeff = ", 1 / max_val
        for barK in dimlessK:
            # Для точности используем i0e(x) = exp(-abs(x)) * i0(x)
            # Возможно использование асимптотики I_0(x) ~ exp(x)/sqrt(2pi*x)
            oneToQeff = (2 * barK) / (Qs * (1 + barK ** 2))
            oneToQeff += 2 * s * barK / (Qg * (1 + (s * barK) ** 2))
            TFcriteria.append(oneToQeff)
        plt.plot(dimlessK, TFcriteria, '.', label=str(R))
        plt.plot(root_for_max, max_val, 'o')
        Qeffs.append([R, root_for_max, 1.0 / max_val])
    plt.legend()
    plt.xlim(0, 50)
    plt.savefig(path + "/qeff_2F_hydro.png")
    return Qeffs


def zGas(surfDensStar, surfDensGas, resolution):
    '''Вычисление вертикального масштаба газового диска в угловых секундах.'''
    Gz = 0.00432 # гравитационная постоянная в нужных единицах
    return (sound_vel ** 2) / (math.pi * Gz * ( surfDensStar + surfDensGas)) / resolution


def zStar(surfDensStar, surfDensGas, resolution, sigmaZ):
    '''Вычисление вертикального масштаба звездного диска в угловых секундах.'''
    Gz = 0.00432 # гравитационная постоянная в нужных единицах
    return (sigmaZ ** 2) / (math.pi * Gz * (surfDensStar + surfDensGas)) / resolution


def plotVerticalScale(surfDensStar, surfDensGas, resolution, sigmaZ, r_ma, path):
    '''Рисуем картинку с вертикальным масштабом вдоль оси галактики.'''
    hzGas = [zGas(R[1], R[2], resolution) / 2 for R in zip(r_ma, surfDensStar, surfDensGas)]
    hzStar = [zStar(R[1], R[2], resolution, R[3]) / 2 for R in zip(r_ma, surfDensStar, surfDensGas, sigmaZ)]

    @run_once
    def plot17():
        plt.figure(17)
        plt.xlabel(r"$R,''$")
        plt.ylabel(r"$h, ''$")
        plt.plot(r_ma, hzGas, 'o-', label='$h_z^g$')
        mamamax = max(hzStar)
        if mamamax > 200:
            plt.ylim(0,200)
        else:
            plt.ylim(0, mamamax)
        plt.plot(r_ma, hzStar, 'o-', label='$h_z^s$')
        plt.legend(loc='upper left')
        plt.ylim(-1)
        plt.savefig(path + "/vert_scale.png")

    plot17()


def findTwoFluidWithDiscQeffs(r_arcs, poly_gas, gas_density, star_density, sigma, path, resolution, hzStar, hzGas,
                              kmax):
    '''Двухжидкостный с диском. Для каждого R из r_arcs строим график 1/Qeff по формуле (14) измененной в
    кинематическое приближение и учитывая конечную толщину диска.
    Затем находим максимум через решение производной и получаем значение Qeff для каждого R.'''
    Qeffs = []
    plt.figure(15)
    plt.xlabel(r"$\bar{k}$")
    plt.ylabel(r"$\frac{1}{Q_eff}$")
    plt.axhline(y=1)
    for R in r_arcs:
        print "#!!!!!!!!!!!!# ========================================="
        print "#!!!!!!!!!!!!# With disc: Compute 1/Qeff for R =", R
        ind = r_arcs.index(R)
        Qs = Qstar(R, poly_gas, star_density[ind], sigma[ind], resolution)
        Qg = Qgas(R, poly_gas, gas_density[ind], resolution)
        dimlessK = [dimlessWavenumber(k, R, sigma[ind], poly_gas, resolution) for k in arange(0.01, kmax, 0.01)]
        s = sound_vel / sigma[ind]
        print "#!!!!!!!!!!!!# Star density ", star_density[ind], " gas density ", gas_density[ind], " sigma ", sigma[
                                                                                                               ind]
        print "#!!!!!!!!!!!!# Qs ", Qs, " Qg ", Qg, " s ", s, "kappa", epicyclicFreq(poly_gas, R, resolution), " hS ", hzStar[ind], " hG ", hzGas[ind]
        TFcriteria = []
        root_for_max = solveDerivTwoFluidWithDiscQeff(Qs, Qg, s, hzStar[ind], hzGas[ind],
            epicyclicFreq(poly_gas, R, resolution), sigma[ind], 0.00001, kmax)
        if root_for_max == -1:
            krange = arange(0,kmax,0.1)
            original = [twoFluidWithDiscQeffs(Qs, Qg, s, hzStar[ind], hzGas[ind], epicyclicFreq(poly_gas, R, resolution),
                sigma[ind], x) for x in krange]
            root_for_max = krange[original.index(max(original))]
        max_val = twoFluidWithDiscQeffs(Qs, Qg, s, hzStar[ind], hzGas[ind], epicyclicFreq(poly_gas, R, resolution),
            sigma[ind], root_for_max)
        print "#!!!!!!!!!!!!# Max 1/Qeff ", max_val, " and Qeff = ", 1 / max_val
        for barK in dimlessK:
            k = barK * epicyclicFreq(poly_gas, R, resolution) / sigma[ind]
            expStarMultipl = (1 - math.exp(-hzStar[ind] * k)) / (k * hzStar[ind])
            expGasMultipl = (1 - math.exp(-hzGas[ind] * k)) / (k * hzGas[ind])
            # Для точности используем i0e(x) = exp(-abs(x)) * i0(x)
            # Возможно использование асимптотики I_0(x) ~ exp(x)/sqrt(2pi*x)
            oneToQeff = expStarMultipl * (1 - scipy.special.i0e(barK ** 2)) * 2 / (Qs * barK)
            oneToQeff += expGasMultipl * 2 * s * barK / (Qg * (1 + (s * barK) ** 2))
            TFcriteria.append(oneToQeff)
        plt.plot(dimlessK, TFcriteria, '.', label=str(R))
        plt.plot(root_for_max, max_val, 'o')
        Qeffs.append([R, root_for_max, 1.0 / max_val])
    plt.legend()
    plt.xlim(0, 50)
    plt.savefig(path + "/qeff_2Fwith_disc.png")
    return Qeffs


def derivTwoFluidQeff(dimlK, Qs, Qg, s):
    '''Производная по \bar{k} от левой части (13) для того, чтобы найти максимум. Коррекция за ассимптотику производится
    с помощью встроенных функций бесселя, нормированных на exp.'''
    part1 = (1 - i0e(dimlK ** 2)) / (-dimlK ** 2)
    part2 = (2 * dimlK * i0e(dimlK ** 2) - 2 * dimlK * i1e(dimlK ** 2)) / dimlK
    part3 = (1 - (dimlK * s) ** 2) / (1 + (dimlK * s) ** 2) ** 2
    return 2 * (part1 + part2) / Qs + 2 * s * part3 / Qg


def derivTwoFluidWithDiscQeff(dimlK, Qs, Qg, s, hs, hg, kappa, sigma):
    '''Производная по \bar{k} от левой части модифицированного (14) для того, чтобы найти максимум.
    Коррекция за ассимптотику производится с помощью встроенных функций бесселя, нормированных на exp.'''
    part1d = (1 - i0e(dimlK ** 2)) / (-dimlK ** 2) + (2 * dimlK * i0e(dimlK ** 2) - 2 * dimlK * i1e(dimlK ** 2)) / dimlK
    part1d *= 2 / Qs
    part2 = 2 * (1 - i0e(dimlK ** 2)) / (dimlK * Qs)
    part3d = (1 - (dimlK * s) ** 2) / (1 + (dimlK * s) ** 2) ** 2
    part3d *= 2 * s / Qg
    part4 = 2 * s * dimlK / (1 + (dimlK * s) ** 2) / Qg
    exp1s = (1 - math.exp(-dimlK * kappa * hs / sigma)) / (dimlK * kappa * hs / sigma)
    exp1g = (1 - math.exp(-dimlK * kappa * hg / sigma)) / (dimlK * kappa * hg / sigma)
    eds = math.exp(-dimlK * kappa * hs / sigma) / dimlK - (1 - math.exp(-dimlK * kappa * hs / sigma)) / (
        (dimlK ** 2) * kappa * hs / sigma)
    edg = math.exp(-dimlK * kappa * hg / sigma) / dimlK - (1 - math.exp(-dimlK * kappa * hg / sigma)) / (
        (dimlK ** 2) * kappa * hg / sigma)
    return part2 * eds + part1d * exp1s + part4 * edg + part3d * exp1g


def solveDerivTwoFluidQeff(Qs, Qg, s, eps, kmax):
    '''Решение уравнения deriv(13) = 0 для нахождения максимума исходной функции. Запускается brentq на исходной сетке,
    в случае если на концах сетки разные знаки функции (промежуток содержит корень),
    затем выбираются лучшие корни, после чего ищется, какой их них дает максимум. Возвращается только этот корень.'''
    grid = arange(0.1, kmax, kmax / 100)
    args = [Qs, Qg, s]
    signs = [derivTwoFluidQeff(x, *args) for x in grid]
    signs = map(lambda x: x / abs(x), signs)
    roots = []
    for i in range(0, signs.__len__() - 1):
        if signs[i] * signs[i + 1] < 0:
            roots.append(brentq(lambda x: derivTwoFluidQeff(x, *args), grid[i], grid[i + 1], xtol=eps))
    original = [twoFluidQeff(Qs, Qg, s, x) for x in roots]
    return roots[original.index(max(original))]


def solveDerivTwoFluidHydroQeff(Qs, Qg, s, eps, kmax):
    '''Решение уравнения deriv(9) = 0 для нахождения максимума исходной функции. Запускается brentq на исходной сетке,
    в случае если на концах сетки разные знаки функции (промежуток содержит корень),
    затем выбираются лучшие корни, после чего ищется, какой их них дает максимум. Возвращается только этот корень.'''
    grid = arange(0.1, kmax, kmax / 100)
    args = [Qs, Qg, s]
    signs = [derivTwoFluidHydroQeff(x, *args) for x in grid]
    signs = map(lambda x: x / abs(x), signs)
    roots = []
    for i in range(0, signs.__len__() - 1):
        if signs[i] * signs[i + 1] < 0:
            roots.append(brentq(lambda x: derivTwoFluidHydroQeff(x, *args), grid[i], grid[i + 1], xtol=eps))
    original = [twoFluidHydroQeff(Qs, Qg, s, x) for x in roots]
    return roots[original.index(max(original))]


def derivTwoFluidHydroQeff(dimlK, Qs, Qg, s):
    '''Производная по \bar{k} от левой части (9) для того, чтобы найти максимум. Коррекция за ассимптотику производится
    с помощью встроенных функций бесселя, нормированных на exp.'''
    part1 = (1 - dimlK ** 2) / (1 + dimlK ** 2) ** 2
    part3 = (1 - (dimlK * s) ** 2) / (1 + (dimlK * s) ** 2) ** 2
    return (2 * part1 / Qs) + (2 * s * part3 / Qg)


def twoFluidQeff(Qs, Qg, s, dimlK):
    '''Возвращает соcчитанное значение (13).'''
    return (1 - i0e(dimlK ** 2)) * 2 / (Qs * dimlK) + 2 * s * dimlK / (Qg * (1 + (s * dimlK) ** 2))


def twoFluidHydroQeff(Qs, Qg, s, dimlK):
    '''Возвращает соcчитанное значение (9).'''
    return  2 * dimlK / (Qs * (1 + dimlK ** 2)) + 2 * s * dimlK / (Qg * (1 + (s * dimlK) ** 2))


def twoFluidWithDiscQeffs(Qs, Qg, s, hs, hg, kappa, sigma, dimlK):
    '''Возвращает соcчитанное значение (14).'''
    k = dimlK * kappa / sigma
    expStarMultipl = (1 - math.exp(-hs * k)) / (k * hs)
    expGasMultipl = (1 - math.exp(-hg * k)) / (k * hg)
    oneToQeff = expStarMultipl * (1 - scipy.special.i0e(dimlK ** 2)) * 2 / (Qs * dimlK)
    oneToQeff += expGasMultipl * 2 * s * dimlK / (Qg * (1 + (s * dimlK) ** 2))
    return oneToQeff


def solveDerivTwoFluidWithDiscQeff(Qs, Qg, s, hs, hg, kappa, sigma, eps, kmax):
    '''Решение уравнения deriv(14) = 0 для нахождения максимума исходной функции. Запускается brentq на исходной сетке,
    в случае если на концах сетки разные знаки функции (промежуток содержит корень),
    затем выбираются лучшие корни, после чего ищется, какой их них дает максимум. Возвращается только этот корень.'''
    grid = arange(0.1, kmax, kmax / 100)
    args = [Qs, Qg, s, hs, hg, kappa, sigma]
    signs = [derivTwoFluidWithDiscQeff(x, *args) for x in grid]
    signs = map(lambda x: x / abs(x), signs)
    roots = []
    for i in range(0, signs.__len__() - 1):
        if signs[i] * signs[i + 1] < 0:
            roots.append(brentq(lambda x: derivTwoFluidWithDiscQeff(x, *args), grid[i], grid[i + 1], xtol=eps))
    original = [twoFluidWithDiscQeffs(Qs, Qg, s, hs, hg, kappa, sigma, x) for x in roots]
    if original.__len__() > 0:
        return roots[original.index(max(original))]
    else:
        return -1


def findOneFluidQeffs(r_arcs, poly_gas, gas_density, star_density, sigma, path, resolution, kmax):
    '''Для каждого R из r_arcs строим график 1/Qeff для критерия Тумре (4), используя
    безразмерное волновое число. Затем находим максимум через решение производной и получаем значение Qeff.'''
    Qeffs = []
    plt.figure(16)
    plt.xlabel(r"$R,''$")
    plt.ylabel(r"$\frac{1}{Q_eff}$")
    plt.axhline(y=1)
    for R in r_arcs:
        print "#!!!!!!!!!!!!# ========================================="
        print "#!!!!!!!!!!!!# Simple one fluid : Compute 1/Qeff for R =", R
        ind = r_arcs.index(R)
        Qs = Qstar(R, poly_gas, star_density[ind], sigma[ind], resolution)
        Qg = Qgas(R, poly_gas, gas_density[ind], resolution)
        s = sound_vel / sigma[ind]
        print "#!!!!!!!!!!!!# Star density ", star_density[ind], " gas density ", gas_density[ind], " sigma ", sigma[
                                                                                                               ind]
        print "#!!!!!!!!!!!!# Qs ", Qs, " Qg ", Qg, " s ", s, "kappa", epicyclicFreq(poly_gas, R, resolution)
        Qeffs.append(Qg)
    plt.plot(r_arcs, map(lambda x: 1 / x, Qeffs), '-')
    plt.legend()
    plt.xlim(-10, 300)
    plt.savefig(path + "/qeff_1F.png")
    return Qeffs