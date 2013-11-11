__author__ = 'amarch'
# -*- coding: utf-8 -*-

from velocityEllipsoidReconstr import *
from rotateAndFitStarRC import *
from instabCriteriaSolution import *
from plotFinal import *
import sys
import time
import os
import shutil
from main import *



def correctGasData(r_g1, v_g1, dv_g1):
    '''Функция, куда убраны все операции подгонки с данными по газу.'''

    r_g = r_g1
    v_g = v_g1
    dv_g = dv_g1

    #Если необходимо выпрямить апроксимацию на краю - можно добавить несколько последних точек,
    #это должно помочь сгладить. Или обрезать по upperBord.

    #    upperBord = 200
    #    r_g, v_g, dv_g = zip(*(filter(lambda x: x[0] < upperBord, zip(r_g, v_g, dv_g))))
    #    r_g = list(r_g)
    #    v_g = list(v_g)
    #    dv_g = list(dv_g)

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


    #    correction = lstqBend(r_g, v_g)
#    correction = 4952/math.sin(36*math.pi/180)
#    v_g = map(lambda x: abs(x - correction), v_g)
#    r_g, v_g, dv_g = map(list, zip(*sorted(zip(r_g, v_g, dv_g))))

    #    add_points = 5
    #    r_points = [32]
    #    v_points = [285]
    #    dv_points = [1]
    #
    #    r_g = r_g + [i[0] + scale * i[1] for i in zip(r_points * add_points, range(1, add_points + 1))]
    #    v_g = v_g + v_points * add_points
    #    dv_g = dv_g + dv_points * add_points
    #
    #    add_points = 54
    #    r_points = [46]
    #    v_points = [268]
    #    dv_points = [1]
    #
    #    r_g = r_g + [i[0] + scale * i[1] for i in zip(r_points * add_points, range(1, add_points + 1))]
    #    v_g = v_g + v_points * add_points
    #    dv_g = dv_g + dv_points * add_points

    return r_g, v_g, dv_g


def correctStarData(r_ma1, v_ma1, dv_ma1):
    '''Корректировка данных по звездам.'''

    r_ma = r_ma1
    v_ma = v_ma1
    dv_ma = dv_ma1

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

    add_points = 70
    r_points = [36]
    v_points = [190]
    dv_points = [2]

    r_ma = r_ma + [i[0] + scale * i[1] for i in zip(r_points * add_points, range(1, add_points + 1))]
    v_ma = v_ma + v_points * add_points
    dv_ma = dv_ma + dv_points * add_points

    return r_ma, v_ma, dv_ma


def correctSigmaLosMaj(r_ma1, sig_los_ma1, dsig_los_ma1):
    '''Корректируем данные по дисперсии скоростей вдоль главной оси. '''

    # Если не сошлось - надо исправить начальное приближение гауссианы ниже:
    x0 = array([0, 100, 5, 100])

    # на случай если данные из разных источников в одном файле
    r_ma, sig_los_ma, dsig_los_ma = map(list, zip(*sorted(zip(r_ma1, sig_los_ma1, dsig_los_ma1))))

    # Можно обрезать в случае плохих краев
    #    r_ma = r_ma[1:-1]
    #    sig_los_ma = sig_los_ma[1:-1]
    #    dsig_los_ma  = dsig_los_ma[1:-1]

    #    #Если необходимо выпрямить апроксимацию на краю - можно добавить несколько последних точек,
    #    #это должно помочь сгладить.
    #
    #    multiplate = 10
    #    addition_points = 1
    #    r_points = heapq.nlargest(addition_points, r_ma)
    #    sig_points = []
    #    dsig_points = []
    #    for po in r_points:
    #        sig_points.append(sig_los_ma[r_ma.index(po)])
    #        dsig_points.append(dsig_los_ma[r_ma.index(po)])
    #    r_ma = r_ma + [i[0] + scale * i[1] for i in
    #                   zip(r_points * multiplate, arange(1, 3 * (multiplate * addition_points) + 1, 3))]
    #    sig_los_ma = sig_los_ma + sig_points * multiplate
    #    dsig_los_ma = dsig_los_ma + dsig_points * multiplate

    add_points = 87
    r_points = [20.4]
    v_points = [238]
    dv_points = [1]

    #    Экспоненциальные точки
    r_ma = r_ma + [i[0] + i[1] for i in zip(r_points * add_points, arange(1, add_points + 1, 1))]
    sig_los_ma = sig_los_ma + [60 * math.exp(-x / 55.0) for x in
                               [i[0] + i[1] for i in zip(r_points * add_points, arange(1, add_points + 1, 1))]]
    dsig_los_ma = dsig_los_ma + dv_points * add_points

    return r_ma, sig_los_ma, dsig_los_ma, x0


def correctSigmaLosMin(r_ma1, sig_los_ma1, dsig_los_ma1):
    '''Корректируем данные по дисперсии скоростей вдоль главной оси. '''

    r_ma, sig_los_ma, dsig_los_ma = map(list, zip(*sorted(zip(r_ma1, sig_los_ma1, dsig_los_ma1))))

    # Можно обрезать в случае плохих краев
    #    r_ma = r_ma[1:-1]
    #    sig_los_ma = sig_los_ma[1:-1]
    #    dsig_los_ma = dsig_los_ma[1:-1]

    # Если не сошлось - надо исправить начальное приближение гауссианы ниже:
    x0 = array([0, 10, 5, 10])

    #Если необходимо выпрямить апроксимацию на краю - можно добавить несколько последних точек,
    #это должно помочь сгладить.
    #    multiplate = 10
    #    addition_points = 1
    #    r_points = heapq.nlargest(addition_points, r_ma)
    #    sig_points = []
    #    dsig_points = []
    #    for po in r_points:
    #        sig_points.append(sig_los_ma[r_ma.index(po)])
    #        dsig_points.append(dsig_los_ma[r_ma.index(po)])
    #    r_ma = r_ma + [i[0] + scale * i[1] for i in
    #                   zip(r_points * multiplate, arange(1, 5 * (multiplate * addition_points) + 1, 5))]
    #    sig_los_ma = sig_los_ma + sig_points * multiplate
    #    dsig_los_ma = dsig_los_ma + dsig_points * multiplate
    add_points = 87
    r_points = [28.54]
    v_points = [238]
    dv_points = [1]

    #    Экспоненциальные точки
    r_ma = r_ma + [i[0] + i[1] for i in zip(r_points * add_points, arange(1, add_points + 1, 1))]
    sig_los_ma = sig_los_ma + [65 * math.exp(-x / 187.0) for x in
                               [i[0] + i[1] for i in zip(r_points * add_points, arange(1, add_points + 1, 1))]]
    dsig_los_ma = dsig_los_ma + dv_points * add_points


    return r_ma, sig_los_ma, dsig_los_ma, x0


startTime = time.time()


if __name__ == "__main__":
    plt.rcParams.update({'font.size': 16})

    path = '/home/amarch/Documents/RotationCurves/Diploma/TwoFluidInstAllDataFromSotn17Feb/Sample/RC/U3546_N2273'
    name = 'U3546_N2273'
    incl = 55
    scale = 1
    resolution = 130 #pc/arcsec
    h_disc = 21.1  # R-band
    M_R = 11.13
    M_B = 12.56
    mu0_c_R = 19.49
    r_eff_bulge = 2.2
    pol_degree_star = 15
    pol_degree_gas = 8
    sig_pol_deg = 10
    sig_pol_deg_mi = 10
    Rmin = 10
    Rmax = 106
    gas_corr_by_incl = True
    M_to_L = mass_to_light(M_B - M_R)
    di = 3
    monte_carlo_realizations = 1
    peculiarities = [90,100]
    maxDisc = 2.5
    sig_wings = 25.0 # откуда крылья для дисперсий фитировать
    use_minor = True # используется ли дисперсия по малой оси

    if not os.path.exists(path+'/EQUAL_BELL/'):
        os.makedirs(path+'/EQUAL_BELL/')
    else:
        for f in os.listdir(path+'/EQUAL_BELL/'):
            os.remove(path+'/EQUAL_BELL/'+f)
    shutil.copy2(path+'/v_stars_ma.dat', path+'/EQUAL_BELL/v_stars_ma.dat')
    shutil.copy2(path+'/v_gas_ma.dat', path+'/EQUAL_BELL/v_gas_ma.dat')
    shutil.copy2(path+'/gas_density.dat', path+'/EQUAL_BELL/gas_density.dat')
    if os.path.exists(path+'/v_stars_mi.dat'):
        shutil.copy2(path+'/v_stars_mi.dat', path+'/EQUAL_BELL/v_stars_mi.dat')

    #EQUAL и Белл
    mainf(PATH=path+'/EQUAL_BELL',
        NAME=name,
        INCL=incl,
        SCALE=scale,
        RESOLUTION=resolution,
        H_DISC=h_disc,
        MR=M_R,
        MB=M_B,
        MU0=mu0_c_R,
        R_EFF_B=r_eff_bulge,
        DEG_STAR=pol_degree_star,
        DEG_GAS=pol_degree_gas,
        SIG_MA_DEG=sig_pol_deg,
        SIG_MI_DEG=sig_pol_deg_mi,
        RMIN=Rmin,
        RMAX=Rmax,
        GAS_CORR=gas_corr_by_incl,
        M_TO_L=M_to_L,
        DI=di,
        MONTE_CARLO=monte_carlo_realizations,
        CORRECTION_GAS=correctGasData,
        CORRECTION_STAR=correctStarData,
        CORRECTION_SIG_MA=correctSigmaLosMaj,
        CORRECTION_SIG_MI=correctSigmaLosMin,
        SURF_DENS_STAR=surfaceDensityStarR,
        METHOD='EQUAL',
        PECULIARITIES=peculiarities,
        SIG_WINGS = sig_wings,
        USE_MINOR = use_minor,
        RUN=1)

    renameFilesByMethod(path+'/EQUAL_BELL/', 'EQUAL_BELL')


    if not os.path.exists(path+'/HALF_MAX/'):
        os.makedirs(path+'/HALF_MAX/')
    else:
        for f in os.listdir(path+'/HALF_MAX/'):
            os.remove(path+'/HALF_MAX/'+f)
    shutil.copy2(path+'/v_stars_ma.dat', path+'/HALF_MAX/v_stars_ma.dat')
    shutil.copy2(path+'/v_gas_ma.dat', path+'/HALF_MAX/v_gas_ma.dat')
    shutil.copy2(path+'/gas_density.dat', path+'/HALF_MAX/gas_density.dat')
    if os.path.exists(path+'/v_stars_mi.dat'):
        shutil.copy2(path+'/v_stars_mi.dat', path+'/HALF_MAX/v_stars_mi.dat')

    #HALF и Макс. диск
    mainf(PATH=path+'/HALF_MAX',
        NAME=name,
        INCL=incl,
        SCALE=scale,
        RESOLUTION=resolution,
        H_DISC=h_disc,
        MR=M_R,
        MB=M_B,
        MU0=mu0_c_R,
        R_EFF_B=r_eff_bulge,
        DEG_STAR=pol_degree_star,
        DEG_GAS=pol_degree_gas,
        SIG_MA_DEG=sig_pol_deg,
        SIG_MI_DEG=sig_pol_deg_mi,
        RMIN=Rmin,
        RMAX=Rmax,
        GAS_CORR=gas_corr_by_incl,
        M_TO_L=maxDisc,
        DI=di,
        MONTE_CARLO=monte_carlo_realizations,
        CORRECTION_GAS=correctGasData,
        CORRECTION_STAR=correctStarData,
        CORRECTION_SIG_MA=correctSigmaLosMaj,
        CORRECTION_SIG_MI=correctSigmaLosMin,
        SURF_DENS_STAR=surfaceDensityStarR,
        METHOD='HALF',
        PECULIARITIES=peculiarities,
        SIG_WINGS = sig_wings,
		USE_MINOR = use_minor,
		RUN=2)

    renameFilesByMethod(path+'/HALF_MAX/', 'HALF_MAX')

    if not os.path.exists(path+'/HALF_BELL/'):
        os.makedirs(path+'/HALF_BELL/')
    else:
        for f in os.listdir(path+'/HALF_BELL/'):
            os.remove(path+'/HALF_BELL/'+f)
    shutil.copy2(path+'/v_stars_ma.dat', path+'/HALF_BELL/v_stars_ma.dat')
    shutil.copy2(path+'/v_gas_ma.dat', path+'/HALF_BELL/v_gas_ma.dat')
    shutil.copy2(path+'/gas_density.dat', path+'/HALF_BELL/gas_density.dat')
    if os.path.exists(path+'/v_stars_mi.dat'):
        shutil.copy2(path+'/v_stars_mi.dat', path+'/HALF_BELL/v_stars_mi.dat')

    #HALF и Белл
    mainf(PATH=path+'/HALF_BELL',
        NAME=name,
        INCL=incl,
        SCALE=scale,
        RESOLUTION=resolution,
        H_DISC=h_disc,
        MR=M_R,
        MB=M_B,
        MU0=mu0_c_R,
        R_EFF_B=r_eff_bulge,
        DEG_STAR=pol_degree_star,
        DEG_GAS=pol_degree_gas,
        SIG_MA_DEG=sig_pol_deg,
        SIG_MI_DEG=sig_pol_deg_mi,
        RMIN=Rmin,
        RMAX=Rmax,
        GAS_CORR=gas_corr_by_incl,
        M_TO_L=M_to_L,
        DI=di,
        MONTE_CARLO=monte_carlo_realizations,
        CORRECTION_GAS=correctGasData,
        CORRECTION_STAR=correctStarData,
        CORRECTION_SIG_MA=correctSigmaLosMaj,
        CORRECTION_SIG_MI=correctSigmaLosMin,
        SURF_DENS_STAR=surfaceDensityStarR,
        METHOD='HALF',
        PECULIARITIES=peculiarities,
        SIG_WINGS = sig_wings,
		USE_MINOR = use_minor,
		RUN=3)

    renameFilesByMethod(path+'/HALF_BELL/', 'HALF_BELL')

    if not os.path.exists(path+'/EQUAL_MAX/'):
        os.makedirs(path+'/EQUAL_MAX/')
    else:
        for f in os.listdir(path+'/EQUAL_MAX/'):
            os.remove(path+'/EQUAL_MAX/'+f)
    shutil.copy2(path+'/v_stars_ma.dat', path+'/EQUAL_MAX/v_stars_ma.dat')
    shutil.copy2(path+'/v_gas_ma.dat', path+'/EQUAL_MAX/v_gas_ma.dat')
    shutil.copy2(path+'/gas_density.dat', path+'/EQUAL_MAX/gas_density.dat')
    if os.path.exists(path+'/v_stars_mi.dat'):
        shutil.copy2(path+'/v_stars_mi.dat', path+'/EQUAL_MAX/v_stars_mi.dat')

    #EQUAL и Макс диск
    mainf(PATH=path+'/EQUAL_MAX',
        NAME=name,
        INCL=incl,
        SCALE=scale,
        RESOLUTION=resolution,
        H_DISC=h_disc,
        MR=M_R,
        MB=M_B,
        MU0=mu0_c_R,
        R_EFF_B=r_eff_bulge,
        DEG_STAR=pol_degree_star,
        DEG_GAS=pol_degree_gas,
        SIG_MA_DEG=sig_pol_deg,
        SIG_MI_DEG=sig_pol_deg_mi,
        RMIN=Rmin,
        RMAX=Rmax,
        GAS_CORR=gas_corr_by_incl,
        M_TO_L=maxDisc,
        DI=di,
        MONTE_CARLO=monte_carlo_realizations,
        CORRECTION_GAS=correctGasData,
        CORRECTION_STAR=correctStarData,
        CORRECTION_SIG_MA=correctSigmaLosMaj,
        CORRECTION_SIG_MI=correctSigmaLosMin,
        SURF_DENS_STAR=surfaceDensityStarR,
        METHOD='EQUAL',
        PECULIARITIES=peculiarities,
        SIG_WINGS = sig_wings,
		USE_MINOR = use_minor,
		RUN=4)

    renameFilesByMethod(path+'/EQUAL_MAX/', 'EQUAL_MAX')

#    plt.show()

    finishTime = time.time()
    print '#!!!!!!!!!!!!# Time total: ', (finishTime - startTime), 's'
    print '#!!!!!!!!!!!!# THE END'
