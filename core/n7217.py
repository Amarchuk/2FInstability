__author__ = 'amarch'
# -*- coding: utf-8 -*-

import time
import shutil

from main import *


def correctGasData(r_g1, v_g1, dv_g1):
    '''Функция, куда убраны все операции подгонки с данными по газу.'''

    r_g = r_g1
    v_g = v_g1
    dv_g = dv_g1

    #Если необходимо выпрямить апроксимацию на краю - можно добавить несколько последних точек,
    #это должно помочь сгладить. Или обрезать по upperBord.

    upperBord = 200
    r_g, v_g, dv_g = zip(*(filter(lambda x: x[0] < upperBord, zip(r_g, v_g, dv_g))))
    r_g = list(r_g)
    v_g = list(v_g)
    dv_g = list(dv_g)

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




    r_g = r_g + [i[0] + scale * i[1] for i in zip(r_points * add_points, range(1, add_points + 1))]
    v_g = v_g + v_points * add_points
    dv_g = dv_g + dv_points * add_points

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

    add_points = 50
    r_points = [75]
    v_points = [221]
    dv_points = [5]

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
    r_ma = r_ma[1:-1]
    sig_los_ma = sig_los_ma[1:-1]
    dsig_los_ma = dsig_los_ma[1:-1]

    #    #Если необходимо выпрямить апроксимацию на краю - можно добавить несколько последних точек,
    #    #это должно помочь сгладить.
    #
    multiplate = 10
    addition_points = 1
    r_points = heapq.nlargest(addition_points, r_ma)
    sig_points = []
    dsig_points = []
    for po in r_points:
        sig_points.append(sig_los_ma[r_ma.index(po)])
        dsig_points.append(dsig_los_ma[r_ma.index(po)])
    r_ma = r_ma + [i[0] + scale * i[1] for i in
                   zip(r_points * multiplate, arange(1, 3 * (multiplate * addition_points) + 1, 3))]
    sig_los_ma = sig_los_ma + sig_points * multiplate
    dsig_los_ma = dsig_los_ma + dsig_points * multiplate

    return r_ma, sig_los_ma, dsig_los_ma, x0


def correctSigmaLosMin(r_ma1, sig_los_ma1, dsig_los_ma1):
    '''Корректируем данные по дисперсии скоростей вдоль главной оси. '''

    r_ma, sig_los_ma, dsig_los_ma = map(list, zip(*sorted(zip(r_ma1, sig_los_ma1, dsig_los_ma1))))

    # Можно обрезать в случае плохих краев
    r_ma = r_ma[1:-1]
    sig_los_ma = sig_los_ma[1:-1]
    dsig_los_ma = dsig_los_ma[1:-1]

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
    #    r_ma = r_ma + [i[0] + scale * i[1] for i in zip(r_points * multiplate, arange(1, 5*(multiplate * addition_points) + 1, 5))]
    #    sig_los_ma = sig_los_ma + sig_points * multiplate
    #    dsig_los_ma = dsig_los_ma + dsig_points * multiplate

    return r_ma, sig_los_ma, dsig_los_ma, x0


startTime = time.time()

if __name__ == "__main__":
    plt.rcParams.update({'font.size': 16})

    path = '/home/amarch/Documents/RotationCurves/Diploma/TwoFluidInstAllDataFromSotn17Feb/Sample/RC/U11914_N7217'
    name = 'U11914_N7217'
    incl = 30 # according to Silchenko et al. 2011
    scale = 1
    resolution = 80 #pc/arcsec
    h_disc = 36.8  # R-band
    M_R = 10.38
    M_B = 11.47
    mu0_c_R = 19.91
    r_eff_bulge = 26.2
    pol_degree_star = 25
    pol_degree_gas = 25
    sig_pol_deg = 10
    sig_pol_deg_mi = 12
    Rmin = 20
    Rmax = 115
    M_to_L = 2.2
    #Два диска и данные в полосе I
    mu_1_I = 17.4
    h_1 = 12.5
    mu_2_I = 18.3
    h_2 = 35.8
    gas_corr_by_incl = False
    di = 5
    monte_carlo_realizations = 2
    peculiarities = [26, 35, 70, 80]
    maxDisc = 10.0
    sig_wings = r_eff_bulge # откуда крылья для дисперсий фитировать
    use_minor = False # используется ли дисперсия по малой оси

    #Используем толстый диск в I
    h_disc = h_2


    if not os.path.exists(path+'/EQUAL_IBAND/'):
        os.makedirs(path+'/EQUAL_IBAND/')
    else:
        for f in os.listdir(path+'/EQUAL_IBAND/'):
            os.remove(path+'/EQUAL_IBAND/'+f)
    shutil.copy2(path+'/v_stars_ma.dat', path+'/EQUAL_IBAND/v_stars_ma.dat')
    shutil.copy2(path+'/v_gas_ma.dat', path+'/EQUAL_IBAND/v_gas_ma.dat')
    shutil.copy2(path+'/gas_density.dat', path+'/EQUAL_IBAND/gas_density.dat')
    if os.path.exists(path+'/v_stars_mi.dat'):
        shutil.copy2(path+'/v_stars_mi.dat', path+'/EQUAL_IBAND/v_stars_mi.dat')

    #EQUAL и I-band
    mainf(PATH=path+'/EQUAL_IBAND',
        NAME=name,
        INCL=incl,
        SCALE=scale,
        RESOLUTION=resolution,
        H_DISC=h_disc,
        MR=M_R,
        MB=M_B,
        MU0=mu_2_I,
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
        SURF_DENS_STAR=surfaceDensityStarForTwoDiscs,
        METHOD='EQUAL',
        PECULIARITIES=peculiarities,
        H_DISC_2=h_1,
        MU0_2=mu_1_I,
        SIG_WINGS = sig_wings, 		USE_MINOR = use_minor, 		RUN=1)

    renameFilesByMethod(path+'/EQUAL_IBAND/', 'EQUAL_IBAND')


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
        MU0=mu_2_I,
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
        SURF_DENS_STAR=surfaceDensityStarForTwoDiscs,
        METHOD='HALF',
        PECULIARITIES=peculiarities,
        H_DISC_2=h_1,
        MU0_2=mu_1_I,
        SIG_WINGS = sig_wings, 		USE_MINOR = use_minor, 		RUN=2)

    renameFilesByMethod(path+'/HALF_MAX/', 'HALF_MAX')

    if not os.path.exists(path+'/HALF_IBAND/'):
        os.makedirs(path+'/HALF_IBAND/')
    else:
        for f in os.listdir(path+'/HALF_IBAND/'):
            os.remove(path+'/HALF_IBAND/'+f)
    shutil.copy2(path+'/v_stars_ma.dat', path+'/HALF_IBAND/v_stars_ma.dat')
    shutil.copy2(path+'/v_gas_ma.dat', path+'/HALF_IBAND/v_gas_ma.dat')
    shutil.copy2(path+'/gas_density.dat', path+'/HALF_IBAND/gas_density.dat')
    if os.path.exists(path+'/v_stars_mi.dat'):
        shutil.copy2(path+'/v_stars_mi.dat', path+'/HALF_IBAND/v_stars_mi.dat')

    #HALF и I-band
    mainf(PATH=path+'/HALF_IBAND',
        NAME=name,
        INCL=incl,
        SCALE=scale,
        RESOLUTION=resolution,
        H_DISC=h_disc,
        MR=M_R,
        MB=M_B,
        MU0=mu_2_I,
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
        SURF_DENS_STAR=surfaceDensityStarForTwoDiscs,
        METHOD='HALF',
        PECULIARITIES=peculiarities,
        H_DISC_2=h_1,
        MU0_2=mu_1_I,
        SIG_WINGS = sig_wings, 		USE_MINOR = use_minor, 		RUN=3)

    renameFilesByMethod(path+'/HALF_IBAND/', 'HALF_IBAND')

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
        MU0=mu_2_I,
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
        SURF_DENS_STAR=surfaceDensityStarForTwoDiscs,
        METHOD='EQUAL',
        PECULIARITIES=peculiarities,
        H_DISC_2=h_1,
        MU0_2=mu_1_I,
        SIG_WINGS = sig_wings, 		USE_MINOR = use_minor, 		RUN=4)

    renameFilesByMethod(path+'/EQUAL_MAX/', 'EQUAL_MAX')

    if not os.path.exists(path+'/AD_MAX/'):
        os.makedirs(path+'/AD_MAX/')
    else:
        for f in os.listdir(path+'/AD_MAX/'):
            os.remove(path+'/AD_MAX/'+f)
    shutil.copy2(path+'/v_stars_ma.dat', path+'/AD_MAX/v_stars_ma.dat')
    shutil.copy2(path+'/v_gas_ma.dat', path+'/AD_MAX/v_gas_ma.dat')
    shutil.copy2(path+'/gas_density.dat', path+'/AD_MAX/gas_density.dat')
    if os.path.exists(path+'/v_stars_mi.dat'):
        shutil.copy2(path+'/v_stars_mi.dat', path+'/AD_MAX/v_stars_mi.dat')

    #AD и Макс диск
    mainf(PATH=path+'/AD_MAX',
        NAME=name,
        INCL=incl,
        SCALE=scale,
        RESOLUTION=resolution,
        H_DISC=h_disc,
        MR=M_R,
        MB=M_B,
        MU0=mu_2_I,
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
        SURF_DENS_STAR=surfaceDensityStarForTwoDiscs,
        METHOD='AD',
        PECULIARITIES=peculiarities,
        H_DISC_2=h_1,
        MU0_2=mu_1_I,
        SIG_WINGS = sig_wings, 		USE_MINOR = use_minor, 		RUN=5)

    renameFilesByMethod(path+'/AD_MAX/', 'AD_MAX')

    if not os.path.exists(path+'/AD_IBAND/'):
        os.makedirs(path+'/AD_IBAND/')
    else:
        for f in os.listdir(path+'/AD_IBAND/'):
            os.remove(path+'/AD_IBAND/'+f)
    shutil.copy2(path+'/v_stars_ma.dat', path+'/AD_IBAND/v_stars_ma.dat')
    shutil.copy2(path+'/v_gas_ma.dat', path+'/AD_IBAND/v_gas_ma.dat')
    shutil.copy2(path+'/gas_density.dat', path+'/AD_IBAND/gas_density.dat')
    if os.path.exists(path+'/v_stars_mi.dat'):
        shutil.copy2(path+'/v_stars_mi.dat', path+'/AD_IBAND/v_stars_mi.dat')

    #AD и I-band
    mainf(PATH=path+'/AD_IBAND',
        NAME=name,
        INCL=incl,
        SCALE=scale,
        RESOLUTION=resolution,
        H_DISC=h_disc,
        MR=M_R,
        MB=M_B,
        MU0=mu_2_I,
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
        SURF_DENS_STAR=surfaceDensityStarForTwoDiscs,
        METHOD='AD',
        PECULIARITIES=peculiarities,
        H_DISC_2=h_1,
        MU0_2=mu_1_I,
        SIG_WINGS = sig_wings, 		USE_MINOR = use_minor, 		RUN=6)

    renameFilesByMethod(path+'/AD_IBAND/', 'AD_IBAND')

##
#    #Логгирование в файл
#    sys.stdout = Tee(path + "/log_" + name + '.txt', 'w')
#
#    # Работа с фотометрией в I полосе - два диска.
#    poly_star, poly_gas, star_data, gas_data = bendStarRC(correctGasData, correctStarData, path, incl, 0.0, False,
#        pol_degree_star, pol_degree_gas, name,
#        scale, gas_corr_by_incl, False)
#    h_disc *= scale
#    R1, R2 = correctDistanceInterval(path, scale)
#    R2 = 121
#    evaluateSigLosWingsExpScale(path, r_eff_bulge)
#    sigLosGaussParams, sigMajData = fitGaussSigLosMaj(correctSigmaLosMaj, path, scale, incl)
#    sigLosPolyParams = fitPolySigLosMaj(correctSigmaLosMaj, path, scale, incl, sig_pol_deg, False, min(Rmax, R2))
#    sigLosSinhParams = fitSechSigLosMaj(correctSigmaLosMaj, path, scale, incl)
#    sigLosGaussParamsMi, sigMiData = fitGaussSigLosMin(correctSigmaLosMin, path, scale, incl)
#    sigLosPolyParamsMi = fitPolySigLosMin(correctSigmaLosMin, path, scale, incl, sig_pol_deg_mi, False, min(Rmax, R2))
#    eval_SigPhi_to_sigR(poly_star, R1, R2, (R2 - R1) / 1000.0, path)
#    evalEpyciclicFreq(poly_gas, arange(R1 + 2, R2, 0.1), path, resolution, h_disc)
#    #M_to_L = mass_to_light_Iband(M_B - M_R)
#    print '#!!!!!!!!!!!!# Mass-to-light ratio in I band (M/L) = ', M_to_L
#    plotSurfDensForTwoDiscs(M_to_L, h_1, mu_1_I, h_2, mu_2_I, 0, Rmax, 0.1, path)
#    gas_sf_data = surfaceDensityGas(path)
#
#    r_surfd_gas = gas_sf_data[0]
#    r_star_and_gas = list(arange(Rmin, Rmax, 0.1)) + r_surfd_gas
#    r_star_and_gas.sort()
#    #    r_star_and_gas = filter(lambda x: ((x <= Rmax) & (x >= Rmin)), r_star_and_gas)
#    #    r_surfd_gas = filter(lambda x: ((x <= min(Rmax, R2)) & (x >= max(Rmin, R1, r_eff_bulge))), r_surfd_gas)
#    r_star_and_gas = filter(lambda x: x > r_eff_bulge, r_star_and_gas)
#    r_surfd_gas = filter(lambda x: x > r_eff_bulge, r_surfd_gas)
#
#    h_kin, sigR2 = asymmetricDriftEvaluation(r_star_and_gas, h_disc, path, poly_star, poly_gas, 90)
#    sigZ2, sigPhi2 = velosityEllipsoid(h_disc, r_star_and_gas, sigR2, path, incl, sigLosPolyParams, poly_star)
#
#
#    # Решаем гравнеустойчивость для точек, где есть данные по газовой плотности
#    star_density = [surfaceDensityStarForTwoDiscs(M_to_L, h_1, mu_1_I, h_2, mu_2_I, R) for R in r_surfd_gas]
#    gas_density = [gas_sf_data[1][gas_sf_data[0].index(R)] for R in r_surfd_gas]
#    sigma_corr_gas = [math.sqrt(sigR2Evaluation(R, h_disc, h_kin, poly_star, poly_gas)) for R in r_surfd_gas]
#    Qeffs = findTwoFluidQeffs(r_surfd_gas, poly_gas, gas_density, star_density, sigma_corr_gas, path, resolution, 60.0)
#    hydroQeffs = findTwoFluidHydroQeffs(r_surfd_gas, poly_gas, gas_density, star_density, sigma_corr_gas, path,
#        resolution, 60.0)
#    hzGas = [zGas(R[1], R[2], resolution) / 2 for R in zip(r_surfd_gas, star_density, gas_density)]
#    sigmaZgas = [math.sqrt(sigZ2Evaluation(R, h_disc, h_kin, poly_star, poly_gas, incl, sigLosPolyParams)) for R in
#                 r_surfd_gas]
#    hzStar = [zStar(R[1], R[2], resolution, R[3]) / 2 for R in zip(r_surfd_gas, star_density, gas_density, sigmaZgas)]
#    plotVerticalScale(star_density, gas_density, resolution, sigmaZgas, r_surfd_gas, path)
#    discQeffs = findTwoFluidWithDiscQeffs(r_surfd_gas, poly_gas, gas_density, star_density, sigma_corr_gas, path,
#        resolution, hzStar, hzGas, 60.0)
#    Qeffs1F = findOneFluidQeffs(r_surfd_gas, poly_gas, gas_density, star_density, sigma_corr_gas, path, resolution,
#        60.0)
#
#    # Смотрим, как отразится уменьшение толщины диска в два раза.
#    hzStar = [hzs / 2 for hzs in hzStar]
#    discQeffs_3 = findTwoFluidWithDiscQeffs(r_surfd_gas, poly_gas, gas_density, star_density, sigma_corr_gas, path,
#        resolution, hzStar, hzGas, 60.0)
#    # Смотрим, какие результаты в случае однородно толстого диска 0.2h
#    hzStar = [0.1 * h_disc] * r_surfd_gas.__len__()
#    discQeffs_4 = findTwoFluidWithDiscQeffs(r_surfd_gas, poly_gas, gas_density, star_density, sigma_corr_gas, path,
#        resolution, hzStar, hzGas, 60.0)
#
#
#    # То же для другого угла - чтобы понять зависимость от угла
#    incl = incl + di
#
#    poly_star1, poly_gas1, star_data1, gas_data1 = bendStarRC(correctGasData, correctStarData, path, incl, 0.0, False,
#        pol_degree_star, pol_degree_gas, name,
#        scale, gas_corr_by_incl, False)
#    sigLosPolyParams1 = fitPolySigLosMaj(correctSigmaLosMaj, path, scale, incl, sig_pol_deg, False, min(Rmax, R2))
#    eval_SigPhi_to_sigR(poly_star1, R1, R2, 0.1, path)
#    evalEpyciclicFreq(poly_gas1, arange(R1 + 2, R2, 0.1), path, resolution, h_disc)
#    h_kin_1, sigR2_1 = asymmetricDriftEvaluation(r_star_and_gas, h_disc, path, poly_star1, poly_gas1, 90)
#    sigZ2_1, sigPhi2_1 = velosityEllipsoid(h_disc, r_star_and_gas, sigR2_1, path, incl, sigLosPolyParams1, poly_star1)
#    sigma_corr_gas_1 = [math.sqrt(sigR2Evaluation(R, h_disc, h_kin_1, poly_star1, poly_gas1)) for R in r_surfd_gas]
#    Qeffs_1 = findTwoFluidQeffs(r_surfd_gas, poly_gas1, gas_density, star_density, sigma_corr_gas_1, path, resolution,
#        60.0)
#    hydroQeffs_1 = findTwoFluidHydroQeffs(r_surfd_gas, poly_gas1, gas_density, star_density, sigma_corr_gas_1, path,
#        resolution, 60.0)
#    sigmaZgas = [math.sqrt(sigZ2Evaluation(R, h_disc, h_kin_1, poly_star1, poly_gas1, incl, sigLosPolyParams1)) for R in
#                 r_surfd_gas]
#    hzStar = [zStar(R[1], R[2], resolution, R[3]) / 2 for R in zip(r_surfd_gas, star_density, gas_density, sigmaZgas)]
#    discQeffs_1 = findTwoFluidWithDiscQeffs(r_surfd_gas, poly_gas1, gas_density, star_density, sigma_corr_gas_1, path,
#        resolution, hzStar, hzGas, 60.0)
#    Qeffs1F_1 = findOneFluidQeffs(r_surfd_gas, poly_gas1, gas_density, star_density, sigma_corr_gas_1, path, resolution,
#        60.0)
#
#    # То же для другого угла
#    incl = incl - 2 * di
#
#    poly_star2, poly_gas2, star_data2, gas_data2 = bendStarRC(correctGasData, correctStarData, path, incl, 0.0, False,
#        pol_degree_star, pol_degree_gas, name,
#        scale, gas_corr_by_incl, False)
#    sigLosPolyParams2 = fitPolySigLosMaj(correctSigmaLosMaj, path, scale, incl, sig_pol_deg, False, min(Rmax, R2))
#    eval_SigPhi_to_sigR(poly_star2, R1, R2, 0.1, path)
#    evalEpyciclicFreq(poly_gas2, arange(R1 + 2, R2, 0.1), path, resolution, h_disc)
#    h_kin_2, sigR2_2 = asymmetricDriftEvaluation(r_star_and_gas, h_disc, path, poly_star2, poly_gas2, 90)
#    sigZ2_2, sigPhi2_2 = velosityEllipsoid(h_disc, r_star_and_gas, sigR2_2, path, incl, sigLosPolyParams2, poly_star2)
#    sigma_corr_gas_2 = [math.sqrt(sigR2Evaluation(R, h_disc, h_kin_2, poly_star2, poly_gas2)) for R in r_surfd_gas]
#    Qeffs_2 = findTwoFluidQeffs(r_surfd_gas, poly_gas2, gas_density, star_density, sigma_corr_gas_2, path, resolution,
#        60.0)
#    hydroQeffs_2 = findTwoFluidHydroQeffs(r_surfd_gas, poly_gas2, gas_density, star_density, sigma_corr_gas_2, path,
#        resolution, 60.0)
#    sigmaZgas = [math.sqrt(sigZ2Evaluation(R, h_disc, h_kin_2, poly_star2, poly_gas2, incl, sigLosPolyParams2)) for R in
#                 r_surfd_gas]
#    hzStar = [zStar(R[1], R[2], resolution, R[3]) / 2 for R in zip(r_surfd_gas, star_density, gas_density, sigmaZgas)]
#    discQeffs_2 = findTwoFluidWithDiscQeffs(r_surfd_gas, poly_gas2, gas_density, star_density, sigma_corr_gas_2, path,
#        resolution, hzStar, hzGas, 60.0)
#    Qeffs1F_2 = findOneFluidQeffs(r_surfd_gas, poly_gas2, gas_density, star_density, sigma_corr_gas_2, path, resolution,
#        60.0)
#
#    # Монте-Карло реализации в количестве monte_carlo_realizations штук.
#
#    incl = incl + di
#    sigR2_list = [sigR2]
#    sigZ2_list = [sigZ2]
#    sigPhi2_list = [sigPhi2]
#    Qeffs_list = [zip(*Qeffs)[2]]
#    hydroQeffs_list = [zip(*hydroQeffs)[2]]
#    discQeffs_list = [zip(*discQeffs)[2]]
#    Qeffs1F_list = [Qeffs1F]
#    MC_iter = 1
#
##    while MC_iter < monte_carlo_realizations:
##        MC_iter += 1
##        print '#!!!!!!!!!!!!# Monte-Carlo iterration number ', MC_iter
##        poly_star_mc, poly_gas_mc, star_data_mc, gas_data_mc = bendStarRC(correctGasData, correctStarData, path, incl,
##            0.0, False, pol_degree_star, pol_degree_gas, name,
##            scale, gas_corr_by_incl, True)
##        sigLosPolyParams_mc = fitPolySigLosMaj(correctSigmaLosMaj, path, scale, incl, sig_pol_deg, True, min(Rmax, R2))
##        eval_SigPhi_to_sigR(poly_star_mc, R1, R2, 0.1, path)
##        evalEpyciclicFreq(poly_gas_mc, arange(R1 + 2, R2, 0.1), path, resolution, h_disc)
##        h_kin_mc, sigR2_mc = asymmetricDriftEvaluation(r_star_and_gas, h_disc, path, poly_star_mc, poly_gas_mc, 90)
##        sigZ2_mc, sigPhi2_mc = velosityEllipsoid(h_disc, r_star_and_gas, sigR2, path, incl, sigLosPolyParams_mc,
##            poly_star_mc)
##        sigma_corr_gas_mc = [math.sqrt(sigR2Evaluation(R, h_disc, h_kin_mc, poly_star_mc, poly_gas_mc)) for R in
##                             r_surfd_gas]
##        Qeffs_mc = findTwoFluidQeffs(r_surfd_gas, poly_gas_mc, gas_density, star_density, sigma_corr_gas_mc, path,
##            resolution, 60.0)
##        hydroQeffs_mc = findTwoFluidHydroQeffs(r_surfd_gas, poly_gas_mc, gas_density, star_density, sigma_corr_gas_mc,
##            path,
##            resolution, 60.0)
##        sigmaZgas_mc = [
##        math.sqrt(sigZ2Evaluation(R, h_disc, h_kin_mc, poly_star_mc, poly_gas_mc, incl, sigLosPolyParams_mc)) for R in
##        r_surfd_gas]
##        hzStar_mc = [zStar(R[1], R[2], resolution, R[3]) / 2 for R in
##                     zip(r_surfd_gas, star_density, gas_density, sigmaZgas_mc)]
##        discQeffs_mc = findTwoFluidWithDiscQeffs(r_surfd_gas, poly_gas_mc, gas_density, star_density, sigma_corr_gas_mc,
##            path,
##            resolution, hzStar_mc, hzGas, 60.0)
##        Qeffs1F_mc = findOneFluidQeffs(r_surfd_gas, poly_gas_mc, gas_density, star_density, sigma_corr_gas_mc, path,
##            resolution,
##            60.0)
##        sigR2_list.append(sigR2_mc)
##        sigZ2_list.append(sigZ2_mc)
##        sigPhi2_list.append(sigPhi2_mc)
##        Qeffs_list.append(zip(*Qeffs_mc)[2])
##        hydroQeffs_list.append(zip(*hydroQeffs_mc)[2])
##        discQeffs_list.append(zip(*discQeffs_mc)[2])
##        Qeffs1F_list.append(Qeffs1F_mc)
#
#    plotFinalPics(path, poly_star, poly_gas, di, star_data, gas_data, incl, resolution, h_disc, r_eff_bulge,
#        sigMajData, sigLosGaussParams, sigLosPolyParams, sigLosSinhParams, r_surfd_gas,
#        zip(Qeffs1F, Qeffs1F_1, Qeffs1F_2) + Qeffs1F_list,
#        zip(zip(*hydroQeffs)[2], zip(*hydroQeffs_1)[2], zip(*hydroQeffs_2)[2]) + hydroQeffs_list,
#        zip(zip(*Qeffs)[2], zip(*Qeffs_1)[2], zip(*Qeffs_2)[2]) + Qeffs_list,
#        zip(zip(*discQeffs)[2], zip(*discQeffs_1)[2], zip(*discQeffs_2)[2], zip(*discQeffs_3)[2], zip(*discQeffs_4)[2])
#        + discQeffs_list,
#        r_star_and_gas,
#        zip(sigR2, sigR2_1, sigR2_2) + sigR2_list,
#        zip(sigPhi2, sigPhi2_1, sigPhi2_2) + sigPhi2_list,
#        zip(sigZ2, sigZ2_1, sigZ2_2) + sigZ2_list,
#        hzStar, peculiarities, 1)

#    plt.show()


    finishTime = time.time()
    print '#!!!!!!!!!!!!# Time total: ', (finishTime - startTime), 's'
    print '#!!!!!!!!!!!!# THE END'
