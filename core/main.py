__author__  =  'amarch'
# -*- coding: utf-8 -*-


import sys

from core.rotateAndFitStarRC import *
from instabCriteriaSolution import *
from plotFinal import *


def mainf(**kwargs):
    ''''''
    path  =  kwargs.get('PATH')
    path = kwargs.get('PATH')
    name = kwargs.get('NAME')
    incl = kwargs.get('INCL')
    scale = kwargs.get('SCALE')
    resolution = kwargs.get('RESOLUTION')
    h_disc = kwargs.get('H_DISC')
    M_R = kwargs.get('MR')
    M_B = kwargs.get('MB')
    mu0_c = kwargs.get('MU0')
    r_eff_bulge = kwargs.get('R_EFF_B')
    pol_degree_star = kwargs.get('DEG_STAR')
    pol_degree_gas = kwargs.get('DEG_GAS')
    sig_pol_deg = kwargs.get('SIG_MA_DEG')
    sig_pol_deg_mi = kwargs.get('SIG_MI_DEG')
    Rmin = kwargs.get('RMIN')
    Rmax = kwargs.get('RMAX')
    gas_corr_by_incl = kwargs.get('GAS_CORR')
    M_to_L = kwargs.get('M_TO_L')
    di = kwargs.get('DI')
    monte_carlo_realizations = kwargs.get('MONTE_CARLO')
    correctGasData = kwargs.get('CORRECTION_GAS')
    correctStarData = kwargs.get('CORRECTION_STAR')
    correctSigmaLosMaj = kwargs.get('CORRECTION_SIG_MA')
    correctSigmaLosMin = kwargs.get('CORRECTION_SIG_MI')
    surfaceDensityStar = kwargs.get('SURF_DENS_STAR')
    method = kwargs.get('METHOD')
    peculiarities = kwargs.get('PECULIARITIES')
    run = kwargs.get('RUN')
    sig_wings = kwargs.get('SIG_WINGS')
    use_minor = kwargs.get('USE_MINOR')

    if surfaceDensityStar.__name__ == 'surfaceDensityStarForTwoDiscs':
        h_disc_2 = kwargs.get('H_DISC_2')
        mu0_c_2 = kwargs.get('MU0_2')



    #Логгирование в файл
    sys.stdout  =  Tee(path + "/log_" + name + '.txt', 'w')

    # Работа с фотометрией в I полосе.
    poly_star, poly_gas, star_data, gas_data  =  bendStarRC(correctGasData, correctStarData, path, incl, 0.0, False,
        pol_degree_star, pol_degree_gas, name,
        scale, gas_corr_by_incl, False)
    h_disc *=  scale
    R1, R2  =  correctDistanceInterval(path, scale)
    R1  =  max(R1, r_eff_bulge)
    R2  =  min(R2, Rmax)
    if os.path.exists(path+'/v_stars_mi.dat') and use_minor:
        evaluateAllSigLosWingsExpScale(path, sig_wings)
    else:
        evaluateSigLosWingsExpScale(path, sig_wings)
    sigLosGaussParams, sigMajData  =  fitGaussSigLosMaj(correctSigmaLosMaj, path, scale, incl)
    sigLosPolyParams  =  fitPolySigLosMaj(correctSigmaLosMaj, path, scale, incl, sig_pol_deg, False, min(Rmax, R2))
    sigLosSinhParams  =  fitSechSigLosMaj(correctSigmaLosMaj, path, scale, incl)
    if os.path.exists(path+'/v_stars_mi.dat') and use_minor:
        sigLosGaussParamsMi, sigMiData  =  fitGaussSigLosMin(correctSigmaLosMin, path, scale, incl)
        sigLosPolyParamsMi  =  fitPolySigLosMin(correctSigmaLosMin, path, scale, incl, sig_pol_deg_mi, False, min(Rmax, R2))

    eval_SigPhi_to_sigR(poly_star, R1, R2, (R2 - R1) / 1000.0, path)
    evalEpyciclicFreq(poly_gas, arange(R1 + 2, R2, 0.1), path, resolution, h_disc)
    print '#!!!!!!!!!!!!# Mass-to-light ratio in (M/L)  =  ', M_to_L


    if surfaceDensityStar.__name__ == 'surfaceDensityStarForTwoDiscs':
        plotSurfDensForTwoDiscs(M_to_L, h_disc, mu0_c, h_disc_2, mu0_c_2, 0, Rmax, 0.1, path)
    else:
        plotSurfDens(M_to_L, h_disc, mu0_c, 0, Rmax, 0.1, path, surfaceDensityStar)
    gas_sf_data  =  surfaceDensityGas(path)

    r_surfd_gas = gas_sf_data[0]
    r_star_and_gas = list(arange(Rmin, Rmax, 0.1)) + r_surfd_gas
    r_star_and_gas.sort()
    r_star_and_gas = filter(lambda x: x > r_eff_bulge, r_star_and_gas)
    r_surfd_gas = filter(lambda x: x > r_eff_bulge, r_surfd_gas)

#    ratioSVEfromSigma(r_star_and_gas, h_disc, path, poly_star, sigLosPolyParams, sigLosPolyParamsMi, 100, incl)

    if method == 'HALF':
        sigR2, sigPhi2, sigZ2 = simpleSVEfromSigma(r_star_and_gas, h_disc, path, poly_star, sigMajData,
            sigLosPolyParams, 0.5, Rmax, incl)
    if method == 'EQUAL':
        sigR2, sigPhi2, sigZ2 = simpleSVEwhenPhiEqualsZ(r_star_and_gas, h_disc, path, poly_star, sigMajData,
            sigLosPolyParams, 0.5, Rmax, incl)
    if method == 'AD':
        h_kin, sigR2 = asymmetricDriftEvaluation(r_star_and_gas, h_disc, path, poly_star, poly_gas, Rmax)
        sigZ2, sigPhi2 = velosityEllipsoid(h_disc,r_star_and_gas, sigR2, path, incl, sigLosPolyParams, poly_star)

    if os.path.exists(path+'/v_stars_mi.dat') and use_minor:
        ratioZtoR(r_star_and_gas, r_eff_bulge, path, poly_star, poly_gas, sigLosPolyParams, sigLosPolyParamsMi, 45, incl)

    # Решаем гравнеустойчивость для точек, где есть данные по газовой плотности
    if surfaceDensityStar.__name__ == 'surfaceDensityStarForTwoDiscs':
        star_density = [surfaceDensityStarForTwoDiscs(M_to_L, h_disc, mu0_c, h_disc_2, mu0_c_2, R) for R in r_surfd_gas]
    else:
        star_density = [surfaceDensityStar(M_to_L, h_disc, R, mu0_c) for R in r_surfd_gas]

    gas_density = [gas_sf_data[1][gas_sf_data[0].index(R)] for R in r_surfd_gas]
    sigma_corr_gas = [math.sqrt(sigR2[r_star_and_gas.index(R)]) for R in r_surfd_gas]
    Qeffs = findTwoFluidQeffs(r_surfd_gas, poly_gas, gas_density, star_density, sigma_corr_gas, path, resolution, 60.0)
    hydroQeffs = findTwoFluidHydroQeffs(r_surfd_gas, poly_gas, gas_density, star_density, sigma_corr_gas, path,
        resolution, 60.0)
    hzGas = [zGas(R[1], R[2], resolution) / 2 for R in zip(r_surfd_gas, star_density, gas_density)]
    sigmaZgas = [math.sqrt(sigZ2[r_star_and_gas.index(R)]) for R in r_surfd_gas]
    hzStar = [zStar(R[1], R[2], resolution, R[3]) / 2 for R in zip(r_surfd_gas, star_density, gas_density, sigmaZgas)]
    plotVerticalScale(star_density, gas_density, resolution, sigmaZgas, r_surfd_gas, path)
    discQeffs = findTwoFluidWithDiscQeffs(r_surfd_gas, poly_gas, gas_density, star_density, sigma_corr_gas, path,
        resolution, hzStar, hzGas, 60.0)
    Qeffs1F = findOneFluidQeffs(r_surfd_gas, poly_gas, gas_density, star_density, sigma_corr_gas, path, resolution,
        60.0)
    set5533Null()

    # Смотрим, как отразится уменьшение толщины диска в два раза.
    hzStar = [hzs / 2 for hzs in hzStar]
    discQeffs_3 = findTwoFluidWithDiscQeffs(r_surfd_gas, poly_gas, gas_density, star_density, sigma_corr_gas, path,
        resolution, hzStar, hzGas, 60.0)
    # Смотрим, какие результаты в случае однородно толстого диска 0.1h
    hzStar = [0.1 * h_disc] * r_surfd_gas.__len__()
    discQeffs_4 = findTwoFluidWithDiscQeffs(r_surfd_gas, poly_gas, gas_density, star_density, sigma_corr_gas, path,
        resolution, hzStar, hzGas, 60.0)

    # То же для другого угла - чтобы понять зависимость от угла
    r_ma = zip(*star_data)[0]
    poly_star1 = bendPolinom(poly_star, min(r_ma), max(r_ma), incl, incl+di)
    r_ma = zip(*gas_data)[0]
    poly_gas1 = bendPolinom(poly_gas, min(r_ma), max(r_ma), incl, incl+di)
    eval_SigPhi_to_sigR(poly_star1, R1, R2, 0.1, path)
    evalEpyciclicFreq(poly_gas1, arange(R1 + 2, R2, 0.1), path, resolution, h_disc)
    if method == 'HALF':
        sigR2_1, sigPhi2_1, sigZ2_1 = simpleSVEfromSigma(r_star_and_gas, h_disc, path, poly_star, sigMajData,
            sigLosPolyParams, 0.5, Rmax, incl)
    if method == 'EQUAL':
        sigR2_1, sigPhi2_1, sigZ2_1 = simpleSVEwhenPhiEqualsZ(r_star_and_gas, h_disc, path, poly_star, sigMajData,
            sigLosPolyParams, 0.5, Rmax, incl)
    if method == 'AD':
        h_kin_1, sigR2_1 = asymmetricDriftEvaluation(r_star_and_gas, h_disc, path, poly_star, poly_gas, Rmax)
        sigZ2_1, sigPhi2_1 = velosityEllipsoid(h_disc,r_star_and_gas, sigR2, path, incl, sigLosPolyParams, poly_star)
    sigma_corr_gas_1 = [math.sqrt(sigR2_1[r_star_and_gas.index(R)]) for R in r_surfd_gas]
    Qeffs_1 = findTwoFluidQeffs(r_surfd_gas, poly_gas1, gas_density, star_density, sigma_corr_gas_1, path, resolution,
        60.0)
    hydroQeffs_1 = findTwoFluidHydroQeffs(r_surfd_gas, poly_gas1, gas_density, star_density, sigma_corr_gas_1, path,
        resolution, 60.0)
    sigmaZgas = [math.sqrt(sigZ2_1[r_star_and_gas.index(R)]) for R in r_surfd_gas]
    hzStar = [zStar(R[1], R[2], resolution, R[3]) / 2 for R in zip(r_surfd_gas, star_density, gas_density, sigmaZgas)]
    discQeffs_1 = findTwoFluidWithDiscQeffs(r_surfd_gas, poly_gas1, gas_density, star_density, sigma_corr_gas_1, path,
        resolution, hzStar, hzGas, 60.0)
    Qeffs1F_1 = findOneFluidQeffs(r_surfd_gas, poly_gas1, gas_density, star_density, sigma_corr_gas_1, path, resolution,
        60.0)
    set5533Null()

    r_ma = zip(*star_data)[0]
    poly_star2 = bendPolinom(poly_star, min(r_ma), max(r_ma), incl, incl-di)
    r_ma = zip(*gas_data)[0]
    poly_gas2 = bendPolinom(poly_gas, min(r_ma), max(r_ma), incl, incl-di)
    eval_SigPhi_to_sigR(poly_star2, R1, R2, 0.1, path)
    evalEpyciclicFreq(poly_gas2, arange(R1 + 2, R2, 0.1), path, resolution, h_disc)
    if method == 'HALF':
        sigR2_2, sigPhi2_2, sigZ2_2 = simpleSVEfromSigma(r_star_and_gas, h_disc, path, poly_star, sigMajData,
            sigLosPolyParams, 0.5, Rmax, incl)
    if method == 'EQUAL':
        sigR2_2, sigPhi2_2, sigZ2_2 = simpleSVEwhenPhiEqualsZ(r_star_and_gas, h_disc, path, poly_star, sigMajData,
            sigLosPolyParams, 0.5, Rmax, incl)
    if method == 'AD':
        h_kin_2, sigR2_2 = asymmetricDriftEvaluation(r_star_and_gas, h_disc, path, poly_star, poly_gas, Rmax)
        sigZ2_2, sigPhi2_2 = velosityEllipsoid(h_disc,r_star_and_gas, sigR2, path, incl, sigLosPolyParams, poly_star)
    sigma_corr_gas_2 = [math.sqrt(sigR2_2[r_star_and_gas.index(R)]) for R in r_surfd_gas]
    Qeffs_2 = findTwoFluidQeffs(r_surfd_gas, poly_gas2, gas_density, star_density, sigma_corr_gas_2, path, resolution,
        60.0)
    hydroQeffs_2 = findTwoFluidHydroQeffs(r_surfd_gas, poly_gas2, gas_density, star_density, sigma_corr_gas_2, path,
        resolution, 60.0)
    sigmaZgas = [math.sqrt(sigZ2_2[r_star_and_gas.index(R)]) for R in r_surfd_gas]
    hzStar = [zStar(R[1], R[2], resolution, R[3]) / 2 for R in zip(r_surfd_gas, star_density, gas_density, sigmaZgas)]
    discQeffs_2 = findTwoFluidWithDiscQeffs(r_surfd_gas, poly_gas2, gas_density, star_density, sigma_corr_gas_2, path,
        resolution, hzStar, hzGas, 60.0)
    Qeffs1F_2 = findOneFluidQeffs(r_surfd_gas, poly_gas2, gas_density, star_density, sigma_corr_gas_2, path, resolution,
        60.0)
    set5533Null()

    # Монте-Карло реализации в количестве monte_carlo_realizations штук.

    sigR2_list = [sigR2]
    sigZ2_list = [sigZ2]
    sigPhi2_list = [sigPhi2]
    Qeffs_list = [zip(*Qeffs)[2]]
    hydroQeffs_list = [zip(*hydroQeffs)[2]]
    discQeffs_list = [zip(*discQeffs)[2]]
    Qeffs1F_list = [Qeffs1F]
    MC_iter = 1

#    while MC_iter < monte_carlo_realizations:
#        MC_iter += 1
#        print '#!!!!!!!!!!!!# Monte-Carlo iterration number ', MC_iter
#        poly_star_mc, poly_gas_mc, star_data_mc, gas_data_mc = bendStarRC(correctGasData, correctStarData, path, incl,
#            0.0, False, pol_degree_star, pol_degree_gas, name, scale, gas_corr_by_incl, True)
#        sigLosPolyParams_mc = fitPolySigLosMaj(correctSigmaLosMaj, path, scale, incl, sig_pol_deg, True, min(Rmax, R2))
#        eval_SigPhi_to_sigR(poly_star_mc, R1, R2, 0.1, path)
#        evalEpyciclicFreq(poly_gas_mc, arange(R1 + 2, R2, 0.1), path, resolution, h_disc)
#        sigR2_mc, sigPhi2_mc, sigZ2_mc = SVEfunction(r_star_and_gas, h_disc, path, poly_star, sigMajData,
#            sigLosPolyParams, 0.5, 71, incl)
#        sigma_corr_gas_mc = [math.sqrt(sigR2_mc[r_star_and_gas.index(R)]) for R in r_surfd_gas]
#        Qeffs_mc = findTwoFluidQeffs(r_surfd_gas, poly_gas_mc, gas_density, star_density, sigma_corr_gas_mc, path,
#            resolution, 60.0)
#        hydroQeffs_mc = findTwoFluidHydroQeffs(r_surfd_gas, poly_gas_mc, gas_density, star_density, sigma_corr_gas_mc,
#            path,
#            resolution, 60.0)
#        sigmaZgas_mc = [math.sqrt(sigZ2_mc[r_star_and_gas.index(R)]) for R in r_surfd_gas]
#        hzStar_mc = [zStar(R[1], R[2], resolution, R[3]) / 2 for R in
#                     zip(r_surfd_gas, star_density, gas_density, sigmaZgas_mc)]
#        discQeffs_mc = findTwoFluidWithDiscQeffs(r_surfd_gas, poly_gas_mc, gas_density, star_density, sigma_corr_gas_mc,
#            path,
#            resolution, hzStar_mc, hzGas, 60.0)
#        Qeffs1F_mc = findOneFluidQeffs(r_surfd_gas, poly_gas_mc, gas_density, star_density, sigma_corr_gas_mc, path,
#            resolution,
#            60.0)
#        sigR2_list.append(sigR2_mc)
#        sigZ2_list.append(sigZ2_mc)
#        sigPhi2_list.append(sigPhi2_mc)
#        Qeffs_list.append(zip(*Qeffs_mc)[2])
#        hydroQeffs_list.append(zip(*hydroQeffs_mc)[2])
#        discQeffs_list.append(zip(*discQeffs_mc)[2])
#        Qeffs1F_list.append(Qeffs1F_mc)
#
    plotFinalPics(path, poly_star, poly_gas, di, star_data, gas_data, incl, resolution, h_disc, r_eff_bulge,
        sigMajData, sigLosGaussParams, sigLosPolyParams, sigLosSinhParams, r_surfd_gas,
        zip(Qeffs1F, Qeffs1F_1, Qeffs1F_2) + Qeffs1F_list,
        zip(zip(*hydroQeffs)[2], zip(*hydroQeffs_1)[2], zip(*hydroQeffs_2)[2]) + hydroQeffs_list,
        zip(zip(*Qeffs)[2], zip(*Qeffs_1)[2], zip(*Qeffs_2)[2]) + Qeffs_list,
        zip(zip(*discQeffs)[2], zip(*discQeffs_1)[2], zip(*discQeffs_2)[2], zip(*discQeffs_3)[2], zip(*discQeffs_4)[2])
        + discQeffs_list,
        r_star_and_gas,
        zip(sigR2, sigR2_1, sigR2_2) + sigR2_list,
        zip(sigPhi2, sigPhi2_1, sigPhi2_2) + sigPhi2_list,
        zip(sigZ2, sigZ2_1, sigZ2_2) + sigZ2_list,
        hzStar, peculiarities, (run-1)*5)