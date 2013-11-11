__author__ = 'amarch'
# -*- coding: utf-8 -*-

from numpy import *
import matplotlib.pyplot as plt
import scipy
import math
from velocityEllipsoidReconstr import *



# Отрисовка чистовых картинок. Разбито на три блока:
# 1) наблюдательные данные - кривые вращения с приближенияит и дисперсии с приближениями
# 2) извлекаемые профили - эллипсоид скоростей SVE и звездный профиль толщины
# 3) критерий неустойчивости - одножидкостный, гидродинамическое и кинетическое приближения двухж. и с учетом толшины

def plotFinalPics(path, poly_star, poly_gas, delta_i, star_data, gas_data, incl, resolution, h_disc, r_eff_bulge,
                  sigLosMajData, sigLosGaussParams, sigLosPolyParams, sigLosSinhParams, r_gas, Q1F, Q2FHydro, Q2FKinem,
                  Q2F_withDisc, r_gas_and_star, sigR2, sigPhi2, sigZ2, hzStar, peculiarities, num):
    '''Главная функция, откуда строятся уже все остальные картинки.'''
    plotFinalRC(path, poly_star, poly_gas, delta_i, star_data, gas_data, incl, resolution, h_disc, r_eff_bulge, num)
    plotFinalSigLosMaj(path, sigLosMajData, sigLosGaussParams, sigLosPolyParams, sigLosSinhParams, num)
    plotFinalInstabCriteria(path, r_gas, Q1F, Q2FHydro, Q2FKinem, Q2F_withDisc, peculiarities, num, sigLosMajData)
    plotFinalSVE(path, incl, r_gas_and_star, sigLosMajData, sigR2, sigPhi2, sigZ2, delta_i, num, r_gas)
    plotFinalVertScale(path, r_gas, hzStar, num)


def plotFinalRC(path, poly_star, poly_gas, delta_i, star_data, gas_data, incl, resolution, h_disc, r_eff_bulge, num):
    '''Профили скоростей для звезд и газа - наблюдения и приближения + отклонение по углу наклона на delta_i.
    Также сохранение всех данных в файл для анализа и чтобы не пересчитывать.'''

    def resize(x):
        '''Перевод угловых секунд в килопарсеки.'''
        return x * resolution / 1000

    outRC = open(path + '/finalRC.dat', 'w+')

    r_ma, v_ma, dv_ma = map(list, zip(*star_data))
    r_g, v_g, dv_g = map(list, zip(*gas_data))

    fig = plt.figure(20+num)
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax2.errorbar(r_ma, v_ma, yerr=dv_ma, fmt='.', marker='.', mew=0, label=r'$\bar{v}_{\varphi}(R)$', color='red')
    ax2.errorbar(r_g, v_g, yerr=dv_g, fmt='.', marker='.', mew=0, label='$v_c(R)$')
    xx = arange(min(r_ma), max(r_ma), 0.1)
    ax2.plot(xx, poly_star(xx), color='black', label='$approximations$')
    yy = arange(min(r_g), max(r_g), 0.1)
    ax2.plot(yy, poly_gas(yy), color='black')
    corr_by_incl = math.sin((incl + delta_i) * math.pi / 180) / math.sin(incl * math.pi / 180)
    ax2.plot(xx, corr_by_incl * poly_star(xx), '--', color='black')
    ax2.plot(yy, corr_by_incl * poly_gas(yy), '--', color='black')
    corr_by_incl = math.sin((incl - delta_i) * math.pi / 180) / math.sin(incl * math.pi / 180)
    ax2.plot(xx, corr_by_incl * poly_star(xx), '--', color='black')
    ax2.plot(yy, corr_by_incl * poly_gas(yy), '--', color='black')

    ax2.set_xlim(0)
    x1, x2 = ax2.get_xlim()
    ax1.set_xlim(resize(x1), resize(x2))

    ax2.set_ylim(0, max(v_g) + 200)
    ax2.axvline(x=h_disc, ymin=0, ymax=0.1, ls='-', color='black')
    ax2.axvline(x=r_eff_bulge, ymin=0, ymax=0.1, ls='--', color='black')
    ax1.set_ylabel("$v(R),\ km/s$")
    ax2.set_xlabel("$R,\ arcsec$")
    ax1.set_xlabel("$R,\ kpc$")
    ax2.legend(loc='upper right', numpoints=1).draw_frame(False)
    ax2.text(0.1, 0.9,'(a)', ha='center', va='center', transform= ax2.transAxes)

    outRC.write('#Star points (r, v, dv)\n')
    for r in zip(r_ma, v_ma, dv_ma):
        outRC.write(str(r[0]) + '  ' + str(r[1]) + '  ' + str(r[2]) + '\n')

    outRC.write('\n#Star and gas polinom (r, p_star, p_gas)\n')
    for r in xx:
        outRC.write(str(r) + '  ' + str(poly_star(r)) + '  ' + str(poly_gas(r)) + '\n')

    outRC.write('\n#Gas points (r, v, dv)\n')
    for r in zip(r_g, v_g, dv_g):
        outRC.write(str(r[0]) + '  ' + str(r[1]) + '  ' + str(r[2]) + '\n')

    fig.savefig(path + "/finalRC.png")
    fig.savefig(path + "/finalRC.eps", format='eps')
    outRC.close()


def plotFinalSigLosMaj(path, sigLosMajData, sigLosGaussParams, sigLosPolyParams, sigLosSinhParams, num):
    '''Дисперсия скоростей вдоль главной оси - наблюдения + все три приближения'''

    outSig = open(path + '/finalSimaLosMaj.dat', 'w+')

    r_ma, sig_los_ma, dsig_los_ma = map(list, zip(*sigLosMajData))
    f = plt.figure(21+num)
    ax = f.add_subplot(111)
#    ax.plot(r_ma, sig_los_ma, '.', color='b', label='$\mathrm{maj}$')
    ax.plot(r_ma, sig_los_ma, '.', color='b')
#    plt.plot(r_ma, fsech(r_ma, sigLosSinhParams), '-', label='$sech^2$')
#    plt.plot(r_ma, gauss(r_ma, sigLosGaussParams), '-', label='$gaussian$')
    ax.plot(list(arange(-0.1, -max(r_ma), -0.1)[::-1]) + list(arange(0, max(r_ma), 0.1)),
        list(sigLosPolyParams(arange(0.1, max(r_ma), 0.1))[::-1]) + list(sigLosPolyParams(arange(0, max(r_ma), 0.1))),
        '-')
    ax.errorbar(r_ma, sig_los_ma, yerr=dsig_los_ma, fmt=None, marker=None, mew=0, ecolor='b')
    ax.set_ylabel("$\sigma_{\mathrm{los}},\ km/s$")
    ax.set_xlabel("$R,\ arcsec$")
    ax.text(0.05, 0.9,'(b)', ha='center', va='center', transform= ax.transAxes)
#    ax.legend(loc='upper right', numpoints=1).draw_frame(False)

    outSig.write('#Sigma LOS Maj (r, sig, dsig, gauss, sech, poly)\n')
    for r in zip(r_ma, sig_los_ma, dsig_los_ma, fsech(r_ma, sigLosSinhParams), gauss(r_ma, sigLosGaussParams),
        list(sigLosPolyParams(arange(0.1, max(r_ma), 0.1))[::-1]) + list(sigLosPolyParams(arange(0, max(r_ma), 0.1)))):
        outSig.write(str(r[0]) + '  ' + str(r[1]) + '  ' + str(r[2]) + '  ' + str(r[3]) + '  ' + str(r[4]) + '  ' + str(
            r[5]) + '\n')

    f.savefig(path + "/finalSigLosMaj.png")
    f.savefig(path + "/finalSigLosMaj.eps", format='eps')
    outSig.close()


def plotFinalInstabCriteria(path, r_gas, Q1F, Q2FHydro, Q2FKinem, Q2F_withDisc, peculiarities, num, sigLosMajData):
    '''Главная картинка: строим отношение поверхностной плотности газа к критической для разных приближений.
    Везде учитываются нужные поправки - 1.6 или 0.67. Также вариация по углу наклона.'''

    outSigma = open(path + '/finalSigma.dat', 'w+')

    length = r_gas.__len__()
    Q1F_0, Q1F_1, Q1F_2 = zip(*(Q1F[0:length]))
    Q1F_list = Q1F[length:]
    Q2FHydro_0, Q2FHydro_1, Q2FHydro_2 = zip(*(Q2FHydro[0:length]))
    Q2FHydro_list = Q2FHydro[length:]
    Q2FKinem_0, Q2FKinem_1, Q2FKinem_2 = zip(*(Q2FKinem[0:length]))
    Q2FKinem_list = Q2FKinem[length:]
    Q2F_withDisc_0, Q2F_withDisc_1, Q2F_withDisc_2, Q2F_withDisc_3, Q2F_withDisc_4 = zip(*(Q2F_withDisc[0:length]))
    Q2F_withDisc_list = Q2F_withDisc[length:]

    Q2FKinem_med = []
    Q2FHydro_med = []
    Q2F_withDisc_med = []
    Q1F_med = []

    Q2FKinem_std = []
    Q2FHydro_std = []
    Q2F_withDisc_std = []
    Q1F_std = []

    for q in zip(*Q2FKinem_list):
        Q2FKinem_med.append(scipy.mean(map(lambda x: 1.6 / x, q)))
        Q2FKinem_std.append(scipy.std(map(lambda x: 1.6 / x, q)))

    for q in zip(*Q2FHydro_list):
        Q2FHydro_med.append(scipy.mean(map(lambda x: 1.6 / x, q)))
        Q2FHydro_std.append(scipy.std(map(lambda x: 1.6 / x, q)))

    for q in zip(*Q2F_withDisc_list):
        Q2F_withDisc_med.append(scipy.mean(map(lambda x: 1.6 / x, q)))
        Q2F_withDisc_std.append(scipy.std(map(lambda x: 1.6 / x, q)))

    for q in zip(*Q1F_list):
        Q1F_med.append(scipy.mean(map(lambda x: 1 / (0.67 * x), q)))
        Q1F_std.append(scipy.std(map(lambda x: 1 / (0.67 * x), q)))

    f = plt.figure(22+num)
    ax = f.add_subplot(111)

#    for pec in peculiarities:
#        ax.axvline(x=pec, ymin=0, ymax=0.05, ls='--', color = 'black')

    ax.axhline(y=1, ls='--', color="black")
    r_ma, sig_los_ma, dsig_los_ma = map(list, zip(*sigLosMajData))
    ax.axvline(x=max(r_ma), ls='-.', color="black")

    ax.plot(r_gas, map(lambda x: 1 / (0.67 * x), Q1F_1), '--', color='blue')
    ax.plot(r_gas, map(lambda x: 1.6 / x, Q2FHydro_1), '--', color='red')
    ax.plot(r_gas, map(lambda x: 1.6 / x, Q2FKinem_1), '--', color='green')
    ax.plot(r_gas, map(lambda x: 1.6 / x, Q2F_withDisc_1), '--', color='black')

    ax.plot(r_gas, map(lambda x: 1 / (0.67 * x), Q1F_2), '--', color='blue')
    ax.plot(r_gas, map(lambda x: 1.6 / x, Q2FHydro_2), '--', color='red')
    ax.plot(r_gas, map(lambda x: 1.6 / x, Q2FKinem_2), '--', color='green')
    ax.plot(r_gas, map(lambda x: 1.6 / x, Q2F_withDisc_2), '--', color='black')

    ax.plot(r_gas, map(lambda x: 1 / (0.67 * x), Q1F_0), 's-', label='$1Fluid$', color='blue')
    ax.plot(r_gas, map(lambda x: 1.6 / x, Q2FHydro_0), 's-', label='$2F\ Hydrodynamics$', color='red')
    ax.plot(r_gas, map(lambda x: 1.6 / x, Q2FKinem_0), 's-', label='$2F\ Kinetics$', color='green')
    ax.plot(r_gas, map(lambda x: 1.6 / x, Q2F_withDisc_0), 's-', label='$2F\ with\ h_{z}$', color='black')
#    ax.plot(r_gas, map(lambda x: 1.6 / x, Q2F_withDisc_3), 'x-', label='$2F\ thickness\ h_z^{star}/2$', color='grey')
    ax.plot(r_gas, map(lambda x: 1.6 / x, Q2F_withDisc_4), 's-', label='$2F\ h_{z}^{\mathrm{s}}=0.1h_{\mathrm{d}}$', color='grey')



    for i in range(0,peculiarities.__len__(),2):
        ax.plot([peculiarities[i],peculiarities[i+1]], [0.0, 0.0], linewidth=15, linestyle='-',c="red")

#    ax.plot(r_gas, Q1F_med, 's', color='blue')
#    ax.plot(r_gas, Q2FHydro_med, 's', color='red')
#    ax.plot(r_gas, Q2FKinem_med, 's', color='green')
#    ax.plot(r_gas, Q2F_withDisc_med, 's', color='black')

#    ax.errorbar(r_gas, Q1F_med, yerr=Q1F_std, fmt=None, marker=None, mew=0, ecolor='blue')
#    ax.errorbar(r_gas, Q2FHydro_med, yerr=Q2FHydro_std, fmt=None, marker=None, mew=0, ecolor='red')
#    ax.errorbar(r_gas, Q2FKinem_med, yerr=Q2FKinem_std, fmt=None, marker=None, mew=0, ecolor='green')
#    ax.errorbar(r_gas, Q2F_withDisc_med, yerr=Q2F_withDisc_std, fmt=None, marker=None, mew=0, ecolor='black')

    outSigma.write(
        '#SigmaGas/SigmaCR (1)r_gas  (2)Q1F (3)Q1F_med (4)Q1F_std (5)Q2FHydro (6)Q2FHydro_med (7)Q2FHydro_std (8)Q2FKinet (9)Q2FKinet_med (10)Q2FKinet_std (11)Q2F_withDisc (12)Q2F_withDisc_med (13)Q2F_withDisc_std)\n')
    for r in zip(r_gas, map(lambda x: 1 / (0.67 * x), Q1F_0), Q1F_med, Q1F_std, map(lambda x: 1.6 / x, Q2FHydro_0),
        Q2FHydro_med, Q2FHydro_std, map(lambda x: 1.6 / x, Q2FKinem_0), Q2FKinem_med, Q2FKinem_std,
        map(lambda x: 1.6 / x, Q2F_withDisc_0), Q2F_withDisc_med, Q2F_withDisc_std):
        outSigma.write(
            str(r[0]) + '  ' + str(r[1]) + '  ' + str(r[2]) + '  ' + str(r[3]) + '  ' + str(r[4]) + '  ' + str(
                r[5]) + '  ' + str(r[6]) + '  ' + str(r[7]) + '  ' + str(r[8]) + '  ' + str(r[9]) + '  ' + str(
                r[10]) + '  ' + str(r[11]) + '  ' + str(r[12]) + '\n')

    ax.set_ylabel("$\Sigma_{\mathrm{g}}/\Sigma_{\mathrm{g}}^{\mathrm{cr}}$")
    ax.set_xlabel("$R,\ arcsec$")
    ax.set_xlim(0, max(r_gas) + 80)
    ax.set_ylim(0, 2.3)
#    ax.legend(loc='upper right', numpoints=1).draw_frame(False)
    ax.legend(loc='upper right', numpoints=1).draw_frame(True)
    ax.text(0.9, 0.9,'(d)', ha='center', va='center', transform= ax.transAxes)

    f.savefig(path + "/finalQeffs.png")
    f.savefig(path + "/finalQeffs.eps", format='eps')
    outSigma.close()


def plotFinalSVE(path, incl, r_gas_and_star, sigLosMajData, sigR2, sigPhi2, sigZ2, delta_i, num, r_gas):
    '''Отрисовка полученного ранее SVE для соответсвующего r_gas_and_star. Пунктиром показана зависимость от наклона.
    Также нарисовано теоретически восстановленное значение sigLosMaj и сравнено с точками. Также вариация по углу.'''
    r_ma, sig_los_ma, dsig_los_ma = zip(*sigLosMajData)
    incl_rad = incl * math.pi / 180

    outSVE = open(path + '/finalSVE.dat', 'w+')
    length = r_gas_and_star.__len__()

    sigR2_0, sigR2_1, sigR2_2 = zip(*(sigR2[0:length]))
    sigR2_list = sigR2[length:]
    sigZ2_0, sigZ2_1, sigZ2_2 = zip(*(sigZ2[0:length]))
    sigZ2_list = sigZ2[length:]
    sigPhi2_0, sigPhi2_1, sigPhi2_2 = zip(*(sigPhi2[0:length]))
    sigPhi2_list = sigPhi2[length:]

    sigR2_med = []
    sigZ2_med = []
    sigPhi2_med = []

    sigR2_std = []
    sigZ2_std = []
    sigPhi2_std = []

    for s in zip(*sigR2_list)[1::35]:
        sigR2_med.append(scipy.mean(map(math.sqrt, s)))
        sigR2_std.append(scipy.std(map(math.sqrt, s)))

    for s in zip(*sigZ2_list)[5::35]:
        sigZ2_med.append(scipy.mean(map(math.sqrt, s)))
        sigZ2_std.append(scipy.std(map(math.sqrt, s)))

    for s in zip(*sigPhi2_list)[10::35]:
        sigPhi2_med.append(scipy.mean(map(math.sqrt, s)))
        sigPhi2_std.append(scipy.std(map(math.sqrt, s)))

    f = plt.figure(23+num)
    ax = f.add_subplot(111)



    ax.plot(r_gas_and_star, map(math.sqrt, sigR2_1), '--', color='blue')
    ax.plot(r_gas_and_star, map(math.sqrt, sigZ2_1), '--', color='red')
    ax.plot(r_gas_and_star, map(math.sqrt, sigPhi2_1), '--', color='green')

    ax.plot(r_gas_and_star, map(math.sqrt, sigR2_2), '--', color='blue')
    ax.plot(r_gas_and_star, map(math.sqrt, sigZ2_2), '--', color='red')
    ax.plot(r_gas_and_star, map(math.sqrt, sigPhi2_2), '--', color='green')

    ax.plot(r_gas_and_star, map(math.sqrt, sigR2_0), '-', label='$\sigma_{R}$', color='blue')
    ax.plot(r_gas_and_star, map(math.sqrt, sigZ2_0), '-', label='$\sigma_{z}$', color='red')
    ax.plot(r_gas_and_star, map(math.sqrt, sigPhi2_0), '-', label=r'$\sigma_{\varphi}$', color='green')

#    ax.plot(r_gas_and_star[1::35], sigR2_med, 's', color='blue')
#    ax.plot(r_gas_and_star[5::35], sigZ2_med, 's', color='red')
#    ax.plot(r_gas_and_star[10::35], sigPhi2_med, 's', color='green')
#
#    ax.errorbar(r_gas_and_star[1::35], sigR2_med, yerr=sigR2_std, fmt=None, marker=None, mew=0, ecolor='blue')
#    ax.errorbar(r_gas_and_star[5::35], sigZ2_med, yerr=sigZ2_std, fmt=None, marker=None, mew=0, ecolor='red')
#    ax.errorbar(r_gas_and_star[10::35], sigPhi2_med, yerr=sigPhi2_std, fmt=None, marker=None, mew=0, ecolor='green')

    ax.plot(map(abs, r_ma), sig_los_ma, 'o', label=r'$\sigma_{\mathrm{los}}^{\mathrm{maj}}$',color='cyan')
    ax.errorbar(map(abs, r_ma), sig_los_ma, yerr=dsig_los_ma, fmt=None, marker=None, mew=0, color='cyan')

    sig_maj_reestablish = [math.sqrt(x[0] * (math.sin(incl_rad) ** 2) + x[1] * (math.cos(incl_rad) ** 2)) for x in
                           zip(sigPhi2_0, sigZ2_0)]
    ax.plot(r_gas_and_star, sig_maj_reestablish, '--', color='black')

    outSVE.write('#SVE points (r, sigmaR, sigmaZ, sigmaPhi, sigmaMajRestored)\n')
    for r in zip(r_gas_and_star, map(math.sqrt, sigR2_1), map(math.sqrt, sigZ2_1), map(math.sqrt, sigPhi2_1),
        sig_maj_reestablish):
        outSVE.write(str(r[0]) + '  ' + str(r[1]) + '  ' + str(r[2]) + '  ' + str(r[3]) + '  ' + str(r[4]) + '\n')

    outSVE.write('\n#SVE monte-carlo points sigmaR (r, sigmaR_med, sigmaR_std)\n')
    for r in zip(r_gas_and_star[1::35], sigR2_med, sigR2_std):
        outSVE.write(str(r[0]) + '  ' + str(r[1]) + '  ' + str(r[2]) + '\n')

    outSVE.write('\n#SVE monte-carlo points sigmaZ (r, sigmaZ_med, sigmaZ_std)\n')
    for r in zip(r_gas_and_star[5::35], sigZ2_med, sigZ2_std):
        outSVE.write(str(r[0]) + '  ' + str(r[1]) + '  ' + str(r[2]) + '\n')

    outSVE.write('\n#SVE monte-carlo points sigmaPhi (r, sigmaPhi_med, sigmaPhi_std)\n')
    for r in zip(r_gas_and_star[10::35], sigPhi2_med, sigPhi2_std):
        outSVE.write(str(r[0]) + '  ' + str(r[1]) + '  ' + str(r[2]) + '\n')

    ax.set_xlim(0)
    ax.set_ylabel("$\sigma,\ km/s$")
    ax.set_xlabel("$R,\ arcsec$")
    ax.legend(loc='upper right', numpoints=1).draw_frame(False)
    mamamax = max(sigR2_0+sigR2_1+sigR2_2+sigZ2_0+sigZ2_1+sigZ2_2)
    mimin = min(sigR2_0+sigR2_1+sigR2_2+sigZ2_0+sigZ2_1+sigZ2_2)
    if mamamax > 40000:
        ax.set_ylim(math.sqrt(mimin), 250)
    else:
        ax.set_ylim(math.sqrt(mimin), math.sqrt(mamamax))
    ax.text(0.5, 0.9,'(c)', ha='center', va='center', transform= ax.transAxes)
    f.savefig(path + "/finalSVE.png")
    f.savefig(path + "/finalSVE.eps", format='eps')
    outSVE.close()


def plotFinalVertScale(path, r_gas, hzStar, num):
    '''Рисуем картинку с вертикальным звездным масштабом вдоль оси галактики - только для точек с газовой плотностью.'''

    outSig = open(path + '/finalHz.dat', 'w+')

    plt.figure(24+num)
    plt.plot(r_gas, hzStar, 'o-')
    plt.xlim(min(r_gas) - 10, max(r_gas) + 10)
    if max(hzStar) > 200:
        plt.ylim(0,200)
    else:
        plt.ylim(0,max(hzStar))
    plt.ylabel("$h_z^{s},\ arcsec$")
    plt.xlabel("$R,\ arcsec$")
    plt.savefig(path + "/finalHz.png")
    plt.savefig(path + "/finalHz.eps", format='eps')

    outSig.write('#R h_z^s\n')
    for r in zip(r_gas, hzStar):
        outSig.write(str(r[0]) + '  ' + str(r[1]) + '\n')

    outSig.close()




def plotPics(path, poly_star, poly_gas, delta_i, star_data, gas_data, incl, resolution, h_disc, r_eff_bulge,
             sigLosMajData, sigLosGaussParams, sigLosPolyParams, sigLosSinhParams, r_gas, Q1F, Q2FHydro, Q2FKinem,
             Q2F_withDisc, r_gas_and_star, sigR2, sigPhi2, sigZ2, hzStar):
    '''Промежуточная фугкция, чтобы не запускать все до конца - строит промежуточные картинки.'''