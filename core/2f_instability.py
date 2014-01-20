__author__ = 'amarch'
# -*- coding: utf-8 -*-
from Galaxy import *
from utils import plotutils as plotutils
from VelocityDispersion import *


if __name__ == "__main__":
    ngc338 = Galaxy(name="NGC 338 (UGC 624)", path="../data/ngc338", incl=64.0, delta_incl=7.5,
                    description="Inclination according to Zasov 2012, photometry in I band.",
                    resolution=311, image="../data/ngc338/ngc338_SDSS.jpeg")
    ngc338.add_param("SDSS DR9 link",
                     "http://skyserver.sdss3.org/dr9/en/tools/explore/obj.asp?ra=15.15164775&dec=30.66902519")
    ngc338.add_img("2MASS image JHK", "../data/ngc338/ngc338_JHK.jpg")
    ngc338.add_img("Image with HI surf. dens.", "../data/ngc338/ugc624.gif")
    ngc338.add_img("SDSS DR9 whole image", "../data/ngc338/sdss_dr9_whole.jpg")
    ngc338.print_info(0)

    stars_ma_rc = RotationCurve(ngc338.path + "/v_stars_ma.dat", "Stars MA RC NGC338")
    stars_ma_rc.print_info(0)

    ngc338.initialize_rc_handler()
    ngc338.rc_handler.set_stellar_ma_rc(RotationCurve("../data/ngc338/v_stars_ma.dat", "Stars MA RC NGC338"))
    ngc338.rc_handler.set_stellar_mi_rc(RotationCurve("../data/ngc338/v_stars_mi.dat", "Stars MI RC NGC338"))
    ngc338.rc_handler.add_gas_rc(RotationCurve("../data/ngc338/v_gas_ma.dat", "Gas MA RC NGC338"), "v_gas_ma")
    ngc338.rc_handler.add_gas_rc(RotationCurve("../data/ngc338/v_gas_Katkov.dat", "Gas Katkov RC NGC338"),
                                 "v_gas_Katkov")
    ngc338.rc_handler.add_gas_rc(RotationCurve("../data/ngc338/v_gas_WSRT.dat", "Gas WSRT RC NGC338"), "v_gas_WSRT")
    ngc338.rc_handler.add_gas_rc(RotationCurve("../data/ngc338/v_gas_Court.dat", "Gas Court RC NGC338"), "v_gas_Court")

    ngc338.handle_rcs(zero_point_star=(4759.845, 0.0), zero_point_gas=(4759.845, 0.0), gas_name="v_gas_ma",
                      star_poly_deg=25, gas_poly_deg=25, star_fake_points=(((49, 238, 5), 50),),
                      gas_fake_points=(((32, 285, 1), 5), ((46, 268, 1), 54)))

    ngc338.initialize_sig_los_handler()
    ngc338.sig_handler.set_sig_ma(VelocityDispersion("../data/ngc338/v_stars_ma.dat", "Sig MA los NGC338"))
    ngc338.sig_handler.set_sig_mi(VelocityDispersion("../data/ngc338/v_stars_mi.dat", "Sig MI los NGC338"))
    ngc338.handle_sig_los(sig_ma_deg=10, sig_mi_deg=8,
                   sig_ma_fake_points=[((15.0, 84.0, 1.0), 90, 114.0)],
                   sig_mi_fake_points=[((42.7500, 87.70, 10.28), 10, 10000)])

    ngc338.plot_sig_los()
    plt.show()

    ngc338.print_info(0)
    # ngc338.sig_handler.plot_two_in_one()
    # plt.show()
    #
    # ngc338.sig_handler.sig_ma.plot('$\sigma_{los}^{maj}$')
    # plt.show()

    # ngc1167 = Galaxy(name="NGC 1167 (UGC 2487)", path="../data/ngc1167", incl=36.0, delta_incl=2.0,
    #                 description="Photometry in R band.",
    #                 resolution=330, image="../data/ngc1167/ngc1167_SDSS.jpeg")
    # ngc1167.add_param("SDSS DR9 link",
    #                  "http://skyserver.sdss3.org/dr9/en/tools/explore/obj.asp?ra=45.42639624&dec=35.20561448")
    # ngc1167.add_img("2MASS image JHK", "../data/ngc1167/ngc1167_JHK.jpg")
    # ngc1167.add_img("Image with HI surf. dens.", "../data/ngc1167/ugc2487.gif")
    # ngc1167.print_info(0)
    #
    # ngc1167.initialize_handler()
    # ngc1167.rc_handler.set_stellar_ma_rc(RotationCurve("../data/ngc1167/v_stars_ma.dat", "Stars MA RC NGC1167"))
    # ngc1167.rc_handler.set_stellar_mi_rc(RotationCurve("../data/ngc1167/v_stars_mi.dat", "Stars MI RC NGC1167"))
    # ngc1167.rc_handler.add_gas_rc(RotationCurve("../data/ngc1167/v_gas_ma.dat", "Gas MA RC NGC1167"), "v_gas_ma")
    # ngc1167.rc_handler.add_gas_rc(RotationCurve("../data/ngc1167/v_gas_noord.dat", "Gas Noord RC NGC1167"),
    #                              "v_gas_noord")
    # ngc1167.rc_handler.add_gas_rc(RotationCurve("../data/ngc1167/v_gas_WSRT.dat", "Gas WSRT RC NGC1167"), "v_gas_WSRT")
    # ngc1167.rc_handler.add_gas_rc(RotationCurve("../data/ngc1167/v_gas_struve.dat", "Gas Struve RC NGC1167"), "v_gas_struve")
    #
    # ngc1167.handle_rcs(zero_point_star=ngc1167.rc_handler.find_zero_point(ngc1167.rc_handler.star_ma_rc),
    #                    zero_point_gas=ngc1167.rc_handler.find_zero_point(ngc1167.rc_handler.gas_rcs["v_gas_ma"]),
    #                    gas_name="v_gas_ma",
    #                   star_poly_deg=15, gas_poly_deg=8, star_fake_points=(((36, 340, 2), 50),),
    #                   gas_fake_points=())

