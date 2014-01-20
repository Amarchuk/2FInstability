__author__ = 'amarch'
# -*- coding: utf-8 -*-
from Galaxy import *
from utils import plotutils as plotutils
from VelocityDispersion import *


if __name__ == "__main__":
    # ngc338 = Galaxy(name="NGC 338 (UGC 624)", path="../data/ngc338", incl=64.0, delta_incl=7.5,
    #                 description="Inclination according to Zasov 2012, photometry in I band.",
    #                 resolution=311.0, image="../data/ngc338/ngc338_SDSS.jpeg")
    # ngc338.add_param("SDSS DR9 link",
    #                  "http://skyserver.sdss3.org/dr9/en/tools/explore/obj.asp?ra=15.15164775&dec=30.66902519")
    # ngc338.add_img("2MASS image JHK", "../data/ngc338/ngc338_JHK.jpg")
    # ngc338.add_img("Image with HI surf. dens.", "../data/ngc338/ugc624.gif")
    # ngc338.add_img("SDSS DR9 whole image", "../data/ngc338/sdss_dr9_whole.jpg")
    # ngc338.print_info(0)
    #
    # stars_ma_rc = RotationCurve(ngc338.path + "/v_stars_ma.dat", "Stars MA RC NGC338")
    # stars_ma_rc.print_info(0)
    #
    # ngc338.initialize_rc_handler()
    # ngc338.rc_handler.set_stellar_ma_rc(RotationCurve("../data/ngc338/v_stars_ma.dat", "Stars MA RC NGC338"))
    # ngc338.rc_handler.set_stellar_mi_rc(RotationCurve("../data/ngc338/v_stars_mi.dat", "Stars MI RC NGC338"))
    # ngc338.rc_handler.add_gas_rc(RotationCurve("../data/ngc338/v_gas_ma.dat", "Gas MA RC NGC338"), "v_gas_ma")
    # ngc338.rc_handler.add_gas_rc(RotationCurve("../data/ngc338/v_gas_Katkov.dat", "Gas Katkov RC NGC338"),
    #                              "v_gas_Katkov")
    # ngc338.rc_handler.add_gas_rc(RotationCurve("../data/ngc338/v_gas_WSRT.dat", "Gas WSRT RC NGC338"), "v_gas_WSRT")
    # ngc338.rc_handler.add_gas_rc(RotationCurve("../data/ngc338/v_gas_Court.dat", "Gas Court RC NGC338"), "v_gas_Court")
    #
    # ngc338.handle_rcs(zero_point_star=(4759.845, 0.0), zero_point_gas=(4759.845, 0.0), gas_name="v_gas_ma",
    #                   star_poly_deg=25, gas_poly_deg=25, star_fake_points=(((49.0, 238.0, 5.0), 50),),
    #                   gas_fake_points=(((32.0, 285.0, 1.0), 5), ((46.0, 268.0, 1.0), 54)))
    #
    # ngc338.initialize_sig_los_handler()
    # ngc338.sig_handler.set_sig_ma(VelocityDispersion("../data/ngc338/v_stars_ma.dat", "Sig MA los NGC338"))
    # ngc338.sig_handler.set_sig_mi(VelocityDispersion("../data/ngc338/v_stars_mi.dat", "Sig MI los NGC338"))
    # ngc338.handle_sig_los(sig_ma_deg=10, sig_mi_deg=8,
    #                sig_ma_fake_points=[((15.0, 84.0, 1.0), 90, 114.0)],
    #                sig_mi_fake_points=[((42.7500, 87.70, 10.28), 10, 10000)])
    #
    # # ngc338.plot_sig_los()
    # # plt.show()
    # #
    # # ngc338.initialize_sve_handler()
    # # ngc338.sve_handler.plot_sigZ2_to_sigR2()
    # # plt.show()
    #
    # ngc1167 = Galaxy(name="NGC 1167 (UGC 2487)", path="../data/ngc1167", incl=36.0, delta_incl=2.0,
    #                 description="Photometry in R band.",
    #                 resolution=330.0, image="../data/ngc1167/ngc1167_SDSS.jpeg")
    # ngc1167.add_param("SDSS DR9 link",
    #                  "http://skyserver.sdss3.org/dr9/en/tools/explore/obj.asp?ra=45.42639624&dec=35.20561448")
    # ngc1167.add_img("2MASS image JHK", "../data/ngc1167/ngc1167_JHK.jpg")
    # ngc1167.add_img("Image with HI surf. dens.", "../data/ngc1167/ugc2487.gif")
    # ngc1167.print_info(0)
    #
    # ngc1167.initialize_rc_handler()
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
    #                   star_poly_deg=15, gas_poly_deg=8, star_fake_points=(((36.0, 340.0, 2.0), 50),),
    #                   gas_fake_points=())
    #
    # ngc1167.initialize_sig_los_handler()
    # ngc1167.sig_handler.set_sig_ma(VelocityDispersion("../data/ngc1167/v_stars_ma.dat", "Sig MA los NGC1167"))
    # ngc1167.sig_handler.set_sig_mi(VelocityDispersion("../data/ngc1167/v_stars_mi.dat", "Sig MI los NGC1167"))
    # ngc1167.handle_sig_los(sig_ma_deg=10, sig_mi_deg=8,
    #                sig_ma_fake_points=[((7.0, 220.0, 1.0), 90, 43.0)],
    #                sig_mi_fake_points=[((28.0, 143.5, 1.0), 40, 140.0)])
    #

    # ngc1167.plot_sig_los()
    # plt.show()
    #
    # ngc1167.initialize_sve_handler()
    # ngc1167.sve_handler.plot_sigZ2_to_sigR2()
    # plt.show()


    # ngc2273 = Galaxy(name="NGC 2273 (UGC 3546)", path="../data/ngc2273", incl=55.0, delta_incl=3.0,
    #                 description="Photometry in R band.",
    #                 resolution=130.0, image="../data/ngc2273/ngc2273_JHK.jpg")
    # ngc2273.add_img("Image with HI surf. dens.", "../data/ngc2273/ugc3546.gif")
    #
    # ngc2273.initialize_rc_handler()
    # ngc2273.rc_handler.set_stellar_ma_rc(RotationCurve("../data/ngc2273/v_stars_ma.dat", "Stars MA RC ngc2273"))
    # ngc2273.rc_handler.set_stellar_mi_rc(RotationCurve("../data/ngc2273/v_stars_mi.dat", "Stars MI RC ngc2273"))
    # ngc2273.rc_handler.add_gas_rc(RotationCurve("../data/ngc2273/v_gas_ma.dat", "Gas MA RC ngc2273"), "v_gas_ma")
    # ngc2273.rc_handler.add_gas_rc(RotationCurve("../data/ngc2273/v_gas_WSRT.dat", "Gas WSRT RC ngc2273"), "v_gas_WSRT")
    #
    # ngc2273.handle_rcs(zero_point_star=(0,0),
    #                    zero_point_gas=(0,0),
    #                    gas_name="v_gas_ma",
    #                   star_poly_deg=15, gas_poly_deg=8, star_fake_points=(((36.0, 155.0, 2.0), 80),),
    #                   gas_fake_points=())
    #
    # ngc2273.initialize_sig_los_handler()
    # ngc2273.sig_handler.set_sig_ma(VelocityDispersion("../data/ngc2273/v_stars_ma.dat", "Sig MA los ngc2273"))
    # ngc2273.sig_handler.set_sig_mi(VelocityDispersion("../data/ngc2273/v_stars_mi.dat", "Sig MI los ngc2273"))
    # ngc2273.handle_sig_los(sig_ma_deg=10, sig_mi_deg=10,
    #                sig_ma_fake_points=[((21.4, 42.0, 1.0), 87, 100.0)],
    #                sig_mi_fake_points=[((37.54, 71.0, 1.0), 87, 187.0), ((0.0, 120.0, 1.0), 2, 1000.0)])
    #
    # ngc2273.plot_sig_los()
    # plt.show()
    #
    # ngc2273.initialize_sve_handler()
    # ngc2273.sve_handler.plot_sigZ2_to_sigR2()
    # plt.show()

    # ngc2985 = Galaxy(name="NGC 2985 (UGC 5253)", path="../data/ngc2985", incl=36.0, delta_incl=3.0,
    #                 description="Photometry in R band.",
    #                 resolution=102.0, image="../data/ngc2985/ngc2985_JHK.jpg")
    # ngc2985.add_img("Image with HI surf. dens.", "../data/ngc2985/ugc5253.gif")
    #
    # ngc2985.initialize_rc_handler()
    # ngc2985.rc_handler.set_stellar_ma_rc(RotationCurve("../data/ngc2985/v_stars_noord.dat", "Stars MA RC Noord ngc2985"))
    # ngc2985.rc_handler.set_stellar_mi_rc(RotationCurve("../data/ngc2985/v_stars_mi.dat", "Stars MI RC ngc2985"))
    # ngc2985.rc_handler.add_gas_rc(RotationCurve("../data/ngc2985/v_gas_ma.dat", "Gas MA RC ngc2985"), "v_gas_ma")
    # ngc2985.rc_handler.add_gas_rc(RotationCurve("../data/ngc2985/v_gas_Halpha.dat", "Gas Halpha RC ngc2985"), "v_gas_Halpha")
    # ngc2985.rc_handler.add_gas_rc(RotationCurve("../data/ngc2985/v_gas_WSRT.dat", "Gas WSRT RC ngc2985"), "v_gas_WSRT")
    #
    # ngc2985.handle_rcs(zero_point_star=(0,0),
    #                    zero_point_gas=(0,0),
    #                    gas_name="v_gas_ma",
    #                   star_poly_deg=25, gas_poly_deg=25, star_fake_points=(((44.0, 226.0, 1.0), 60),),
    #                   gas_fake_points=(((71.0, 240.0, 1.0), 34),))
    #
    # # ngc2985.plot_rcs()
    # # plt.show()
    #
    # ngc2985.initialize_sig_los_handler()
    # ngc2985.sig_handler.set_sig_ma(VelocityDispersion("../data/ngc2985/s_stars_maN.dat", "Sig MA los ngc2985"))
    # ngc2985.sig_handler.set_sig_mi(VelocityDispersion("../data/ngc2985/s_stars_miN.dat", "Sig MI los ngc2985"))
    # ngc2985.handle_sig_los(sig_ma_deg=8, sig_mi_deg=12,
    #                sig_ma_fake_points=[((20.0, 81.0, 1.0), 150, 265.0)],
    #                sig_mi_fake_points=[((15.0, 100.0, 1.0), 100, 250.0), ])
    #
    # ngc2985.plot_sig_los()
    # plt.show()
    #
    # ngc2985.initialize_sve_handler()
    # ngc2985.sve_handler.plot_sigZ2_to_sigR2()
    # plt.show()


    # ngc3898 = Galaxy(name="NGC 2985 (UGC 6787)", path="../data/ngc3898", incl=70.0, delta_incl=3.0,
    #                 description="Photometry in R band.",
    #                 resolution=92.0, image="../data/ngc3898/ngc3898_SDSS.jpeg")
    # ngc3898.add_img("Image with HI surf. dens.", "../data/ngc3898/ugc6787.gif")
    # ngc3898.add_img("2MASS image JHK", "../data/ngc3898/ngc3898_JHK.jpg")
    #
    # ngc3898.initialize_rc_handler()
    # ngc3898.rc_handler.set_stellar_ma_rc(RotationCurve("../data/ngc3898/v_stars_noord.dat", "Stars MA RC Noord ngc3898"))
    # ngc3898.rc_handler.set_stellar_mi_rc(RotationCurve("../data/ngc3898/v_stars_mi.dat", "Stars MI RC ngc3898"))
    # ngc3898.rc_handler.add_gas_rc(RotationCurve("../data/ngc3898/v_gas_ma.dat", "Gas MA RC ngc3898"), "v_gas_ma")
    # ngc3898.rc_handler.add_gas_rc(RotationCurve("../data/ngc3898/v_gas_Ha.dat", "Gas Halpha RC ngc3898"), "v_gas_Halpha")
    # ngc3898.rc_handler.add_gas_rc(RotationCurve("../data/ngc3898/v_gas_WSRT.dat", "Gas WSRT RC ngc3898"), "v_gas_WSRT")
    #
    # ngc3898.handle_rcs(zero_point_star=(0,0),
    #                    zero_point_gas=(0,0),
    #                    gas_name="v_gas_ma",
    #                   star_poly_deg=12, gas_poly_deg=25,
    #                   star_fake_points=(((16.5, 156.0, 2.0), 126),((8.5, 189.0, 1.0), 2)),
    #                   gas_fake_points=(((5.0, 274.0, 5.0), 5),((189.0, 240.0, 5.0), 30), ((0.0, 0.0, 1.0), 1), ((22.0, 217.0, 1.0), 7)))
    #
    #
    # # ngc3898.plot_rcs()
    # # plt.show()
    #
    # ngc3898.initialize_sig_los_handler()
    # ngc3898.sig_handler.set_sig_ma(VelocityDispersion("../data/ngc3898/s_stars_maN.dat", "Sig MA los ngc3898"))
    # ngc3898.sig_handler.set_sig_mi(VelocityDispersion("../data/ngc3898/s_stars_miN.dat", "Sig MI los ngc3898"))
    # ngc3898.handle_sig_los(sig_ma_deg=8, sig_mi_deg=12,
    #                sig_ma_fake_points=[((8.8, 147.0, 1.0), 120, 360.0)],
    #                sig_mi_fake_points=[((24.0, 183.0, 1.0), 100, 300.0), ])
    #
    # ngc3898.plot_sig_los()
    # plt.show()
    #
    # ngc3898.initialize_sve_handler()
    # ngc3898.sve_handler.plot_sigZ2_to_sigR2()
    # plt.show()

    ngc5533 = Galaxy(name="NGC 5533 (UGC 9133)", path="../data/ngc5533", incl=53.0, delta_incl=2.5,
                    description="Photometry in R band.",
                    resolution=260.0, image="../data/ngc5533/ngc5533_SDSS.jpeg")
    ngc5533.add_img("Image with HI surf. dens.", "../data/ngc5533/ugc9133.gif")
    ngc5533.add_img("2MASS image JHK", "../data/ngc5533/ngc5533_JHK.jpg")

    ngc5533.initialize_rc_handler()
    ngc5533.rc_handler.set_stellar_ma_rc(RotationCurve("../data/ngc5533/v_stars_noord.dat", "Stars MA RC Noord ngc5533"))
    ## ngc5533.rc_handler.set_stellar_mi_rc(RotationCurve("../data/ngc5533/v_stars_mi.dat", "Stars MI RC ngc5533"))
    ngc5533.rc_handler.add_gas_rc(RotationCurve("../data/ngc5533/v_gas_ma.dat", "Gas MA RC ngc5533"), "v_gas_ma")
    ngc5533.rc_handler.add_gas_rc(RotationCurve("../data/ngc5533/v_gas_Court.dat", "Gas Court RC ngc5533"), "v_gas_Court")
    ngc5533.rc_handler.add_gas_rc(RotationCurve("../data/ngc5533/v_gas_WSRT.dat", "Gas WSRT RC ngc5533"), "v_gas_WSRT")
    ngc5533.rc_handler.add_gas_rc(RotationCurve("../data/ngc5533/v_gas_SBV.dat", "Gas SBV RC ngc5533"), "v_gas_SBV")

    ngc5533.handle_rcs(zero_point_star=(0,0),
                       zero_point_gas=(0,0),
                       gas_name="v_gas_ma",
                      star_poly_deg=12, gas_poly_deg=5,
                      star_fake_points=(((32.0, 214.0, 1.0), 60),),
                      gas_fake_points=(((37.6,287.5,0.0),1),))

    # ngc5533.plot_rcs()
    # plt.show()

    ngc5533.initialize_sig_los_handler()
    ngc5533.sig_handler.set_sig_ma(VelocityDispersion("../data/ngc5533/s_stars_maN.dat", "Sig MA los ngc5533"))
    ngc5533.sig_handler.set_sig_mi(VelocityDispersion("../data/ngc5533/s_stars_miN.dat", "Sig MI los ngc5533"))
    ngc5533.handle_sig_los(sig_ma_deg=8, sig_mi_deg=12,
                   sig_ma_fake_points=[((19.7, 106.0, 1.0), 90, 100.0)],
                   sig_mi_fake_points=[((18.7, 123.0, 1.0), 100, 200.0), ])


    ngc5533.plot_sig_los()
    plt.show()

    ngc5533.initialize_sve_handler()
    ngc5533.sve_handler.plot_sigZ2_to_sigR2()
    plt.show()