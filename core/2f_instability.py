__author__ = 'amarch'
# -*- coding: utf-8 -*-
from Galaxy import *

if __name__ == "__main__":
    ngc338 = Galaxy(name="NGC 338 (UGC 624)", path="../data/ngc338", incl=64.0, delta_incl=7.5,
                    description="Inclination according to Zasov 2012, photometry in I band.",
                    resolution="311", image="../data/ngc338/ngc338_SDSS.jpeg")
    ngc338.add_param("SDSS DR9 link",
                     "http://skyserver.sdss3.org/dr9/en/tools/explore/obj.asp?ra=15.15164775&dec=30.66902519")
    ngc338.add_param("2MASS image", "../data/ngc338/ngc338_JHK.jpg")
    ngc338.print_info(0)

    stars_ma_rc = RotationCurve(ngc338.path + "/v_stars_ma.dat", "Stars MA RC NGC338")
    stars_ma_rc.print_info(0)

    ngc338.initialize_handler()
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

    handler = ngc338.rc_handler

    # handler.bended_star_ma_rc.plot("")
    # handler.bended_gas_ma_rc.plot("")
    # xx = arange(min(handler.bended_star_ma_rc.radii()), max(handler.bended_star_ma_rc.radii()) + 25, 0.1)
    # plt.plot(xx, handler.bended_star_ma_rc.poly_fit(xx), '-')
    # plt.plot(xx, handler.bended_gas_ma_rc.poly_fit(xx), '-')
    # plt.plot(zip(*handler.bended_star_ma_rc.fake_data_points)[0],
    #          zip(*handler.bended_star_ma_rc.fake_data_points)[1], 'x')
    # plt.plot(zip(*handler.bended_gas_ma_rc.fake_data_points)[0],
    #          zip(*handler.bended_gas_ma_rc.fake_data_points)[1], 'x')
    # plt.show()

    handler.plot_all_uncorrected_centralize()
    plt.show()