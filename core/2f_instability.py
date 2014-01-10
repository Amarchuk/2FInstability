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
    ngc338.add_param("Image with HI surf. dens.", "../data/ngc338/ugc624.gif")
    ngc338.add_param("SDSS DR9 whole image", "../data/ngc338/sdss_dr9_whole.jpg")
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

    imgs = []
    imgs.append((ngc338.name, ngc338.image))
    imgs.append(("Image with HI surf. dens.", getattr(ngc338, "Image with HI surf. dens.")))
    imgs.append(("2MASS image", getattr(ngc338, "2MASS image")))
    imgs.append(("SDSS DR9 whole image", getattr(ngc338, "SDSS DR9 whole image")))
    ngc338.plot_imgs_in_subplots(imgs)
    plt.show()

