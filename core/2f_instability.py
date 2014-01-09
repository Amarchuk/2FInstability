__author__ = 'amarch'
# -*- coding: utf-8 -*-
from Galaxy import *
from RotationCurve import *

if __name__ == "__main__":

    ngc338 = Galaxy(name="NGC 338 (UGC 624)", path = "../data/ngc338", incl = "64.0", delta_incl = "7.5",
               description = "Inclination according to Zasov 2012, photometry in I band.",
               resolution = "311", image = "../data/ngc338/ngc338_SDSS.jpeg")
    ngc338.add_param("SDSS DR9 link",
                "http://skyserver.sdss3.org/dr9/en/tools/explore/obj.asp?ra=15.15164775&dec=30.66902519")
    ngc338.add_param("2MASS image", "../data/ngc338/ngc338_JHK.jpg")
    ngc338.print_info()

    stars_ma_rc = RotationCurve(ngc338.path + "/v_stars_ma.dat", "Stars MA RC NGC338")
    stars_ma_rc.print_info()

    fig = plt.figure(0)
    ngc338.plot()
    plt.show()