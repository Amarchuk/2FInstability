__author__ = 'amarch'
# -*- coding: utf-8 -*-

from numpy import *
import matplotlib.pyplot as plt
from scipy.optimize import *
import scipy
from scipy.special import *
import math
import inspect
import Tkinter
from PIL import ImageTk, Image
from utils import strutils as infoutils
from utils.bcolors import *
from RotationCurveHandler import *

class Galaxy():

    params = [
        "name", # galaxy name
        "path", # path to galaxy's data folder
        "incl", # inclination in deg
        "delta_incl", # variation in inclination
        "description", #galaxy notes
        "resolution", # resolution in pc/arcsec
        "image" # path to general image of galaxy
    ]

    def __init__(self, **params):
        self.name = None
        self.path = "."
        self.incl = None
        self.delta_incl = None
        self.description = ""
        self.resolution = None
        self.image = None
        self.star_rc = None
        self.gas_rc = None
        for param in self.params:
            setattr(self, param, params.pop(param, getattr(self, param)))

    def print_info(self, indent):
        infoutils.print_header("Galaxy", self.name, self.description, indent)
        for param in self.params:
            if param not in ("name","description"):
                infoutils.print_simple_param(1+indent, param, str(getattr(self, param)))
        if self.star_rc is not None:
            self.star_rc.print_info(indent)
        if self.gas_rc is not None:
            self.gas_rc.print_info(indent)


    def add_param(self, new_param, value=None):
        self.params.append(new_param)
        setattr(self, new_param, value)

    def plot(self):
        im = plt.imread(self.image)
        plt.imshow(im)
        plt.title(self.name)

    def plot_rcs(self):
        plt.title("Rotation curves for " + self.name)
        self.star_rc.plot('$V_{star}(R)$')
        self.gas_rc.plot('$V_{gas}(R)$', color='green')

    def initialize_handler(self):
        self.rc_handler = RotationCurveHandler(self.name, self)

    def handle_rcs(self, zero_point_star=(0,0), zero_point_gas=(0,0), incl=(None, None), gas_name ="",
                   star_poly_deg=0, gas_poly_deg=0,
                   star_fake_points=(), gas_fake_points=()):

        star_incl = incl[0] or self.incl
        gas_incl = incl[1] or self.incl
        self.rc_handler.get_corrected_star_ma_rc(zero_point_star, star_incl)
        for entity in star_fake_points:
            self.rc_handler.bended_star_ma_rc.add_fake_points(entity[0], entity[1])
        self.rc_handler.interpolate_poly_rc(self.rc_handler.bended_star_ma_rc,star_poly_deg)
        self.star_rc = self.rc_handler.bended_star_ma_rc

        self.rc_handler.get_corrected_gas_ma_rc(gas_name, zero_point_gas, gas_incl)
        for entity in gas_fake_points:
            self.rc_handler.bended_gas_ma_rc.add_fake_points(entity[0], entity[1])
        self.rc_handler.interpolate_poly_rc(self.rc_handler.bended_gas_ma_rc, gas_poly_deg)
        self.gas_rc = self.rc_handler.bended_gas_ma_rc


if __name__ == "__main__":
    b = Galaxy(name="NGC1", path = "some/path", incl = "60.0", delta_incl = "0.0",
               description = "Really cool galaxy\n with many-many\n cool observations available")
    b.print_info()
    b.add_param("SDSS Link", "http://link-to-sdss")
    b.print_info()

    # top = Tkinter.Tk()
    # top.title("Title bitches")
    # top.geometry("500x500")
    # # entry = Tkinter.Entry(top)
    # # entry.grid(column=0,row=0,sticky='EW')
    # # button = Tkinter.Button(top,text=u"Click me !")
    # # button.grid(column=1,row=0)
    # # label = Tkinter.Label(top,anchor="w",fg="white",bg="blue")
    # # label.grid(column=0,row=1,columnspan=2,sticky='EW')
    # # top.grid_columnconfigure(0,weight=1)
    # img = ImageTk.PhotoImage(Image.open("../incl_compar.gif"))
    # panel = Tkinter.Label(top, image = img)
    # panel.pack(side = "bottom", fill = "both", expand = "yes")
    # top.mainloop()