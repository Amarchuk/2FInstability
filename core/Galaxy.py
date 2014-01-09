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


class Galaxy():

    name = None
    path = "."
    incl = None
    delta_incl = None
    description = ""
    resolution = None
    image = None

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
        for param in self.params:
            setattr(self, param, params.pop(param, getattr(self, param)))

    def print_info(self):
        infoutils.print_header("Galaxy", self.name, self.description, 0)
        for param in self.params:
            if param not in ("name","description"):
                infoutils.print_simple_param(1, param, getattr(self, param))


    def add_param(self, new_param, value=None):
        self.params.append(new_param)
        setattr(self, new_param, value)

    def plot(self):
        im = plt.imread(self.image)
        plt.imshow(im)
        plt.title(self.name)


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