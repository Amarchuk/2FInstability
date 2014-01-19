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
from VelocityDispersionHandler import *

class Galaxy():

    def __init__(self, **params):
        self.params = [
            "name", # galaxy name
            "path", # path to galaxy's data folder
            "incl", # inclination in deg
            "delta_incl", # variation in inclination
            "description", #galaxy notes
            "resolution", # resolution in pc/arcsec
            "image" # path to general image of galaxy
        ]
        self.name = None
        self.path = "."
        self.incl = None
        self.delta_incl = None
        self.description = ""
        self.resolution = None
        self.image = None
        self.star_rc = None
        self.gas_rc = None
        self.imgs = []
        for param in self.params:
            setattr(self, param, params.pop(param, getattr(self, param)))
        if self.image is not None:
            self.imgs.append((self.name, self.image))

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

    def add_img(self, img_desc, img_path):
        self.imgs.append((img_desc, img_path))
        self.params.append(img_desc)
        setattr(self, img_desc, img_path)

    def plot(self):
        im = plt.imread(self.image)
        plt.imshow(im)
        plt.title(self.name)

    def plot_rcs(self):
        plt.title("Rotation curves for " + self.name)
        self.star_rc.plot('$V_{star}(R)$')
        self.gas_rc.plot('$V_{gas}(R)$', color='green')

    def plot_own_imgs(self):
        self.plot_imgs_in_subplots(self.imgs)

    def plot_imgs_in_subplots(self, imgs):
        '''
        plot images 2 per row,
        if odd, span first img
        '''
        imgs_count = imgs.__len__()
        curr_img_number = 0
        for img_and_descr in imgs:
            img = plt.imread(img_and_descr[1])
            descr = img_and_descr[0]
            if imgs_count%2 == 1 and curr_img_number == 0:
                plt.subplot((imgs_count/2 + imgs_count%2), 1, curr_img_number+1)
                curr_img_number += 2
            else:
                plt.subplot((imgs_count/2 + imgs_count%2), 2, curr_img_number+1)
                curr_img_number += 1
            plt.imshow(img)
            plt.title(descr)
            plt.tick_params(axis='x', which='both',bottom='off',top='off', labelbottom='off')
            plt.tick_params(axis='y', which='both',bottom='off',top='off',  labelleft='off')
            plt.subplots_adjust(left=0.01, right=0.991, top=0.95, bottom=0.01, hspace=0.1)



    def initialize_rc_handler(self):
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

    def initialize_sig_los_handler(self):
        self.sig_handler = VelocityDispersionHandler(self.name, self)

    def handle_sig_los(self, sig_ma_deg=0, sig_mi_deg=0,
                   sig_ma_fake_points=(), sig_mi_fake_points=()):

        for entity in sig_ma_fake_points:
            self.sig_handler.sig_ma.add_fake_exp_points(entity[0], entity[1], expscale=entity[2])
        self.sig_handler.interpolate_poly_sig(self.sig_handler.sig_ma, sig_ma_deg)
        self.sig_los_ma = self.sig_handler.sig_ma

        for entity in sig_mi_fake_points:
            self.sig_handler.sig_mi.add_fake_exp_points(entity[0], entity[1], expscale=entity[2])
        self.sig_handler.sig_mi.expand_minor(self.incl)
        self.sig_handler.interpolate_poly_sig(self.sig_handler.sig_mi, sig_mi_deg)
        self.sig_los_mi = self.sig_handler.sig_mi

    def plot_sig_los(self):
        self.sig_handler.plot_two_in_subplots()