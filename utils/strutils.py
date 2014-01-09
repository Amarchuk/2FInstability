__author__ = 'amarch'
# -*- coding: utf-8 -*-

import numpy
import sys
from bcolors import *


def print_list_summary(base_indent, title, list, description=None):
    """
    Print base statistics of list or tuple
    @param description: optional description
    @param base_indent: int indent in spaces
    @param title: title for list
    @param list: numerical data for info
    """
    indent_level = " " * base_indent * 4
    next_indent_level = " " * (base_indent*4 + 4)
    print indent_level + bcolors.OKGREEN + title + bcolors.ENDC + ":"
    if description != None:
        print indent_level + bcolors.BOLD + "Description: " + bcolors.ENDC +\
              str.replace(description, "\n", "\n" + next_indent_level)
    print next_indent_level + "points count: " + str(list.__len__())
    print next_indent_level + "range: <" + str(numpy.array(list).min()) + " : " + str(numpy.array(list).max()) + ">"
    print next_indent_level + "mean: " + str(numpy.array(list).mean()) + " median: " + str(
        numpy.median(numpy.array(list)))


def print_header(title, name, description=None, indent=0):
    indent_level = " " * indent * 4
    next_indent_level = " " * (indent + 4)
    print indent_level + bcolors.HEADER +  title + ": " + bcolors.ENDC + name
    if description != None:
        print indent_level + bcolors.HEADER + "Description: " + bcolors.ENDC \
              + str.replace(description, "\n", "\n" + next_indent_level)


def print_simple_param(indent, name, value):
    indent_level = " " * indent * 4
    if value != None:
        print indent_level + bcolors.OKGREEN +  name + ": " + bcolors.ENDC + str(value)
    else:
        print indent_level + bcolors.OKGREEN +  name + ": " + bcolors.FAIL + "None" + bcolors.ENDC