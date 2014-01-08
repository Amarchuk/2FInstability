__author__ = 'amarch'
# -*- coding: utf-8 -*-

class bcolors():

    def __init__(self):
        pass

    PURPLE = '\033[95m'
    HEADER = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    BOLD = "\033[1m"
    ENDC = '\033[0m'

    def disable(self):
        self.PURPLE = ''
        self.HEADER = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''