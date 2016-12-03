# -*- coding: utf-8 -*-

import yaml as ym
import numpy as np


class Cfg:

    def __init__(self):
        pass


class Srv:

    def __init__(self):
        pass

    def unpack(self, filename):
        with open(filename) as yml:
            return ym.load(yml)

    def dms2ddd(self, dms):
        degs = dms[0:2]
        mins = dms[3:5]
        secs = dms[6:len(dms)]
        degrees = float(degs) + float(mins) / 60 + float(secs) / 3600
        print 'DMS: ', degs, ':', mins, ":", secs
        return degrees

    def ellipsoid_params(self, ellipsoid_name):
        ell_library = Srv().unpack('ellipsoids.yml')
        elldict = {}
        if ellipsoid_name in ell_library:
            elldict['a'] = ell_library[ellipsoid_name]['a']
            elldict['alpha'] = 1/ell_library[ellipsoid_name]['ralpha']

    def transform_params(self, system_old, system_new):
        pass
