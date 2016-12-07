# -*- coding: utf-8 -*-

import numpy as np
import yaml as ym


class Trans:

    def __init__(self):
        pass

    def ellparam(self, ellname):
        with open('ellipsoids.yml') as ells:
            ellipsoids = ym.load(ells)
        readyels = {}
        if ellipsoids.has_key(ellname):
            readyels['a'] = ellipsoids[ellname]['a']
            readyels['alpha'] = 1/ellipsoids[ellname]['ralpha']
        return readyels

    def transparams(self, ellname1, ellname2):
        stringtofind = ellname1 + '2' + ellname2
        print 'stringtofind', stringtofind
        with open('transes.yml') as trans:
            transs = ym.load(trans)
        reatrans = {}
        if transs.has_key(stringtofind):
            reatrans['m'] = eval(transs[stringtofind]['m'])
            reatrans['dx'] = eval(transs[stringtofind]['dx'])
            reatrans['dy'] = eval(transs[stringtofind]['dy'])
            reatrans['dz'] = eval(transs[stringtofind]['dz'])
            reatrans['wx'] = eval(transs[stringtofind]['wx'])
            reatrans['wy'] = eval(transs[stringtofind]['wy'])
            reatrans['wz'] = eval(transs[stringtofind]['wz'])
        return reatrans

ell1 = 'WGS84'
ell2 = 'PZ90'

print Trans().ellparam(ell1)
print Trans().ellparam(ell2)
print Trans().transparams(ell1, ell2)









