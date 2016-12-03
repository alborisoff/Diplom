# -*- coding: utf-8 -*-

import numpy as np
import yaml as ym


class Cfg:

    def __init__(self):
        pass

    def unpack(self, filename):
        with open(filename) as params:
            return ym.load(params)


class Trans:

    def __init__(self):
        pass

    def ellparam(self, ellname):
        ellipsoids = Cfg().unpack('ellipsoids.yml')
        readyels = {}
        if ellipsoids.has_key(ellname):
            readyels['a'] = ellipsoids[ellname]['a']
            readyels['alpha'] = 1/ellipsoids[ellname]['ralpha']
        return readyels

    def transparams(self, ellname1, ellname2):
        stringtofind = ellname1 + '2' + ellname2
        print 'stringtofind', stringtofind
        transs = Cfg().unpack('transes.yml')
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

    def elltorect(self, ellname, ellcoords):
        ellpars = self.ellparam(ellname)
        a = ellpars['a']
        alpha = ellpars['alpha']
        print 'a = ', a
        print 'alpha = ', alpha
        B = ellcoords[0]
        L = ellcoords[1]
        H = ellcoords[2]
        print 'B = ', B
        print 'L = ', L
        print 'H = ', H
        sinB = np.sin(np.radians(B))
        cosB = np.cos(np.radians(B))
        sinL = np.sin(np.radians(L))
        cosL = np.cos(np.radians(L))
        sinsqB = sinB**2
        e2 = 2*alpha - alpha**2
        subrad = 1 - e2*sinsqB
        print 'e2 = ', e2
        print 'subrad = ', subrad
        N = a/(subrad**0.5)
        print 'N = ', N
        X = (N + H)*cosB*cosL
        Y = (N + H)*cosB*sinL
        Z = (H + N)*sinB - e2*N*sinB
        return [X, Y, X]

    def dms2degrees(self, dms):
        degs = dms[0:2]
        mins = dms[3:5]
        secs = dms[6:len(dms)]
        degrees = float(degs) + float(mins)/60 + float(secs)/3600
        print 'DMS: ', degs, ':', mins, ":", secs
        return degrees

    def helmertrectangle(self, ellname1, ellname2, ell1_coordinates):
        rectangle = self.elltorect(ellname1, ell1_coordinates)
        traparam = self.transparams(ellname1, ellname2)
        print 'rectangle = ', rectangle
        print 'traparam = ', traparam
        m = traparam['m']
        dx = traparam['dx']
        dy = traparam['dy']
        dz = traparam['dz']
        wx = traparam['wx']
        wy = traparam['wy']
        wz = traparam['wz']
        matromegas = np.matrix([[1, wz, -wy],
                                [-wz, 1, wx],
                                [wy, -wx, 1]])
        matrdeltas = np.matrix([[dx],
                                [dy],
                                [dz]])
        oldcoords = np.matrix([[float(rectangle[0])],
                               [float(rectangle[1])],
                               [float(rectangle[2])]])
        newcoords = (1 + m)*matromegas*oldcoords + matrdeltas
        return newcoords

ell1 = 'WGS84'
ell2 = 'Kras40'
ell_B_string = '55 51 6.47082'
ell_L_string = '37 26 11.16388'
ell_H_string = '189.000'
ell_B = Trans().dms2degrees(ell_B_string)
ell_L = Trans().dms2degrees(ell_L_string)
ell_H = float(ell_H_string)
ellcoordes = [ell_B, ell_L, ell_H]

print Trans().elltorect(ell1, ellcoordes)
print Trans().helmertrectangle(ell1, ell2, ellcoordes)