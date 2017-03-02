# -*- coding: utf-8 -*-

import numpy as np
import json


class Srv():

    def __init__(self):
        self.elldict = self.unpack("ellipsoids.json")
        self.transparams = self.unpack("transparams.json")

    def unpack(self, filename):  # Функция для распаковки внешних JSON-файлов в переменные для дальнейшей работы
        return json.loads(open(filename).read())

    def ellparams(self, ell_name):  # Получение полной информации об эллипсоиде по его названию
        for one_ellipsoid in self.elldict:
            if one_ellipsoid["name"] == ell_name:
                a = one_ellipsoid["a"]
                alpha = 1.0/one_ellipsoid["ralpha"]
                e2 = 2*alpha - alpha**2
                e12 = e2/(1 - e2)
                b = a*((1 - e2)**0.5)
                ready_dict = {"name": one_ellipsoid["name"],
                              "a": a,
                              "b": b,
                              "alpha": alpha,
                              "e2": e2,
                              "e12": e12}
                return ready_dict

    # По названию системы координат подгружается эллипсоид, на котором она основана
    def ellipsoidbysystem(self, system_name):
        for one_ellipsoid in self.elldict:
            if system_name in one_ellipsoid["systems"]:
                return one_ellipsoid["name"]

    def n_big(self, ellipsoid, lat):
        a = self.ellparams(ellipsoid)["a"]
        dnmtr = self.bigsqrt(ellipsoid, lat)
        return a/dnmtr

    def bigsqrt(self, ellipsoid, lat):
        e2 = self.ellparams(ellipsoid)["e2"]
        sinqB = (np.sin(lat)) ** 2
        radix = (1 - e2 * sinqB) ** 0.5
        return radix


    # Преобразование эллипсоидальных координат из градусов, минут и секунд
    # в градусы и доли градусов.
    def dms2ddd(self, dms):
        dmssplit = dms.split(' ')
        degrees = float(dmssplit[0]) + float(dmssplit[1]) / 60 + float(dmssplit[2]) / 3600
        return degrees


# Функции для внутреннего преобразования координат
class IntTrans:

    # Преобразование геодезических координат в прямоугольные
    def geodez2geocentr(self, ellipsoid, b_big, l_big, h_big):
        n_big = Srv().n_big(ellipsoid, b_big)
        e2 = Srv().ellparams(ellipsoid)["e2"]
        sinB = np.sin(np.radians(b_big))
        cosB = np.cos(np.radians(b_big))
        sinL = np.sin(np.radians(l_big))
        cosL = np.cos(np.radians(l_big))
        x_big = (n_big + h_big)*cosB*cosL
        y_big = (n_big + h_big)*cosB*sinL
        z_big = sinB*(n_big*(1 - e2) + h_big)
        return [x_big, y_big, z_big]

    def geocentr2geodez(self, ellipsoid, x_big, y_big, z_big):
        b_big, l_big, h_big, b, c = 0, 0, 0, 0, 0
        a = Srv().ellparams(ellipsoid)["a"]
        e2 = Srv().ellparams(ellipsoid)["e2"]
        r = (x_big**2 + y_big**2)**0.5
        if not(x_big == 0 and y_big == 0 and z_big == 0):
            if r == 0:
                b_big = (np.pi*z_big)/abs(2*z_big)
                l_big = 0
                h_big = z_big*np.sin(b_big) - a*Srv().bigsqrt(ellipsoid, np.degrees(b_big))
            else:
                l_big_pre = abs(np.arcsin(y_big/r))
                if x_big > 0 and y_big < 0:
                    l_big = 2*np.pi - l_big_pre
                elif x_big < 0 and y_big < 0:
                    l_big = np.pi + l_big_pre
                elif x_big < 0 and y_big > 0:
                    l_big = np.pi - l_big_pre
                else:
                    l_big = l_big_pre

                if z_big == 0:
                    print "Z == 0!!!"
                    b_big = 0
                    print r, a
                    h_big = r - a
                else:
                    ############## ПОЛНОСТЬЮ ПЕРЕСМОТРЕТЬ ФРАГМЕНТ ###########################################
                    print 'Lets count!'
                    r_big = (x_big**2 + y_big**2 + z_big**2)**0.5
                    c = np.arcsin(z_big/r_big)
                    p = (e2*a)/(2*r_big)
                    s1 = 0
                    d = 1
                    d_norm = np.radians(0.0001/3600)
                    while d >= d_norm:
                        b = c + s1
                        s2 = np.arcsin((p*np.sin(2*b))/(Srv().bigsqrt(ellipsoid, np.degrees(b))))
                        d = abs(s2 - s1)
                        s1 = s2
                        print "!", d, "!", np.degrees(b), s2, s1, d_norm, '...'
                    b_big = b
                    print r*np.cos(b_big), "++"
                    print z_big*np.sin(b_big), "++"
                    print a*Srv().bigsqrt(ellipsoid, np.degrees(b)), "--++--"
                    h_big = r*np.cos(b_big) + z_big*np.sin(b_big) - a*Srv().bigsqrt(ellipsoid, np.degrees(b))
                    ###########################################################################################
        return [np.degrees(b_big), np.degrees(l_big), h_big]


system = "WGS84"
# b = Srv().dms2ddd("55 51 4.8857")
b = Srv().dms2ddd("54 41 12.20557")
l = Srv().dms2ddd("25 17 50.97001")
h = 191.000
print "Params", Srv().ellparams(system)
print b, l, h
print "Geocentr: ", IntTrans().geodez2geocentr(system, b, l, h)
x = IntTrans().geodez2geocentr(system, b, l, h)[0]
y = IntTrans().geodez2geocentr(system, b, l, h)[1]
z = IntTrans().geodez2geocentr(system, b, l, h)[2]
print "X = ", x
print "Y = ", y
print "Z = ", z
print IntTrans().geocentr2geodez(system, x, y, z)
