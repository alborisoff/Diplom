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


    def bigsqrt(self, ellipsoid, lat):
        e2 = self.ellparams(ellipsoid)["e2"]
        sinqB = (np.sin(lat))**2
        radix = (1 - e2*sinqB)**0.5
        return radix

    def n_big(self, ellipsoid, lat):
        a = self.ellparams(ellipsoid)["a"]
        dnmtr = self.bigsqrt(ellipsoid, lat)
        return a/dnmtr

    def m_big(self, ellipsoid, lat):
        a = self.ellparams(ellipsoid)["a"]
        e2 = self.ellparams(ellipsoid)["e2"]
        numerator = a*(1 - e2)
        denominator = (self.bigsqrt(ellipsoid, lat))**3
        return numerator/denominator


    # Преобразование эллипсоидальных координат из градусов, минут и секунд
    # в градусы и доли градусов.
    def dms2ddd(self, dms):
        dmssplit = dms.split(" ")
        print "dmssplit", dmssplit[0], dmssplit[1], dmssplit[2]
        degrees = float(dmssplit[0]) + float(dmssplit[1])/60.0 + float(dmssplit[2])/3600.0
        return degrees

    def ddd2dms(self, ddd):
        degrees = int(ddd)
        rest_degrees = ddd - degrees
        minutes = int(rest_degrees*60)
        rest_minutes = rest_degrees*60 - minutes
        seconds = rest_minutes*60
        return [degrees, minutes, seconds]



# Функции для внутреннего преобразования координат
class IntTrans:

    # Преобразование геодезических координат в прямоугольные
    def geodez2geocentr(self, ellipsoid, b_big, l_big, h_big):
        n_big = Srv().n_big(ellipsoid, np.radians(b_big))
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
                h_big = z_big*np.sin(b_big) - a*Srv().bigsqrt(ellipsoid, b_big)
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
                    r_big = (x_big**2 + y_big**2 + z_big**2)**0.5
                    c = np.arcsin(z_big/r_big)
                    p = (e2*a)/(2.0*r_big)
                    s1 = 0
                    s2 = 0
                    d_norm = np.radians(0.0001/3600)
                    d = 1.0
                    while d >= d_norm:
                        b = c + s1
                        s2 = np.arcsin(p*np.sin(2*b)/Srv().bigsqrt(ellipsoid, b))
                        d = abs(s2 - s1)
                        s1 = s2
                    b_big = b
                    h_big = r*np.cos(b_big) + z_big*np.sin(b_big) - a*Srv().bigsqrt(ellipsoid, b_big)
        return [np.degrees(b_big), np.degrees(l_big), h_big]

    # Преобразование геодезических координат в плоские прямоугольные Гаусса
    def geodez2gaussb(self, ellipsoid, b_big, l_big, h_big):
        a = Srv().ellparams(ellipsoid)["a"]
        e2 = Srv().ellparams(ellipsoid)["e2"]
        e4 = e2**2
        es2 = e2/(1 - e2)
        ro = (360*3600)/(2.0*np.pi)
        b = np.radians(b_big)
        n_big = Srv().n_big(ellipsoid, b)
        m_big = Srv().m_big(ellipsoid, b)
        l = np.radians(l_big)
        sinb = np.sin(b)
        sinb2 = sinb**2
        sinb4 = sinb**4
        sin2b = np.sin(2.0*b)
        cosb = np.cos(b)
        cosb3 = cosb**3
        cosb5 = cosb**5
        cosb7 = cosb**7
        cos2b = np.cos(2.0*b)
        tanb = np.tan(b)
        tanb2 = tanb**2
        tanb4 = tanb**4
        tanb6 = tanb**6
        eta = (es2**0.5)*cosb
        bsec = b_big*3600
        nzone = int((l_big + 6)/6.0)
        l = (l_big - (3.0 + 6.0*(nzone - 1)))/57.29577951
        print "l = ", l
        a_big_x = 1 + (3.0/4.0)*e2 + (45.0/64.0)*e4
        b_big_x = (3.0/4.0)*e2 + (15.0/16.0)*e4
        c_big_x = (15.0/64.0)*e4
        x_big = a*(1.0 - e2)*(a_big_x*(bsec/ro) - (b_big_x/2.0)*sin2b + (c_big_x/4.0)*sinb4)
        a2 = (1.0/2.0)*n_big*sinb*cosb
        a4 = (1.0/24.0)*n_big*sinb*cosb3*(5.0 - tanb2 + 9.0*(eta**2) + 4.0*(eta**4))
        a6 = (1.0/720.0)*n_big*sinb*cosb5*(61.0 - 58.0*tanb2 + tanb4 +
                                           270.0*(eta**2) - 330.0*(eta**2)*tanb2)
        a8 = (1.0/40320.0)*n_big*sinb*cosb7*(1385.0 - 3111.0*tanb2 + 543.0*tanb4 - tanb6)
        b1 = n_big*cosb
        b3 = (1.0/6.0)*n_big*cosb3*(-tanb2 + eta**2)
        b5 = (1.0/120.0)*n_big*cosb5*(5.0 - 18.0*tanb2 + tanb4 - 14.0*(eta**2) - 58.0*(eta**2)*tanb2)
        b7 = (1.0/5040.0)*n_big*cosb7*(61.0 - 479.0*tanb2 + 179.0*tanb4 - tanb6)
        x = x_big + a2*(l**2) + a4*(l**4) + a6*(l**6) + a8*(l**8)
        print "x_big = ", x_big
        print a2*(l**2)
        print a4*(l**4)
        print a6*(l**6)
        print a8*(l**8)
        y = b1*l + b3*(l**3) + b5*(l**5) + b7*(l**7) + 500000
        return [nzone, x, y]

    # Метод преобразования геодезических координат в координаты Гаусса, немного отличающийся
    # от описанного в работе.
    def geodez2gauss(self, ellipsoid, b_big, l_big):
        a = Srv().ellparams(ellipsoid)["a"]
        e2 = Srv().ellparams(ellipsoid)["e2"]
        e4 = e2**2
        e6 = e2**3
        e8 = e2**4
        es2 = e2 / (1 - e2)
        ro = (360*3600)/(2.0*np.pi)
        b = np.radians(b_big)
        n_big = Srv().n_big(ellipsoid, b)
        m_big = Srv().m_big(ellipsoid, b)
        nzone = int((l_big + 6) / 6.0)
        l = (l_big - (3.0 + 6.0*(nzone - 1))) / 57.29577951
        sinb = np.sin(b)
        sin2b = np.sin(2.0*b)
        sin4b = np.sin(4.0*b)
        sin6b = np.sin(6.0*b)
        cosb = np.cos(b)
        cosb2 = cosb**2
        cosb3 = cosb**3
        cosb4 = cosb**4
        cosb5 = cosb**5
        cosb7 = cosb**7
        cos2b = np.cos(2.0*b)
        tanb = np.tan(b)
        tanb2 = tanb**2
        tanb4 = tanb**4
        eta = (es2**0.5)*cosb
        bsec = b_big*3600
        g0 = 1 + (3.0/4.0)*e2 + (45.0/64.0)*e4 + (175.0/256.0)*e6 + (11025.0/16384.0)*e8
        g1 = - (3.0/8.0)*e2 - (15.0/32.0)*e4 - (525.0/1024.0)*e6 - (2205.0/4096.0)*e8
        g2 = (15.0/256.0)*e4 + (105.0/1024.0)*e6 + (2205.0/16384.0)*e8
        g3 = - (35.0/3072.0)*e6 - (315.0/12288.0)*e8
        x_big = a*(1 - e2)*(g0*(bsec/ro) + g1*sin2b + g2*sin4b + g3*sin6b)
        x_part1 = (n_big*cosb*sinb*(l**2))/(2.0)
        x_part21 = cosb2*(l**2)/12.0
        x_part22 = 5.0 - tanb2 + 9.0*(eta**2)
        x_part31 = ((l**4)*cosb4)/360.0
        x_part32 = 61.0 - 58.0*tanb2 + tanb4
        x = x_big + x_part1*(1.0 + x_part21*x_part22 + x_part31*x_part32)
        y_part1 = n_big*cosb*l
        y_part21 = ((l**2)*cosb2)/6.0
        y_part22 = 1.0 - tanb2 + eta**2
        y_part31 = ((l**4)*cosb4)/120.0
        y_part32 = 5.0 - 18.0*tanb2 + tanb4 + 14*(eta**2) - 58*(eta**2)*tanb2
        y_pre = y_part1*(1.0 + y_part21*y_part22 + y_part31*y_part32)
        y = y_pre + nzone*1000000.0 + 500000.0
        return [x, y]


class ExtTrans:

    def __init__(self):
        pass

    def geocentrhelmert(self, x_big, y_big, z_big):
        oldmatrix = np.matrix([[x_big],
                               [y_big],
                               [z_big]])






system = "Krasovsky1940"
b = Srv().dms2ddd("55 50 54.80512")
l = Srv().dms2ddd("37 27 21.87885")
h = 188.232
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
print "*******"
print Srv().ddd2dms(b)
print Srv().ddd2dms(IntTrans().geocentr2geodez(system, x, y, z)[0])
print "*******"
print Srv().ddd2dms(l)
print Srv().ddd2dms(IntTrans().geocentr2geodez(system, x, y, z)[1])
print "*******"
print h
print IntTrans().geocentr2geodez(system, x, y, z)[2]
print "*******"
print IntTrans().geodez2gauss(system, b, l), "."
print "*******"
print IntTrans().geodez2gaussb(system, b, l, h)