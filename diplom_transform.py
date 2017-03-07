# -*- coding: utf-8 -*-

import numpy as np
import json


class Srv:

    def __init__(self):
        self.elldict = self.unpack("ellipsoids.json")
        self.transparams = self.unpack("transform.json")

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
        sinqb = (np.sin(lat))**2
        radix = (1 - e2*sinqb)**0.5
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

    def ellbysystem(self, system):
        for oneellipsoid in self.elldict:
            if system in oneellipsoid["systems"]:
                return self.ellparams(oneellipsoid["name"])
    
    # Получение параметров трансформирования. Параметры старого и нового эллипсоидов плюс семь ключей.
    def transformparams(self, oldsys, newsys):
        readydict = {}
        for onesetparam in self.transparams:
            if oldsys == onesetparam["oldsys"] and newsys == onesetparam["newsys"]:
                a_old = self.ellbysystem(oldsys)["a"]
                alpha_old = self.ellbysystem(oldsys)["alpha"]
                e2_old = self.ellbysystem(oldsys)["e2"]
                a_new = self.ellbysystem(newsys)["a"]
                alpha_new = self.ellbysystem(newsys)["alpha"]
                e2_new = self.ellbysystem(newsys)["e2"]
                dx = float(onesetparam["dx"])
                dy = float(onesetparam["dy"])
                dz = float(onesetparam["dz"])
                wx = np.radians(float(onesetparam["wx"])/3600)
                wy = np.radians(float(onesetparam["wy"])/3600)
                wz = np.radians(float(onesetparam["wz"])/3600)
                m = float(onesetparam["m"])/(10**6)
                readydict = {"a_old": a_old,
                             "alpha_old": alpha_old,
                             "e2_old": e2_old,
                             "a_new": a_new,
                             "alpha_new": alpha_new,
                             "e2_new": e2_new,
                             "dx": dx,
                             "dy": dy,
                             "dz": dz,
                             "wx": wx,
                             "wy": wy,
                             "wz": wz,
                             "m": m
                             }
        return readydict
                

# Функции для внутреннего преобразования координат
class IntTrans:

    def __init__(self):
        pass

    # Преобразование геодезических координат в геоцентрические
    def geodez2geocentr(self, ellipsoid, b_big, l_big, h_big):
        n_big = Srv().n_big(ellipsoid, np.radians(b_big))
        e2 = Srv().ellparams(ellipsoid)["e2"]
        sinb = np.sin(np.radians(b_big))
        cosb = np.cos(np.radians(b_big))
        sinl = np.sin(np.radians(l_big))
        cosl = np.cos(np.radians(l_big))
        x_big = (n_big + h_big)*cosb*cosl
        y_big = (n_big + h_big)*cosb*sinl
        z_big = sinb*(n_big*(1 - e2) + h_big)
        return [x_big, y_big, z_big]

    # Итеративный способ преобразования геоцентрических координат в геодезические
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
                if x_big > 0 > y_big:
                    l_big = 2*np.pi - l_big_pre
                elif x_big < 0 and y_big < 0:
                    l_big = np.pi + l_big_pre
                elif y_big > 0 > x_big:
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

    # Метод преобразования геодезических координат в координаты Гаусса, немного отличающийся
    # от описанного в работе.
    def geodez2gauss(self, ellipsoid, b_big, l_big, h_big):
        a = Srv().ellparams(ellipsoid)["a"]
        e2 = Srv().ellparams(ellipsoid)["e2"]
        e4 = e2**2
        e6 = e2**3
        e8 = e2**4
        es2 = e2 / (1 - e2)
        ro = (360*3600)/(2.0*np.pi)
        b = np.radians(b_big)
        n_big = Srv().n_big(ellipsoid, b)
        nzone = int((l_big + 6) / 6.0)
        l = (l_big - (3.0 + 6.0*(nzone - 1))) / 57.29577951
        sinb = np.sin(b)
        sin2b = np.sin(2.0*b)
        sin4b = np.sin(4.0*b)
        sin6b = np.sin(6.0*b)
        cosb = np.cos(b)
        cosb2 = cosb**2
        cosb4 = cosb**4
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
        x_part1 = (n_big*cosb*sinb*(l**2))/2.0
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
        y = y_pre + 500000.0
        return [nzone, x, y, h_big]

    # Преобразование координат Гаусса-Крюгера в геодезические по методу, предложенному Серапинасом
    def gauss2geodez_sepns(self, ellipsoid, nzone, x, y):
        ys = y - 500000
        a = Srv().ellparams(ellipsoid)["a"]
        e2 = Srv().ellparams(ellipsoid)["e2"]
        es2 = e2 / (1 - e2)
        s0 = a*(1 - e2)
        s2 = (3.0/2.0)*e2*s0
        s4 = (5.0/4.0)*e2*s2
        s6 = (7.0/6.0)*e2*s4
        s8 = (9.0/8.0)*e2*s6
        c_big_0 = s0 + (1.0/2.0)*s2 + (3.0/8.0)*s4 + (5.0/16.0)*s6 + (35.0/128.0)*s8
        c_big_2 = (1.0/4.0)*(s2 + s4 + (15.0/16.0)*s6 + (7.0/8.0)*s8)
        c_big_4 = (1.0/32.0)*(s4 + (3.0/2.0)*s6 + (7.0/4.0)*s8)
        c_big_6 = (1.0/96.0)*((1.0/2.0)*s6 + s8)
        d_big_2 = (c_big_2/c_big_0)*(1 + c_big_4/c_big_0 - (c_big_2**2)/(2*(c_big_0**2)))
        d_big_4 = (c_big_2**2)/(c_big_0**2) - (c_big_4/c_big_0)
        d_big_6 = c_big_6/c_big_0 - ((3*c_big_2)/c_big_0)*(c_big_4/c_big_0 - ((c_big_2**2)/(2*(c_big_0**2))))
        beta = x/c_big_0
        sin2beta = np.sin(2*beta)
        sin4beta = np.sin(4*beta)
        sin6beta = np.sin(6*beta)
        b_big_x = beta + d_big_2*sin2beta + d_big_4*sin4beta + d_big_6*sin6beta
        n_big_x = Srv().n_big(ellipsoid, b_big_x)
        m_big_x = Srv().m_big(ellipsoid, b_big_x)
        tanx = np.tan(b_big_x)
        tanx2 = tanx**2
        tanx4 = tanx**4
        etax2 = es2*(np.cos(b_big_x)**2)
        b_big_part1 = (tanx*ys**2)/(2*m_big_x*n_big_x)  # Выражение перед большой и страшной скобкой
        b_big_part21 = (ys**2)/(12*(n_big_x**2))
        b_big_part22 = 5 + 3*tanx2 + etax2 - 9*etax2*tanx2
        b_big_part31 = (ys**4)/(360*(n_big_x**4))
        b_big_part32 = 61 + 90*tanx2 + 45*tanx4
        b_big_part2 = b_big_part21*b_big_part22
        b_big_part3 = b_big_part31*b_big_part32
        b_big = b_big_x - b_big_part1*(1 - b_big_part2 + b_big_part3)
        l_part1 = ys/(n_big_x*np.cos(b_big_x))
        l_part21 = (ys**2)/(6*(n_big_x**2))
        l_part22 = 1 + 2*tanx2 + etax2
        l_part31 = (ys**4)/(120*(n_big_x**4))
        l_part32 = 5 + 28*tanx2 + 24*tanx4 + 6*etax2 + 8*etax2*tanx2
        l_part2 = l_part21*l_part22
        l_part3 = l_part31*l_part32
        l = l_part1*(1 - l_part2 + l_part3)
        l_big_0 = np.radians(6*nzone - 3)
        l_big = l_big_0 + l
        return [np.degrees(b_big), np.degrees(l_big)]


# Внешние преобразования координат
class ExtTrans:

    def __init__(self):
        pass

    # Способ Гельмерта для преобразования геоцентрических координат
    def geocentrhelmert(self, x_big, y_big, z_big, oldsys, newsys):
        transpars = Srv().transformparams(oldsys, newsys)
        dx = transpars["dx"]
        dy = transpars["dy"]
        dz = transpars["dz"]
        wx = transpars["wx"]
        wy = transpars["wy"]
        wz = transpars["wz"]
        m = transpars["m"]
        print dx, dy, dz, wx, wy, wz, m
        oldmatrix = np.matrix([[x_big],
                               [y_big],
                               [z_big]])
        deltamatrix = np.matrix([[dx],
                                 [dy],
                                 [dz]])
        omegamatrix = np.matrix([[1, wz, -wy],
                                 [-wz, 1, wx],
                                 [wy, -wx, 1]])
        newmatrix = deltamatrix + (1 + m)*omegamatrix*oldmatrix
        result = np.squeeze(np.asarray(newmatrix))
        return result

    # Способ преобразования геодезических координат, предложенный в ГОСТ
    def geodezgost(self, b_big, l_big, h_big, oldsys, newsys):
        transpars = Srv().transformparams(oldsys, newsys)
        dx = transpars["dx"]
        dy = transpars["dy"]
        dz = transpars["dz"]
        wx = transpars["wx"]
        wy = transpars["wy"]
        wz = transpars["wz"]
        m = transpars["m"]
        a_old = transpars["a_old"]
        e2_old = transpars["e2_old"]
        a_new = transpars["a_new"]
        e2_new = transpars["e2_new"]
        sinb = np.sin(np.radians(b_big))
        sinb2 = sinb**2
        sinl = np.sin(np.radians(l_big))
        cosb = np.cos(np.radians(b_big))
        cos2b = np.cos(np.radians(2*b_big))
        cosl = np.cos(np.radians(l_big))
        tgb = np.tan(np.radians(b_big))
        delta_a = a_new - a_old
        delta_e2 = e2_new - e2_old
        a = (a_old + a_new)/2
        e2 = (e2_old + e2_new)/2
        ro = (360 * 3600) / (2.0 * np.pi)
        n_big = a/((1 - e2*(sinb**2))**0.5)
        m_big = (a*(1 - e2))/((1 - e2*(sinb**2))**(3.0/2.0))
        delta_b_part1 = ro/(n_big + m_big)
        delta_b_part2 = (n_big/a)*e2*sinb*cosb*delta_a
        delta_b_part3 = ((n_big**2)/(a**2) + 1)*n_big*sinb*cosb*((delta_e2**2)/2)
        delta_b_part4 = sinb*(dx*cosl + dy*sinl)
        delta_b_part5 = dz*cosb
        delta_b_part6 = wx*sinl*(1 - e2*cos2b)
        delta_b_part7 = wy*cosl*(1 + e2*cos2b)
        delta_b_part8 = ro*m*e2*sinb*cosb
        delta_b_sec = delta_b_part1*(delta_b_part2 + delta_b_part3 - delta_b_part4 + delta_b_part5) - \
                  delta_b_part6 + delta_b_part7 - delta_b_part8
        delta_l_part1 = (ro/((n_big + h_big)*cosb))*(- dx*sinl + dy*cosl)
        delta_l_part2 = tgb*(1 - e2)*(wx*cosl + wy*sinl)
        delta_l_sec = delta_l_part1 + delta_l_part2 - wz
        delta_h_part1 = -(a/n_big)*delta_a
        delta_h_part2 = n_big*sinb2*(delta_e2/2)
        delta_h_part3 = cosb*(dx*cosl + dy*sinl)
        delta_h_part4 = dz*sinb
        delta_h_part5 = n_big*e2*sinb*cosb*((wx/ro)*sinl - (wy/ro)*cosl)
        delta_h_part6 = m*(((a**2)/n_big) + h_big)
        delta_h = delta_h_part1 + delta_h_part2 + delta_h_part3 + delta_h_part4 - delta_h_part5 + delta_h_part6
        delta_b = np.radians(delta_b_sec/3600.0)
        delta_l = np.radians(delta_l_sec/3600.0)
        new_b_big = np.radians(b_big) + delta_b
        new_l_big = np.radians(l_big) + delta_l
        new_h_big = h_big + delta_h
        return [np.degrees(new_b_big),
                np.degrees(new_l_big),
                new_h_big]

    # Стандартный способ Молоденского для преобразования геодезических координат
    def geosezmolodstandard(self, b_big, l_big, h_big, oldsys, newsys):
        transpars = Srv().transformparams(oldsys, newsys)
        dx = transpars["dx"]
        dy = transpars["dy"]
        dz = transpars["dz"]
        a_old = transpars["a_old"]
        e2_old = transpars["e2_old"]
        alpha_old = transpars["alpha_old"]
        a_new = transpars["a_new"]
        e2_new = transpars["e2_new"]
        alpha_new = transpars["alpha_new"]
        delta_a = a_new - a_old
        a = (a_old + a_new) / 2
        e2 = (e2_old + e2_new) / 2
        alpha = (alpha_old + alpha_new)/2
        delta_alpha = alpha_new - alpha_old
        sinb = np.sin(np.radians(b_big))
        sinb2 = sinb ** 2
        sinl = np.sin(np.radians(l_big))
        cosb = np.cos(np.radians(b_big))
        cosl = np.cos(np.radians(l_big))
        n_big = Srv().n_big(Srv().ellipsoidbysystem(oldsys), np.radians(b_big))
        m_big = Srv().m_big(Srv().ellipsoidbysystem(oldsys), np.radians(b_big))
        delta_b_part1 = 1.0/(m_big + h_big)
        delta_b_part2 = -dx*sinb*cosl
        delta_b_part3 = dy*sinb*sinl
        delta_b_part4 = dz*cosb
        delta_b_part5 = ((n_big*e2*sinb*cosb)/a)*delta_a
        delta_b_part6 = sinb*cosb*((m_big/(1 - alpha)) + n_big*(1 - alpha))*delta_alpha
        delta_b = delta_b_part1*(delta_b_part2 - delta_b_part3 + delta_b_part4 + delta_b_part5 + delta_b_part6)
        delta_l_part1 = 1.0/((n_big + h_big)*cosb)
        delta_l_part2 = -dx*sinl
        delta_l_part3 = dy*cosl
        delta_l = delta_l_part1*(delta_l_part2 + delta_l_part3)
        delta_h_part1 = dx*cosb*cosl
        delta_h_part2 = dy*cosb*sinl
        delta_h_part3 = dz*sinb
        delta_h_part4 = (a/n_big)*delta_a
        delta_h_part5 = n_big*(1 - alpha)*sinb2*delta_alpha
        delta_h = delta_h_part1 + delta_h_part2 + delta_h_part3 - delta_h_part4 + delta_h_part5
        b_big_new = b_big + np.degrees(delta_b)
        l_big_new = l_big + np.degrees(delta_l)
        h_big_new = h_big + delta_h
        return [b_big_new, l_big_new, h_big_new]

    # Сокращённый способ Молоденского для преобразования геодезических координат
    def geodezmolodshort(self, b_big, l_big, h_big, oldsys, newsys):
        transpars = Srv().transformparams(oldsys, newsys)
        dx = transpars["dx"]
        dy = transpars["dy"]
        dz = transpars["dz"]
        a_old = transpars["a_old"]
        e2_old = transpars["e2_old"]
        alpha_old = transpars["alpha_old"]
        a_new = transpars["a_new"]
        e2_new = transpars["e2_new"]
        alpha_new = transpars["alpha_new"]
        delta_a = a_new - a_old
        a = (a_old + a_new) / 2
        e2 = (e2_old + e2_new) / 2
        alpha = (alpha_old + alpha_new) / 2
        delta_alpha = alpha_new - alpha_old
        sinb = np.sin(np.radians(b_big))
        sin2b = np.sin(np.radians(2*b_big))
        sinb2 = sinb ** 2
        sinl = np.sin(np.radians(l_big))
        cosb = np.cos(np.radians(b_big))
        cosl = np.cos(np.radians(l_big))
        n_big = Srv().n_big(Srv().ellipsoidbysystem(oldsys), np.radians(b_big))
        m_big = Srv().m_big(Srv().ellipsoidbysystem(oldsys), np.radians(b_big))
        delta_b_part1 = 1/m_big
        delta_b_part2 = -dx*sinb*cosl
        delta_b_part3 = dy*sinb*sinl
        delta_b_part4 = dz*cosb
        delta_b_part5 = sin2b*(alpha*delta_a + a*delta_alpha)
        delta_b = delta_b_part1*(delta_b_part2 - delta_b_part3 + delta_b_part4 + delta_b_part5)
        delta_l_part1 = 1/(n_big*cosb)
        delta_l_part2 = -dx*sinl
        delta_l_part3 = dy*cosl
        delta_l = delta_l_part1*(delta_l_part2 + delta_l_part3)
        delta_h_part1 = dx*cosb*cosl
        delta_h_part2 = dy*cosb*sinl
        delta_h_part3 = dz*sinb
        delta_h_part4 = delta_a
        delta_h_part5 = sinb2*(alpha*delta_a + a*delta_alpha)
        delta_h = delta_h_part1 + delta_h_part2 + delta_h_part3 - delta_h_part4 + delta_h_part5
        b_big_new, l_big_new, h_big_new = b_big + np.degrees(delta_b), l_big + np.degrees(delta_l), h_big + delta_h
        return [b_big_new, l_big_new, h_big_new]


system = "WGS84"
b = Srv().dms2ddd("55 50 54.95392")
l = Srv().dms2ddd("37 27 15.10279")
h = 193
print "Params", Srv().ellparams(system)
print b, l, h
print "Geocentr: ", IntTrans().geodez2geocentr(system, b, l, h)
x_big = IntTrans().geodez2geocentr(system, b, l, h)[0]
y_big = IntTrans().geodez2geocentr(system, b, l, h)[1]
z_big = IntTrans().geodez2geocentr(system, b, l, h)[2]
print "X = ", x_big
print "Y = ", y_big
print "Z = ", z_big
print IntTrans().geocentr2geodez(system, x_big, y_big, z_big)
print "*******"
print IntTrans().geodez2gauss(system, b, l, h), "."
nzone = IntTrans().geodez2gauss(system, b, l, h)[0]
x = IntTrans().geodez2gauss(system, b, l, h)[1]
y = IntTrans().geodez2gauss(system, b, l, h)[2]
h = IntTrans().geodez2gauss(system, b, l, h)[3]
print "*******"
print nzone
print x
print y
print h
print "*******"
print IntTrans().gauss2geodez_sepns(system, nzone, x, y)
print "*******"
print Srv().ddd2dms(IntTrans().geocentr2geodez(system, x_big, y_big, z_big)[0])
print Srv().ddd2dms(IntTrans().gauss2geodez_sepns(system, nzone, x, y)[0])
print "*******"
print ExtTrans().geocentrhelmert(x_big, y_big, z_big, "WGS84", "PZ90")
print "*******"
print ExtTrans().geodezgost(b, l, h, "WGS84", "PZ90")
print "*******"
print ExtTrans().geosezmolodstandard(b, l, h, "WGS84", "PZ90")
print "*******"
print ExtTrans().geodezmolodshort(b, l, h, "WGS84", "PZ90")