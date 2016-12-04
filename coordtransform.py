# -*- coding: utf-8 -*-

import yaml as ym
import numpy as np


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
        return degrees

    def ellipsoid_params(self, ellipsoid_name):
        ell_library = Srv().unpack('ellipsoids.yml')
        elldict = {}
        if ellipsoid_name in ell_library:
            elldict['a'] = ell_library[ellipsoid_name]['a']
            elldict['ralpha'] = ell_library[ellipsoid_name]['ralpha']
        return elldict

class Transform:  # Класс для трансформирования координат

    def __init__(self):
        pass

    def blh2xyz(self, B, L, H, ellipsoid_name):  # Пересчёт эллипсоидальных координат B, L, H в геоцентрические X, Y, Z
        # Входные данные должны быть представлены в виде массива [B, L, H] - геодезические широта, долгота и высота.
        # Широта и долгота в виде ггг.гггггг (т.е. градусы и доли градуса), высота в метрах и долях метра.
        # Получаем параметры эллипсоида
        a = float(Srv().ellipsoid_params(ellipsoid_name)['a'])  # Большая полуось
        alpha = 1/float(Srv().ellipsoid_params(ellipsoid_name)['ralpha'])  # Сжатие эллипсоида
        # Синус и косинус широты и долготы
        sinB = np.sin(np.radians(B))
        sinL = np.sin(np.radians(L))
        cosB = np.cos(np.radians(B))
        cosL = np.cos(np.radians(L))
        sinsqB = sinB**2  # Квадрат синуса широты
        esq = 2*alpha - alpha**2  # Квадрат первого эксцентриситета
        subradical = 1 - esq*sinsqB  # Подкоренное значение для вычисления радиуса кривизны первого вертикала
        N = a/(subradical ** 0.5)  # Вычисление радиуса первого вертикала
        # Вычисление X, Y, Z
        X = (N + H) * cosB * cosL
        Y = (N + H) * cosB * sinL
        Z = (H + N) * sinB - esq * N * sinB
        # Выдача результата в виде массива [X, Y, Z]
        return [X, Y, Z]

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

    def ell2gauss(self, B, L, H, ellipsoid_name):
        n = int((6 + L)/6)
        l = (L - (3 + 6*(n - 1)))/57.29577951
        sinB = np.sin(np.radians(B))
        sin2B = np.sin(np.radians(2*B))
        cosB = np.cos(np.radians(B))
        sinL = np.sin(np.radians(L))
        cosL = np.cos(np.radians(L))
        x_claster1 = (l**2)*(109500 - 574700*(sinB**2) + 863700*(sinB**4) - 398600*(sinB**6))
        print 'x_claster1 = ', x_claster1
        x_claster2 = (l**2)*(278194.0 - 830174.0*(sinB**2) + 572434.0*(sinB**4) - 16010.0*(sinB**6) + x_claster1)
        print 'x_claster2 = ', x_claster2
        x_claster3 = (l**2)*(672483.4 - 811219.9*(sinB**2) + 5420.0*(sinB**4) - 10.6*(sinB**6) + x_claster2)
        print 'x_claster3 = ', x_claster3
        x_claster4 = (l**2)*(1594561.25 + 5336.535*(sinB**2) + 26.790*(sinB**4) + 0.149*(sinB**6) + x_claster3)
        x_claster5 = (l**2)*(16002.8900 + 66.9607*(sinB**2) + 0.3515*(sinB**4) - x_claster4)
        x = 6367558.4968*np.radians(B) - sin2B*x_claster5
        return n, x

ellipsoid_name = 'WGS84'
B_dms = '55 51 6.47082'
L_dms = '37 26 11.16388'
H = 189.000
B = Srv().dms2ddd(B_dms)
L = Srv().dms2ddd(L_dms)
print Transform().blh2xyz(B, L, H, ellipsoid_name)
print str(Transform().ell2gauss(B, L, H, ellipsoid_name))






