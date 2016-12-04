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



ellipsoid_name = 'WGS84'
B_dms = '55 51 6.47082'
L_dms = '37 26 11.16388'
H = 189.000
B = Srv().dms2ddd(B_dms)
L = Srv().dms2ddd(L_dms)
print Transform().blh2xyz(B, L, H, ellipsoid_name)






