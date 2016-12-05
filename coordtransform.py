# -*- coding: utf-8 -*-

import yaml as ym
import numpy as np


class Srv:

    def __init__(self):
        pass

    def unpack(self, filename):
        with open(filename) as yml:
            return ym.load(yml)

    # Преобразование эллипсоидальных координат из градусов, минут и секунд
    # в градусы и доли градусов.
    def dms2ddd(self, dms):
        dmssplit = dms.split(' ')
        degrees = float(dmssplit[0]) + float(dmssplit[1]) / 60 + float(dmssplit[2]) / 3600
        return degrees

    # Извлечение из внешнего файла параметров нужного эллипсоида.
    def ellipsoid_params(self, ellipsoid_name):
        ell_library = Srv().unpack('ellipsoids.yml')
        elldict = {}
        if ellipsoid_name in ell_library:
            elldict['a'] = float(ell_library[ellipsoid_name]['a'])
            elldict['alpha'] = 1/float(ell_library[ellipsoid_name]['ralpha'])
        return elldict

    def transform_params(self, oldsystem, newsystem):
        # Определение семи параметров трансформирования координат.
        # По введённым системам координат получаем
        ellnsys = self.unpack('ellsandsystems.yml')
        params = self.unpack('transes.yml')
        stringtotrans = oldsystem + '--' + newsystem
        paramdict = {}
        if stringtotrans in params:
            paramdict['dx'] = eval(params[stringoftrans]['dx'])
            paramdict['dy'] = eval(params[stringoftrans]['dy'])
            paramdict['dz'] = eval(params[stringoftrans]['dz'])
            paramdict['wx'] = eval(params[stringoftrans]['wx'])
            paramdict['wy'] = eval(params[stringoftrans]['wy'])
            paramdict['wz'] = eval(params[stringoftrans]['wz'])
            paramdict['m'] = eval(params[stringoftrans]['m'])
        return paramdict



class Transform:  # Класс для трансформирования координат

    def __init__(self):
        pass

    def blh2xyz(self, B, L, H, ellipsoid_name):  # Пересчёт эллипсоидальных координат B, L, H в геоцентрические X, Y, Z
        # Входные данные должны быть представлены в виде массива [B, L, H] - геодезические широта, долгота и высота.
        # Широта и долгота в виде ггг.гггггг (т.е. градусы и доли градуса), высота в метрах и долях метра.
        # Получаем параметры эллипсоида
        a = Srv().ellipsoid_params(ellipsoid_name)['a']  # Большая полуось
        alpha = Srv().ellipsoid_params(ellipsoid_name)['alpha']  # Сжатие эллипсоида
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

    def ell2gauss(self, B, L, H, ellipsoid_name):
        print 'B, L', B, L
        n = int((6 + L)/6)
        l_part1 = 6*(n - 1)
        l_part2 = 3 + l_part1
        l_part3 = L - l_part2
        l = l_part3/57.29577951
        # Вычисление синуса и косинуса широты и долготы.
        sinB = np.sin(np.radians(B))
        sin2B = np.sin(np.radians(B*2.0))
        cosB = np.cos(np.radians(B))
        sinL = np.sin(np.radians(L))
        cosL = np.cos(np.radians(L))
        # В связи с громоздкостью формул вычисления плоских координат, вычисление будет проходить постепенно,
        # от кластера к кластеру.
        # Поэтапное вычисление x.
        x_part1 = (l**2)*(109500.0 - 574700.0*(sinB**2) + 863700.0*(sinB**4) - 398600.0*(sinB**6))
        x_part2 = (l**2)*(278194.0 - 830174.0*(sinB**2) + 572434.0*(sinB**4) - 16010.0*(sinB**6) + x_part1)
        x_part3 = (l**2)*(672483.4 - 811219.9*(sinB**2) + 5420.0*(sinB**4) - 10.6*(sinB**6) + x_part2)
        x_part4 = (l**2)*(1594561.25 + 5336.535*(sinB**2) + 26.790*(sinB**4) + 0.149*(sinB**6) + x_part3)
        x_part5 = (l**2)*(16002.8900 + 66.9607*(sinB**2) + 0.3515*(sinB**4) - x_part4)
        x = 6367558.4968*np.radians(B) - sin2B*x_part5
        # Поэтапное вычисление y.
        y_part1 = (l**2)*(79690.0 - 866190.0*(sinB**2) + 1730360.0*(sinB**4) - 945460.0*(sinB**6))
        y_part2 = (l**2)*(270806.0 - 1523417.0*(sinB**2) + 1327645.0*(sinB**4) - 21701.0*(sinB**6) + y_part1)
        y_part3 = (l**2)*(1070204.16 - 2136826.66*(sinB**2) + 17.98*(sinB**4) - 11.99*(sinB**6) + y_part2)
        y_part4 = 6378245.0 + 21346.1415*(sinB**2) + 107.1590*(sinB**4) + 0.5977*(sinB**6) + (l**2)*y_part3
        y = (5 + 10*n)*(10**5) + l*cosB*y_part4
        return n, x, y

    def recthelmert(self, X, Y, Z, old_system, new_system):
        omegamatr = np.matrix([[1, wz, -wy],
                               [-wz, 1, wx],
                               [wy, -wx, 1]])
        deltamatr = np.matrix([[dx],
                               [dy],
                               [dz]])
        oldmatr = np.matrix([[X],
                             [Y],
                             [Z]])




ellipsoid_name = 'Kras40'
B_dms = '55 45 48.85872'
L_dms = '37 39 49.14229'
H = 155.508
B = Srv().dms2ddd(B_dms)
L = Srv().dms2ddd(L_dms)
print Transform().blh2xyz(B, L, H, ellipsoid_name)
print str(Transform().ell2gauss(B, L, H, ellipsoid_name))
