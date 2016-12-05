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
        l_claster1 = 6*(n - 1)
        l_claster2 = 3 + l_claster1
        l_claster3 = L - l_claster2
        l = l_claster3/57.29577951
        # Вычисление синуса и косинуса широты и долготы.
        sinB = np.sin(np.radians(B))
        sin2B = np.sin(np.radians(B*2.0))
        cosB = np.cos(np.radians(B))
        sinL = np.sin(np.radians(L))
        cosL = np.cos(np.radians(L))
        # В связи с громоздкостью формул вычисления плоских координат, вычисление будет проходить постепенно,
        # от кластера к кластеру.
        # Поэтапное вычисление x.
        x_claster1 = (l**2)*(109500.0 - 574700.0*(sinB**2) + 863700.0*(sinB**4) - 398600.0*(sinB**6))
        x_claster2 = (l**2)*(278194.0 - 830174.0*(sinB**2) + 572434.0*(sinB**4) - 16010.0*(sinB**6) + x_claster1)
        x_claster3 = (l**2)*(672483.4 - 811219.9*(sinB**2) + 5420.0*(sinB**4) - 10.6*(sinB**6) + x_claster2)
        x_claster4 = (l**2)*(1594561.25 + 5336.535*(sinB**2) + 26.790*(sinB**4) + 0.149*(sinB**6) + x_claster3)
        x_claster5 = (l**2)*(16002.8900 + 66.9607*(sinB**2) + 0.3515*(sinB**4) - x_claster4)
        x = 6367558.4968*np.radians(B) - sin2B*x_claster5
        # Поэтапное вычисление y.
        y_claster1 = (l**2)*(79690.0 - 866190.0*(sinB**2) + 1730360.0*(sinB**4) - 945460.0*(sinB**6))
        y_claster2 = (l**2)*(270806.0 - 1523417.0*(sinB**2) + 1327645.0*(sinB**4) - 21701.0*(sinB**6) + y_claster1)
        y_claster3 = (l**2)*(1070204.16 - 2136826.66*(sinB**2) + 17.98*(sinB**4) - 11.99*(sinB**6) + y_claster2)
        y_claster4 = 6378245.0 + 21346.1415*(sinB**2) + 107.1590*(sinB**4) + 0.5977*(sinB**6) + (l**2)*y_claster3
        y = (5 + 10*n)*(10**5) + l*cosB*y_claster4
        return n, x, y

ellipsoid_name = 'Kras40'
B_dms = '55 45 48.85872'
L_dms = '37 39 49.14229'
H = 155.508
B = Srv().dms2ddd(B_dms)
L = Srv().dms2ddd(L_dms)
print Transform().blh2xyz(B, L, H, ellipsoid_name)
print str(Transform().ell2gauss(B, L, H, ellipsoid_name))






