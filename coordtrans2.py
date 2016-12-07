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

    def ellandsystem(self, system_name):
        ellnsys = self.unpack('ellsandsystems.yml')
        ellipsoid_name = ''
        for ell in ellnsys:
            if system_name == ellnsys[ell]:
                ellipsoid_name = ell
        print 'ellipsoid_name', ellipsoid_name
        return ellipsoid_name

    def transform_params(self, oldsystem, newsystem):
        # Определение семи параметров трансформирования координат.
        # По введённым системам координат получаем
        params = self.unpack('transes.yml')
        print 'params', params
        params_dict = {}
        # По входным названиям систем координат определяем эллипсоиды, на которых эти системы основаны
        old_ellipsoid = self.ellipsoid_params(self.ellandsystem(oldsystem))
        new_ellipsoid = self.ellipsoid_params(self.ellandsystem(newsystem))
        print old_ellipsoid, new_ellipsoid
        params_dict['old_a'] = old_ellipsoid['a']
        params_dict['old_alpha'] = old_ellipsoid['alpha']
        params_dict['new_a'] = new_ellipsoid['a']
        params_dict['new_alpha'] = new_ellipsoid['alpha']
        name_of_transform = oldsystem + '--' + newsystem
        if name_of_transform in params:
            params_dict['dx'] = params[name_of_transform]['dx']
            params_dict['dy'] = params[name_of_transform]['dy']
            params_dict['dz'] = params[name_of_transform]['dz']
            params_dict['wx'] = params[name_of_transform]['wx']
            params_dict['wy'] = params[name_of_transform]['wy']
            params_dict['wz'] = params[name_of_transform]['wz']
            params_dict['m'] = params[name_of_transform]['m']
        return params_dict

    
class IntTrans:  # Класс для трансформирвания координат в пределах одной системы

    def __init__(self):
        pass

    # Пересчёт эллипсоидальных координат B, L, H в геоцентрические X, Y, Z
    def blh2xyz(self, B, L, H, ellipsoid_name):
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
        sinB2 = sinB ** 2  # Квадрат синуса широты
        e2 = 2 * alpha - alpha ** 2  # Квадрат первого эксцентриситета
        subradical = 1 - e2 * sinB2  # Подкоренное значение для вычисления радиуса кривизны первого вертикала
        N = a / (subradical ** 0.5)  # Вычисление радиуса первого вертикала
        # Вычисление X, Y, Z
        X = (N + H) * cosB * cosL
        Y = (N + H) * cosB * sinL
        Z = (H + N) * sinB - e2 * N * sinB
        # Выдача результата в виде массива [X, Y, Z]
        return [X, Y, Z]

    def ell2gauss(self, B, L, H, ellipsoid_name):
        print 'B, L', B, L
        Bsec = B*3600
        n = int((6+L)/6)
        l_part1 = 6*(n - 1)
        l_part2 = 3 + l_part1
        l_part3 = L - l_part2
        l = l_part3/57.29577951
        a = Srv().ellipsoid_params(ellipsoid_name)['a']
        alpha = Srv().ellipsoid_params(ellipsoid_name)['alpha']
        # Вычисление синуса и косинуса широты и долготы.
        sinB = np.sin(np.radians(B))
        sin2B = np.sin(np.radians(2.0*B))
        sinL = np.sin(np.radians(L))
        cosB = np.cos(np.radians(B))
        cosL = np.cos(np.radians(L))
        tgB = np.tan(np.radians(B))
        tgB2 = tgB**2
        sinB2 = sinB**2  # Квадрат синуса широты
        e2 = 2*alpha - alpha**2  # Квадрат первого эксцентриситета
        es2 = e2/(1 - e2)
        subradical = 1 - e2*sinB2  # Подкоренное значение для вычисления радиуса кривизны первого вертикала
        N = a/(subradical**0.5)  # Вычисление радиуса первого вертикала
        M = N*(1 - e2)
        # Постепенный пересчёт эллиптических координат в плоские прямоугольные Гаусса.
        a2 = 0.5*(N*sinB*cosB)
        teta = cosB*(es2**0.5)
        a4 = (1/24)*(N*sinB*(cosB**3))*(5 - (tgB**2) + 9*(teta**2) + 4*(teta**4))
        a6 = (1/720)*(N*sinB*(cosB**5))*(61 - 58*(tgB**2) + tgB**4 + 270*(teta**2) - 330*(teta**2)*(tgB**2))
        a8 = (1/40320)*(N*sinB*(cosB**7))*N*sinB*(cosB**7)*(1385 -  3111*(tgB**2) + 543*(tgB**4) - (tgB**6))
        ro = 206264.806
        b1 = N*cosB
        b3 = (1/6)*N*(cosB**3)*(-(tgB**2) + teta**2)
        b5 = (1/120)*N*(cosB**5)*(5 - 18*(tgB**2) + tgB**4 - 14*(teta**2) - 58*(teta**3)*(tgB**2))
        b7 = (1/5040)*N*(cosB**7)*(61 - 479*(tgB**2) + 179*(tgB**4) - tgB**6)
        Ax = 1 + (3/4)*e2 + (45/64)*(e2**2)
        Bx = (3/4)*e2 + (15/16)*(e2**2)
        Cx = (15/64)*(e2**2)
        X = a*(1 - e2)*(Ax*(Bsec/ro) - (Bx/2)*sin2B + (Cx/4)*(sinB**4))
        x = X + a2*(l**2) + a4*(l**4) + a6*(l**6) + a8*(l**8)
        y = b1*l + b3*(l**3) + b5*(l**5) + b7*(l**7)
        return n, x, y

    def xyz2bl(self, X, Y, Z, ellipsoid_name):
        # С вычислением долготы всё просто...
        ellparams = Srv().ellipsoid_params(ellipsoid_name)
        a = ellparams['a']
        alpha = ellparams['alpha']
        e2 = 2 * alpha - alpha ** 2
        L = np.degrees(np.arctan(Y/X))
        # А вот широту придётся вычислять постепенно.
        R = (X**2 + Y**2)**0.5
        r = (Z**2 + (X**2 + Y**2)*(1 - e2))**0.5
        es2 = e2/(1 - e2)
        b = a*((1 - e2)**0.5)
        big_chisl = r**3 + b*es2*(Z**2)
        big_znaml = r**3 - b*(e2**0.5)*(1 - e2)*(R**2)
        tanB = (Z/R)*(big_chisl/big_znaml)
        B = np.degrees(np.arctan(tanB))
        return [B, L]






class ExtTrans:  # Функции трансформирования координат из одной системы в другую.

    def __init__(self):
        pass

    def recthelmert(self, X, Y, Z, old_system, new_system):  # Пересчёт прямоугольных координат по формуле Гельмерта
        # Загрузка семи параметров пересчёта
        tparameters = Srv().transform_params(old_system, new_system)
        dx = eval(tparameters['dx'])
        dy = eval(tparameters['dy'])
        dz = eval(tparameters['dz'])
        wx = eval(tparameters['wx'])
        wy = eval(tparameters['wy'])
        wz = eval(tparameters['wz'])
        m = eval(tparameters['m'])
        # Подготовка матриц пересчёта
        oldmatrix = np.matrix([[X],
                               [Y],
                               [Z]])
        deltamatrix = np.matrix([[dx],
                                 [dy],
                                 [dz]])
        omegamatrix = np.matrix([[1, wz, -wy],
                                 [-wz, 1, wx],
                                 [wy, -wx, 1]])
        newmatrix = deltamatrix + (1 + m)*omegamatrix*oldmatrix
        result = np.squeeze(np.asarray(newmatrix))
        return [result[0], result[1], result[2]]

    def molodensky(self, B, L, H, old_system, new_system):
        tparameters = Srv().transform_params(old_system, new_system)
        dx = eval(tparameters['dx'])
        dy = eval(tparameters['dy'])
        dz = eval(tparameters['dz'])
        wx = eval(tparameters['wx'])
        wy = eval(tparameters['wy'])
        wz = eval(tparameters['wz'])
        m = eval(tparameters['m'])
        old_a = eval(tparameters['old_a'])
        new_a = eval(tparameters['new_a'])
        old_alpha = eval(tparameters['old_alpha'])
        new_alpha = eval(tparameters['new_alpha'])
        da = new_a - old_a
        old_e2 = 2 * old_alpha - old_alpha ** 2  # Квадрат первого эксцентриситета
        es2 = old_e2 / (1 - old_e2)
        subradical = 1 - old_e2 * sinB2  # Подкоренное значение для вычисления радиуса кривизны первого вертикала
        old_N = a / (subradical ** 0.5)  # Вычисление радиуса первого вертикала
        old_M = old_N * (1 - old_e2)
        ro = 206264.806
        # Вычисление синуса и косинуса широты и долготы.
        sinB = np.sin(np.radians(B))
        sin2B = np.sin(np.radians(2.0 * B))
        sinL = np.sin(np.radians(L))
        cosB = np.cos(np.radians(B))
        cosL = np.cos(np.radians(L))
        tgB = np.tan(np.radians(B))
        tgB2 = tgB ** 2
        sinB2 = sinB ** 2  # Квадрат синуса широты
        # Постепенное вычисление поправок широты, долготы и высоты.
        dB_part1 = ro/(old_M + H)
        dB_part2 = (N/old_a)










print Srv().transform_params('WGS84', 'SK42')

ellipsoid_name = 'WGS84'
prelat = '56 59 9.40754'
prelon = '41 0 12.94567'
H = 146
lat = Srv().dms2ddd(prelat)
lon = Srv().dms2ddd(prelon)
rects = IntTrans().blh2xyz(lat, lon, H, ellipsoid_name)
print rects
print IntTrans().ell2gauss(lat, lon, H, ellipsoid_name)
print 'Transformed:'
rerects = ExtTrans().recthelmert(rects[0], rects[1], rects[2], 'WGS84', 'SK42')
print rerects
againell = IntTrans().xyz2bl(rects[0], rects[1], rects[2], ellipsoid_name)
print 'againell', againell