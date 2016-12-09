# -*- coding: utf-8 -*-
import yaml as ym
import numpy as np


class Const:

    def __init__(self):
        self.ro = 206264.806


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
            elldict['e2'] = 2*elldict['alpha'] - (elldict['alpha'])**2
            elldict['es2'] = elldict['e2']/(1 - elldict['e2'])
        return elldict

    def bignm(self, a, e2, lat):
        sin_bigb = np.sin(np.radians(lat))
        subradical = 1 - e2*(sin_bigb**2)
        readyn = a/(subradical**0.5)
        readym = readyn*(1 - e2)
        return [readyn, readym]

    def ellandsystem(self, system_name):
        ellnsys = self.unpack('ellsandsystems.yml')
        ellipsoid_name = ''
        for ell in ellnsys:
            if system_name == ellnsys[ell]:
                ellipsoid_name = ell
        return ellipsoid_name

    def transform_params(self, oldsystem, newsystem):
        # Определение семи параметров трансформирования координат.
        # По введённым системам координат получаем
        params = self.unpack('transes.yml')
        params_dict = {}
        # По входным названиям систем координат определяем эллипсоиды, на которых эти системы основаны
        old_ellipsoid = self.ellipsoid_params(self.ellandsystem(oldsystem))
        new_ellipsoid = self.ellipsoid_params(self.ellandsystem(newsystem))
        params_dict['old_a'] = old_ellipsoid['a']
        params_dict['old_alpha'] = old_ellipsoid['alpha']
        params_dict['old_e2'] = old_ellipsoid['e2']
        params_dict['new_a'] = new_ellipsoid['a']
        params_dict['new_alpha'] = new_ellipsoid['alpha']
        params_dict['new_e2'] = new_ellipsoid['e2']
        name_of_transform = oldsystem + '--' + newsystem
        if name_of_transform in params:
            params_dict['dx'] = eval(params[name_of_transform]['dx'])
            params_dict['dy'] = eval(params[name_of_transform]['dy'])
            params_dict['dz'] = eval(params[name_of_transform]['dz'])
            params_dict['wx'] = eval(params[name_of_transform]['wx'])
            params_dict['wy'] = eval(params[name_of_transform]['wy'])
            params_dict['wz'] = eval(params[name_of_transform]['wz'])
            params_dict['m'] = eval(params[name_of_transform]['m'])
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
        # Синус и косинус широты и долготы
        sinB = np.sin(np.radians(B))
        sinL = np.sin(np.radians(L))
        cosB = np.cos(np.radians(B))
        cosL = np.cos(np.radians(L))
        e2 = Srv().ellipsoid_params(ellipsoid_name)['e2']  # Квадрат первого эксцентриситета
        N = Srv().bignm(a, e2, B)[0]
        # Вычисление X, Y, Z
        X = (N + H) * cosB * cosL
        Y = (N + H) * cosB * sinL
        Z = (H + N) * sinB - e2 * N * sinB
        # Выдача результата в виде массива [X, Y, Z]
        return [X, Y, Z]

    def ell2gauss(self, B, L, H, ellipsoid_name):
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
        sin2B = np.sin(np.radians(B*2.0))
        cosB = np.cos(np.radians(B))
        tgB = np.tan(np.radians(B))
        sinB2 = sinB**2  # Квадрат синуса широты
        e2 = 2*alpha - alpha**2  # Квадрат первого эксцентриситета
        es2 = e2/(1 - e2)
        subradical = 1 - e2*sinB2  # Подкоренное значение для вычисления радиуса кривизны первого вертикала
        N = a/(subradical**0.5)  # Вычисление радиуса первого вертикала
        # Постепенный пересчёт эллиптических координат в плоские прямоугольные Гаусса.
        a2 = 0.5*(N*sinB*cosB)
        etha = cosB*(es2**0.5)
        a4 = (1/24)*(N*sinB*(cosB**3))*(5 - (tgB**2) + 9*(etha**2) + 4*(etha**4))
        a6 = (1/720)*(N*sinB*(cosB**5))*(61 - 58*(tgB**2) + tgB**4 + 270*(etha**2) - 330*(etha**2)*(tgB**2))
        a8 = (1/40320)*(N*sinB*(cosB**7))*N*sinB*(cosB**7)*(1385 -  3111*(tgB**2) + 543*(tgB**4) - (tgB**6))
        ro = 206264.806
        b1 = N*cosB
        b3 = (1/6)*N*(cosB**3)*(-(tgB**2) + etha**2)
        b5 = (1/120)*N*(cosB**5)*(5 - 18*(tgB**2) + tgB**4 - 14*(etha**2) - 58*(etha**3)*(tgB**2))
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

    # Пересчёт прямоугольных координат по формуле Гельмерта
    def recthelmert(self, bigx_old, bigy_old, bigz_old, old_system, new_system):
        # Загрузка семи параметров пересчёта
        tparameters = Srv().transform_params(old_system, new_system)
        dx = tparameters['dx']
        dy = tparameters['dy']
        dz = tparameters['dz']
        wx = tparameters['wx']
        wy = tparameters['wy']
        wz = tparameters['wz']
        m = tparameters['m']
        # Подготовка матриц пересчёта
        oldmatrix = np.matrix([[bigx_old],
                               [bigy_old],
                               [bigz_old]])
        deltamatrix = np.matrix([[dx],
                                 [dy],
                                 [dz]])
        omegamatrix = np.matrix([[1, wz, -wy],
                                 [-wz, 1, wx],
                                 [wy, -wx, 1]])
        print '***'
        print 'Параметры преобразования'
        print 'Масштабный коэффициент m', m
        print 'Матрица разностей bigx_old, bigy_old, bigz_old'
        print deltamatrix
        print ''
        print 'Матрица разворотов координатных осей'
        print omegamatrix
        print '***'
        newmatrix = deltamatrix + (1 + m)*omegamatrix*oldmatrix
        result = np.squeeze(np.asarray(newmatrix))
        return [result[0], result[1], result[2]]

    def molodensky(self, lat_old, lon_old, hgt_old, oldsystem, newsystem):
        # Актуализировать формулы в соответствии с
        # http://cyberleninka.ru/article/n/metody-transformatsii-geodezicheskih-i-prostranstvennyh-pryamougolnyh-koordinat-ih-algoritmy-parametry-tochnost
        tpar = Srv().transform_params(oldsystem, newsystem)
        delta_a = tpar['new_a'] - tpar['old_a']
        delta_e2 = tpar['new_e2'] - tpar['old_e2']
        sinlatold = np.sin(np.radians(lat_old))
        coslatold = np.cos(np.radians(lat_old))
        coslat2old = np.sin(np.radians(2.0 * coslatold))
        tanlatold = np.tan(np.radians(lat_old))
        sinlatold_v2 = sinlatold**2
        sinlonold = np.sin(np.radians(lon_old))
        coslonold = np.cos(np.radians(lon_old))
        a_mid = (tpar['new_a'] + tpar['old_a'])/2
        e2_mid = (tpar['new_e2'] + tpar['old_e2'])/2
        big_n = Srv().bignm(a_mid, e2_mid, lat_old)[0]
        big_m = Srv().bignm(a_mid, e2_mid, lat_old)[1]
        # Поэтапное вычисление поправки широты.
        delta_b_p1 = (Const().ro/(big_m + hgt_old))
        delta_b_p2 = (big_n/a_mid)*e2_mid*sinlatold*coslatold*delta_a
        delta_b_p3 = (big_n**2/a_mid**2 + 1)*big_n*sinlatold*coslatold*(delta_e2/2)
        delta_b_p4 = sinlatold*(tpar['dx']*coslonold + tpar['dy']*sinlonold) + tpar['dz']*coslatold
        delta_b_p5 = tpar['wx']*sinlonold*(1 - e2_mid*coslat2old)
        delta_b_p6 = tpar['wy']*coslonold*(1 + e2_mid*coslat2old)
        delta_b_p7 = Const().ro*tpar['m']*e2_mid*sinlatold*coslat2old
        delta_b = delta_b_p1*(delta_b_p2 + delta_b_p3 - delta_b_p4) - delta_b_p5*delta_b_p6 - delta_b_p7
        delta_l_p01 = Const().ro/(coslatold*(big_n + hgt_old))
        delta_l_p02 = (-tpar['dx']*sinlonold + tpar['dy']*coslonold)
        delta_l_p11 = delta_l_p01*delta_l_p02
        delta_l_p03 = tanlatold*(1 - e2_mid)*tpar['wx']*coslatold
        delta_l_p04 = tpar['wy']*sinlonold - tpar['wz']
        delta_l = delta_l_p11 + delta_l_p03 + delta_l_p04
        delta_h_p1 = -(a_mid/big_n)*delta_a
        delta_h_p2 = big_n*sinlatold_v2*(delta_e2/2)
        delta_h_p3 = tpar['dx']*coslonold
        delta_h_p4 = tpar['dy']*sinlonold*coslatold
        delta_h_p5 = tpar['dz']*sinlatold
        delta_h_p6 = big_n*e2_mid*sinlatold*coslatold
        delta_h_p7 = (tpar['wx']/Const().ro)*sinlonold - (tpar['wy']/Const().ro)*coslonold
        delta_h_p8 = (a_mid**2/big_n + hgt_old)*tpar['m']
        delta_h = delta_h_p1 + delta_h_p2 + delta_h_p3 + delta_h_p4 + delta_h_p5 - delta_h_p6*delta_h_p7 + delta_h_p8
        lat_new = lat_old + delta_b
        lon_new = lon_old + delta_l
        hgt_new = hgt_old + delta_h
        return [lat_new, lon_new, hgt_new]


ellipsoid_name = 'WGS84'
print 'ВХОДНЫЕ КООРДИНАТЫ (WGS84)'
print 'Введите геоцентрические коордиаты пункта.'
x_big = float(input('X: '))
y_big = float(input('Y: '))
z_big = float(input('Z: '))
print '----------'
latlon = IntTrans().xyz2bl(x_big, y_big, z_big, ellipsoid_name)
print 'Криволинейные координаты пункта'
print 'B = ', latlon[0]
print 'L = ', latlon[1]
print '----------'
print 'Плоские прямоугольные Гаусса'
gauss =  IntTrans().ell2gauss(lat, lon, height, ellipsoid_name)
print 'Зона: ', gauss[0]
print 'x: ', gauss[1]
print 'y: ', gauss[2]
print '=========='
print 'ТРАНСФОРМИРОВАННЫЕ КООРДИНАТЫ (WGS84 => PZ90)'
print 'Геоцентрические по 7-параметрической формуле Гельмерта'
rerects_wgs_pz = ExtTrans().recthelmert(rects[0], rects[1], rects[2], 'WGS84', 'PZ90')
print 'X1: ', rerects_wgs_pz[0]
print 'Y1: ', rerects_wgs_pz[1]
print 'Z1: ', rerects_wgs_pz[2]
print '----------'
print 'Эллипсоидальные по формуле Молоденского'
mldn_wgs_pz = ExtTrans().molodensky(lat, lon, height, 'WGS84', 'PZ90')
print 'B1: ', mldn_wgs_pz[0]
print 'L1: ', mldn_wgs_pz[1]
print 'H1: ', mldn_wgs_pz[2]
print '=========='
print 'ТРАНСФОРМИРОВАННЫЕ КООРДИНАТЫ (PZ90 => SK42)'
print 'Геоцентрические по 7-параметрической формуле Гельмерта'
rerects_pz_sk42 = ExtTrans().recthelmert(rects[0], rects[1], rects[2], 'PZ90', 'SK42')
print 'X1: ', rerects_pz_sk42[0]
print 'Y1: ', rerects_pz_sk42[1]
print 'Z1: ', rerects_pz_sk42[2]
print '----------'
print 'Эллипсоидальные по формуле Молоденского'
mldn_pz_sk42 = ExtTrans().molodensky(lat, lon, height, 'PZ90', 'SK42')
print 'B1: ', mldn_pz_sk42[0]
print 'L1: ', mldn_pz_sk42[1]
print 'H1: ', mldn_pz_sk42[2]


