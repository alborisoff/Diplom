# -*- coding: utf-8 -*-

import numpy as np


class Const:

    def __init__(self):
        pass

    n_group = 1  # Номер группы
    n_my = 3  # Мой номер по журналу

    # Координаты X, Y, Z первой точки
    x_big_1 = 2780905.887 + 10*n_group
    y_big_1 = 2095562.894 + 10*n_group
    z_big_1 = 5326025.900 + 10*n_group

    # Приращения координат второй точки
    delta_x_1_2 = 120.000 + 0.2*n_my
    delta_y_1_2 = 1635.000 + 0.1*n_my
    delta_z_1_2 = 64.000 + 0.3*n_my

    # Приращения координат третьей точки
    delta_x_1_3 = 151.000 + 0.1*n_my
    delta_y_1_3 = 1600.000 + 0.3*n_my
    delta_z_1_3 = 90.000 + 0.2*n_my

    # Координаты второй точки
    x_big_2 = x_big_1 + delta_x_1_2
    y_big_2 = y_big_1 + delta_y_1_2
    z_big_2 = z_big_1 + delta_z_1_2

    # Координаты третьей точки
    x_big_3 = x_big_1 + delta_x_1_3
    y_big_3 = y_big_1 + delta_y_1_3
    z_big_3 = z_big_1 + delta_z_1_3

    # Параметры эллипсоида Красовского
    a = 6378245.0  # Большая полуось
    alpha = 1.0/298.3  # Сжатие
    e_2 = 2*alpha - alpha**2  # Квадрат первого эксцентриситета
    es_2 = e_2/(1 - e_2)  # Квадрат второго эксцентриситета
    b = a*((1 - e_2)**0.5)

    r_earth = 6371000  # Средний радиус Земли, м
    ro = 206264.806  # Количество секунд в радиане


class Calc:  # Функции с вычислениями

    def __init__(self):
        pass

    def geocentrtoellipsoidal(self, x_big, y_big, z_big, b_axis):
        # Вычисление R
        r_big = (x_big**2 + y_big**2)**0.5
        # Постепенное вычисление r
        r_small_p1 = z_big**2
        r_small_p2 = (x_big**2 + y_big**2)
        r_small_p3 = 1 - Const().e_2
        r_small = (r_small_p1 + r_small_p2*r_small_p3)**0.5
        # Постепенное вычисление тангенса широты
        tglat_ch_1 = z_big
        tglat_ch_2 = r_small**3 + b_axis*Const().es_2*(z_big**2)
        tglat_zn_1 = r_big
        tglat_zn_2 = r_small**3 - b_axis*Const().e_2*(1 - Const().e_2)*(r_big**2)
        tglat = (tglat_ch_1*tglat_ch_2)/(tglat_zn_1*tglat_zn_2)  # Вычисленный тангенс широты
        latrad = np.arctan(tglat)  # Широта в радианах
        latdeg = np.degrees(latrad)  # Широта в градусах
        tanlon = y_big/x_big  # Тангенс долготы
        lonrad = np.arctan(tanlon)  # Долгота в радианах
        londeg = np.degrees(lonrad)  # Долгота в градусах
        # Постепенное вычисление высоты
        h_big_p1 = r_big*np.cos(latrad)
        h_big_p2 = z_big*np.sin(latrad)
        h_big_p3 = 1 - Const().e_2*(np.sin(latrad))**2
        h_big = h_big_p1 + h_big_p2 - Const().a*(h_big_p3**0.5)
        # Возвращаем широту в радианах, широту в градусах, долготу в радианах, долготу в градусах и высоту.
        return [latrad, latdeg, lonrad, londeg, h_big]

    def ellipsoidaltogauss(self, latrad, lonrad, a_axis):
        latsec = np.degrees(latrad)*3600
        sinlat = np.sin(latrad)
        sinlat2 = sinlat**2
        sinlat4 = sinlat**4
        sin2lat = np.sin(2.0*latrad)
        coslat = np.cos(latrad)
        coslat3 = coslat**3
        coslat5 = coslat**5
        coslat7 = coslat**7
        tanlat = np.tan(latrad)
        tanlat2 = tanlat**2
        tanlat4 = tanlat**4
        tanlat6 = tanlat**6
        subradical = 1 - Const().e_2*sinlat2
        n_big = a_axis/(subradical**0.5)
        eta = (Const().es_2**0.5)*coslat
        n_zone = int((6 + np.degrees(lonrad)) / 6)
        l_part1 = 6 * (n_zone - 1)
        l_part2 = 3 + l_part1
        l_part3 = np.degrees(lonrad) - l_part2
        l = l_part3 / 57.29577951
        a_big_x = 1 + (3/4)*Const().e_2 + (45/64)*(Const().e_2**2)
        b_big_x = (3/4)*Const().e_2 + (15/16)*(Const().e_2**2)
        c_big_x = (15/64)*(Const().e_2**2)
        x_big_m1 = a_axis
        x_big_m2 = 1 - Const().e_2
        x_big_s31 = a_big_x*(latsec/Const().ro)
        x_big_s32 = (b_big_x/2)*sin2lat
        x_big_s33 = (c_big_x/4)*sinlat4
        x_big = x_big_m1*x_big_m2*(x_big_s31 + x_big_s32 + x_big_s33)
        a2 = (1/2)*n_big*sinlat*coslat
        a4 = (1/24)*n_big*sinlat*coslat3*(5 - tanlat2 + 9*(eta**2) + 4*(eta**4))
        a6 = (1/720)*n_big*sinlat*coslat5*(61 - 58*tanlat2 + tanlat4 + 270*(eta**2) - 330*(eta**2)*tanlat2)
        a8 = (1/40320)*n_big*sinlat*coslat7*(1385 - 3111*tanlat2 + 543*tanlat4 - tanlat6)
        b1 = n_big*coslat
        b3 = (1/6)*n_big*coslat3*(-tanlat2 + (eta**2))
        b5 = (1/120)*n_big*coslat5*(5 - 18*tanlat2 + tanlat4 - 14*(eta**2) - 58*(eta**2)*tanlat2)
        b7 = (1/5040)*n_big*coslat7*(61 - 479*tanlat2 + 179*tanlat4 - tanlat6)
        x_gauss = x_big + a2*(l**2) + a4*(l**4) + a6*(l**6) + a8*(l**8)
        y_gauss = b1*l + b3*(l**3) + b5*(l**5) + b7*(l**7) + 500000
        return [n_zone, x_gauss, y_gauss]
    
    def deltass(self, x_b_o, y_b_o, z_b_o, x_b_n, y_b_n, z_b_n, x_g_o, y_g_o, x_g_n, y_g_n):
        delta_x_b = x_b_n - x_b_o
        delta_y_b = y_b_n - y_b_o
        delta_z_b = z_b_n - z_b_o
        delta_x_g = x_g_n - x_g_o
        delta_y_g = y_g_n - y_g_o
        h_o = self.geocentrtoellipsoidal(x_b_o, y_b_o, z_b_o, Const().b)[4]
        h_n = self.geocentrtoellipsoidal(x_b_n, y_b_n, z_b_n, Const().b)[4]
        delta_h = h_n - h_o
        s_big = (delta_x_b**2 + delta_y_b**2 + delta_z_b**2)**0.5
        s_big_g = (s_big**2 - delta_h**2)**0.5
        s_big_g_k = (delta_x_g**2 + delta_y_g**2)**0.5
        delta_big_s = s_big_g - s_big_g_k
        return [s_big, s_big_g, s_big_g_k, delta_big_s]


# Итоговые вычисления для вывода информации

x_b_1 = Const().x_big_1
y_b_1 = Const().y_big_1
z_b_1 = Const().z_big_1
x_b_2 = Const().x_big_2
y_b_2 = Const().y_big_2
z_b_2 = Const().z_big_2
x_b_3 = Const().x_big_3
y_b_3 = Const().y_big_3
z_b_3 = Const().z_big_3
# Криволинейные координаты первого, второго и третьего пунктов
# В градусах
lat1 = Calc().geocentrtoellipsoidal(x_b_1, y_b_1, z_b_1, Const().b)[1]
lon1 = Calc().geocentrtoellipsoidal(x_b_1, y_b_1, z_b_1, Const().b)[3]
hgt1 = Calc().geocentrtoellipsoidal(x_b_1, y_b_1, z_b_1, Const().b)[4]
lat2 = Calc().geocentrtoellipsoidal(x_b_2, y_b_2, z_b_2, Const().b)[1]
lon2 = Calc().geocentrtoellipsoidal(x_b_2, y_b_2, z_b_2, Const().b)[3]
hgt2 = Calc().geocentrtoellipsoidal(x_b_2, y_b_2, z_b_2, Const().b)[4]
lat3 = Calc().geocentrtoellipsoidal(x_b_3, y_b_3, z_b_3, Const().b)[1]
lon3 = Calc().geocentrtoellipsoidal(x_b_3, y_b_3, z_b_3, Const().b)[3]
hgt3 = Calc().geocentrtoellipsoidal(x_b_3, y_b_3, z_b_3, Const().b)[4]
# В радианах
lat1rad = np.radians(lat1)
lon1rad = np.radians(lon1)
lat2rad = np.radians(lat2)
lon2rad = np.radians(lon2)
lat3rad = np.radians(lat3)
lon3rad = np.radians(lon3)
# Координаты первого, второго и третьего пунктов в проекции Гаусса
n1 = Calc().ellipsoidaltogauss(lat1rad, lon1rad, Const().a)[0]
x1 = Calc().ellipsoidaltogauss(lat1rad, lon1rad, Const().a)[1]
y1 = Calc().ellipsoidaltogauss(lat1rad, lon1rad, Const().a)[2]
n2 = Calc().ellipsoidaltogauss(lat2rad, lon2rad, Const().a)[0]
x2 = Calc().ellipsoidaltogauss(lat2rad, lon2rad, Const().a)[1]
y2 = Calc().ellipsoidaltogauss(lat2rad, lon2rad, Const().a)[2]
n3 = Calc().ellipsoidaltogauss(lat3rad, lon3rad, Const().a)[0]
x3 = Calc().ellipsoidaltogauss(lat3rad, lon3rad, Const().a)[1]
y3 = Calc().ellipsoidaltogauss(lat3rad, lon3rad, Const().a)[2]

# Вторая часть работы
h_big_midi = (hgt2 + hgt3)/2
k = 1 + (h_big_midi/Const().r_earth)
a_new = Const().a*k
b_new = Const().b*k

lat1_new = Calc().geocentrtoellipsoidal(x_b_1, y_b_1, z_b_1, b_new)[1]
lon1_new = Calc().geocentrtoellipsoidal(x_b_1, y_b_1, z_b_1, b_new)[3]
hgt1_new = Calc().geocentrtoellipsoidal(x_b_1, y_b_1, z_b_1, b_new)[4]
lat2_new = Calc().geocentrtoellipsoidal(x_b_2, y_b_2, z_b_2, b_new)[1]
lon2_new = Calc().geocentrtoellipsoidal(x_b_2, y_b_2, z_b_2, b_new)[3]
hgt2_new = Calc().geocentrtoellipsoidal(x_b_2, y_b_2, z_b_2, b_new)[4]
lat3_new = Calc().geocentrtoellipsoidal(x_b_3, y_b_3, z_b_3, b_new)[1]
lon3_new = Calc().geocentrtoellipsoidal(x_b_3, y_b_3, z_b_3, b_new)[3]
hgt3_new = Calc().geocentrtoellipsoidal(x_b_3, y_b_3, z_b_3, b_new)[4]
# В радианах
lat1rad_new = np.radians(lat1_new)
lon1rad_new = np.radians(lon1_new)
lat2rad_new = np.radians(lat2_new)
lon2rad_new = np.radians(lon2_new)
lat3rad_new = np.radians(lat3_new)
lon3rad_new = np.radians(lon3_new)
# Координаты первого, второго и третьего пунктов в проекции Гаусса
n1_new = Calc().ellipsoidaltogauss(lat1rad_new, lon1rad_new, a_new)[0]
x1_new = Calc().ellipsoidaltogauss(lat1rad_new, lon1rad_new, a_new)[1]
y1_new = Calc().ellipsoidaltogauss(lat1rad_new, lon1rad_new, a_new)[2]
n2_new = Calc().ellipsoidaltogauss(lat2rad_new, lon2rad_new, a_new)[0]
x2_new = Calc().ellipsoidaltogauss(lat2rad_new, lon2rad_new, a_new)[1]
y2_new = Calc().ellipsoidaltogauss(lat2rad_new, lon2rad_new, a_new)[2]
n3_new = Calc().ellipsoidaltogauss(lat3rad_new, lon3rad_new, a_new)[0]
x3_new = Calc().ellipsoidaltogauss(lat3rad_new, lon3rad_new, a_new)[1]
y3_new = Calc().ellipsoidaltogauss(lat3rad_new, lon3rad_new, a_new)[2]

s_big_1_2 = Calc().deltass(x_b_1, y_b_1, z_b_1, x_b_2, y_b_2, z_b_2, x1, x2, y1, y2)[0]
s_big_g_1_2 = Calc().deltass(x_b_1, y_b_1, z_b_1, x_b_2, y_b_2, z_b_2, x1, x2, y1, y2)[1]
s_big_g_k_1_2 = Calc().deltass(x_b_1, y_b_1, z_b_1, x_b_2, y_b_2, z_b_2, x1, x2, y1, y2)[2]
delta_big_s_1_2 = Calc().deltass(x_b_1, y_b_1, z_b_1, x_b_2, y_b_2, z_b_2, x1, x2, y1, y2)[3]

s_big_1_3 = Calc().deltass(x_b_1, y_b_1, z_b_1, x_b_3, y_b_3, z_b_3, x1, x3, y1, y3)[0]
s_big_g_1_3 = Calc().deltass(x_b_1, y_b_1, z_b_1, x_b_3, y_b_3, z_b_3, x1, x3, y1, y3)[1]
s_big_g_k_1_3 = Calc().deltass(x_b_1, y_b_1, z_b_1, x_b_3, y_b_3, z_b_3, x1, x3, y1, y3)[2]
delta_big_s_1_3 = Calc().deltass(x_b_1, y_b_1, z_b_1, x_b_3, y_b_3, z_b_3, x1, x3, y1, y3)[3]

s_big_2_3 = Calc().deltass(x_b_2, y_b_2, z_b_2, x_b_3, y_b_3, z_b_3, x2, x3, y2, y3)[0]
s_big_g_2_3 = Calc().deltass(x_b_2, y_b_2, z_b_2, x_b_3, y_b_3, z_b_3, x2, x3, y2, y3)[1]
s_big_g_k_2_3 = Calc().deltass(x_b_2, y_b_2, z_b_2, x_b_3, y_b_3, z_b_3, x2, x3, y2, y3)[2]
delta_big_s_2_3 = Calc().deltass(x_b_2, y_b_2, z_b_2, x_b_3, y_b_3, z_b_3, x2, x3, y2, y3)[3]

# Итоговый вывод информации
print 'ИСХОДНЫЕ ДАННЫЕ'
print 'Информация об эллипсоиде КРАСОВСКОГО:'
print 'Большая полуось: ', Const().a
print 'Сжатие: ', Const().alpha
print 'Квадрат первого эксцентриситета: ', Const().e_2
print 'Квадрат второго эксцентриситета: ', Const().es_2
print '=========='
print 'Номер группы: ', Const().n_group
print 'Номер студента по журналу: ', Const().n_my
print '=========='
print 'ИСХОДНЫЕ ГЕОЦЕНТРИЧЕСКИЕ КООРДИНАТЫ ПУНКТОВ'
print 'Пункт 1'
print 'X1 = ', Const().x_big_1
print 'Y1 = ', Const().y_big_1
print 'Z1 = ', Const().z_big_1
print '----------'
print 'Пункт 2'
print 'X2 = ', Const().x_big_2
print 'Y2 = ', Const().y_big_2
print 'Z2 = ', Const().z_big_2
print '----------'
print 'Пункт 3'
print 'X3 = ', Const().x_big_3
print 'Y3 = ', Const().y_big_3
print 'Z3 = ', Const().z_big_3
print '**********'
print 'ВЫЧИСЛЕННЫЕ КООРДИНАТЫ ПУНКТОВ'
print '=== Криволинейные (в градусах и долях градуса) ==='
print 'Пункт 1'
print 'B1 = ', lat1
print 'L1 = ', lon1
print 'H1 = ', hgt1
print '-----------'
print 'Пункт 2'
print 'B2 = ', lat2
print 'L2 = ', lon2
print 'H2 = ', hgt2
print '-----------'
print 'Пункт 3'
print 'B3 = ', lat3
print 'L3 = ', lon3
print 'H3 = ', hgt3
print '=========='
print '=== Плоские в проекции Гаусса ==='
print 'Пункт 1'
print 'Зона = ', n1
print 'x1 = ', x1
print 'y1 = ', y1
print '-----------'
print 'Пункт 2'
print 'Зона = ', n2
print 'x2 = ', x2
print 'y2 = ', y2
print '-----------'
print 'Пункт 3'
print 'Зона = ', n3
print 'x3 = ', x3
print 'y3 = ', y3
print '**********'
print 'ВЫЧИСЛЕНИЕ КООРДИНАТ НА НОВОМ ЭЛЛИПСОИДЕ'
print 'Информация о НОВОМ эллипсоиде:'
print 'Большая полуось: ', a_new
print 'Малая полуось: ', b_new
print 'Сжатие: ', Const().alpha
print 'Квадрат первого эксцентриситета: ', Const().e_2
print 'Квадрат второго эксцентриситета: ', Const().es_2
print '**********'
print 'ВЫЧИСЛЕННЫЕ КООРДИНАТЫ ПУНКТОВ'
print '=== Криволинейные (в градусах и долях градуса) ==='
print 'Пункт 1'
print 'B1 = ', lat1_new
print 'L1 = ', lon1_new
print 'H1 = ', hgt1_new
print '-----------'
print 'Пункт 2'
print 'B2 = ', lat2_new
print 'L2 = ', lon2_new
print 'H2 = ', hgt2_new
print '-----------'
print 'Пункт 3'
print 'B3 = ', lat3_new
print 'L3 = ', lon3_new
print 'H3 = ', hgt3_new
print '=========='
print '=== Плоские в проекции Гаусса ==='
print 'Пункт 1'
print 'Зона = ', n1_new
print 'x1 = ', x1_new
print 'y1 = ', y1_new
print '-----------'
print 'Пункт 2'
print 'Зона = ', n2_new
print 'x2 = ', x2_new
print 'y2 = ', y2_new
print '-----------'
print 'Пункт 3'
print 'Зона = ', n3_new
print 'x3 = ', x3_new
print 'y3 = ', y3_new
print '**********'
print 'ВЫЧИСЛЕНИЕ РАЗНОСТЕЙ РАССТОЯНИЙ'
print 'Hср = ', h_big_midi
print '----------'
print 'S(1 - 2) = ', s_big_1_2
print 'S(1 - 3) = ', s_big_1_3
print 'S(2 - 3) = ', s_big_2_3
print '----------'
print 'DELTA Sг(1 - 2) = ', s_big_g_1_2
print 'DELTA Sг(1 - 3) = ', s_big_g_1_3
print 'DELTA Sг(2 - 3) = ', s_big_g_1_2
print '----------'
print 'DELTA Sгк(1 - 2) = ', s_big_g_k_1_2
print 'DELTA Sгк(1 - 3) = ', s_big_g_k_1_3
print 'DELTA Sгк(2 - 3) = ', s_big_g_k_2_3
print '----------'
print 'DELTA S (1 - 2) = ', delta_big_s_1_2
print 'DELTA S (1 - 3) = ', delta_big_s_1_3
print 'DELTA S (2 - 3) = ', delta_big_s_2_3
print '***** КОНЕЦ ВЫЧИСЛЕНИЙ *****'
