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
    ro = 206264.806



class Calc:  # Функции с вычислениями

    def __init__(self):
        pass

    def geocentrtoellipsoidal(self, x_big, y_big, z_big):
        # Вычисление R
        r_big = (x_big**2 + y_big**2)**0.5
        # Постепенное вычисление r
        r_small_p1 = z_big**2
        r_small_p2 = (x_big**2 + y_big**2)
        r_small_p3 = 1 - Const().e_2
        r_small = (r_small_p1 + r_small_p2*r_small_p3)**0.5
        # Постепенное вычисление тангенса широты
        tglat_ch_1 = z_big
        tglat_ch_2 = r_small**3 + Const().b*Const().es_2*(z_big**2)
        tglat_zn_1 = r_big
        tglat_zn_2 = r_small**3 - Const().b*Const().e_2*(1 - Const().e_2)*(r_big**2)
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

    def ell2gauss(self, B, L, H):
        Bsec = B*3600
        n = int((6+L)/6)
        l_part1 = 6*(n - 1)
        l_part2 = 3 + l_part1
        l_part3 = L - l_part2
        l = l_part3/57.29577951
        a = Const().a
        alpha = Const().alpha
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

lat1 = Calc().geocentrtoellipsoidal(Const().x_big_1, Const().y_big_1, Const().z_big_1)[1]
lon1 = Calc().geocentrtoellipsoidal(Const().x_big_1, Const().y_big_1, Const().z_big_1)[3]
hgt1 = Calc().geocentrtoellipsoidal(Const().x_big_1, Const().y_big_1, Const().z_big_1)[4]

lat2 = Calc().geocentrtoellipsoidal(Const().x_big_2, Const().y_big_2, Const().z_big_2)[1]
lon2 = Calc().geocentrtoellipsoidal(Const().x_big_2, Const().y_big_2, Const().z_big_2)[2]
hgt2 = Calc().geocentrtoellipsoidal(Const().x_big_2, Const().y_big_2, Const().z_big_2)[4]

lat3 = Calc().geocentrtoellipsoidal(Const().x_big_3, Const().y_big_3, Const().z_big_3)[1]
lon3 = Calc().geocentrtoellipsoidal(Const().x_big_3, Const().y_big_3, Const().z_big_3)[2]
hgt3 = Calc().geocentrtoellipsoidal(Const().x_big_3, Const().y_big_3, Const().z_big_3)[4]

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

n1 = Calc().ell2gauss(np.radians(lat1), np.radians(lon1), hgt1)[0]
x1 = Calc().ell2gauss(np.radians(lat1), np.radians(lon1), hgt1)[1]
y1 = Calc().ell2gauss(np.radians(lat1), np.radians(lon1), hgt1)[2]

n2 = Calc().ell2gauss(np.radians(lat2), np.radians(lon2), hgt2)[0]
x2 = Calc().ell2gauss(np.radians(lat2), np.radians(lon2), hgt2)[1]
y2 = Calc().ell2gauss(np.radians(lat2), np.radians(lon2), hgt2)[2]

n3 = Calc().ell2gauss(np.radians(lat3), np.radians(lon3), hgt3)[0]
x3 = Calc().ell2gauss(np.radians(lat3), np.radians(lon3), hgt3)[1]
y3 = Calc().ell2gauss(np.radians(lat3), np.radians(lon3), hgt3)[2]

print 'Пункт 1'
print 'Зона = ', lat1
print 'x1 = ', lon1
print 'y1 = ', hgt1
print '-----------'
print 'Пункт 2'
print 'Зона = ', lat2
print 'x2 = ', lon2
print 'y2 = ', hgt2
print '-----------'
print 'Пункт 3'
print 'Зона = ', lat3
print 'x3 = ', lon3
print 'y3 = ', hgt3
print '=========='



