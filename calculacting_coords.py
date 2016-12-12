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
print 'B1 = ', Calc().geocentrtoellipsoidal(Const().x_big_1, Const().y_big_1, Const().z_big_1)[1]
print 'L1 = ', Calc().geocentrtoellipsoidal(Const().x_big_1, Const().y_big_1, Const().z_big_1)[3]
print 'H1 = ', Calc().geocentrtoellipsoidal(Const().x_big_1, Const().y_big_1, Const().z_big_1)[4]
print '-----------'
print 'Пункт 2'
print 'B2 = ', Calc().geocentrtoellipsoidal(Const().x_big_2, Const().y_big_2, Const().z_big_2)[1]
print 'L2 = ', Calc().geocentrtoellipsoidal(Const().x_big_2, Const().y_big_2, Const().z_big_2)[3]
print 'H2 = ', Calc().geocentrtoellipsoidal(Const().x_big_2, Const().y_big_2, Const().z_big_2)[4]
print '-----------'
print 'Пункт 3'
print 'B3 = ', Calc().geocentrtoellipsoidal(Const().x_big_3, Const().y_big_3, Const().z_big_3)[1]
print 'L3 = ', Calc().geocentrtoellipsoidal(Const().x_big_3, Const().y_big_3, Const().z_big_3)[3]
print 'H3 = ', Calc().geocentrtoellipsoidal(Const().x_big_3, Const().y_big_3, Const().z_big_3)[4]
print '=========='




