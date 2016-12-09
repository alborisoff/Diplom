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
    es_2 = e_2**2/(1 - e_2)  # Квадрат второго эксцентриситета

class Calc:  # Функции с вычислениями

    def __init__(self):
        pass




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
print 'КООРДИНАТЫ ПУНКТОВ'
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


