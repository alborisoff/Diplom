# -*- utf-8 -*-

def dmstodegrees(dms):
    degrees = dms[0:2]
    minutes = dms[3:5]
    seconds = dms[6:len(dms)]
    print degrees, minutes, seconds
    degros = float(degrees) + float(minutes)/60 + float(seconds)/3600
    return degros

print dmstodegrees('55 37 34.4543223')