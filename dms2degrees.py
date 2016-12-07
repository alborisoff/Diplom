# -*- coding: utf-8 -*-


def dmstodegrees(dms):
    dmssplit = dms.split(' ')
    degrees = float(dmssplit[0]) + float(dmssplit[1])/60 + float(dmssplit[2])/3600
    return degrees


print dmstodegrees('55 37 34.4543223')