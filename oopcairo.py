#!/usr/bin/env python

import cairo
import numpy as np

def convertrgb(rgb_list):
    '''this function converts rgb 0-255 to rgb 0-1 values'''

    rgb = np.array(rgb_list)/255
    return rgb

#setting dimensions and rectangle color
WIDTH, HEIGHT = 500, 500
rectcol = [191, 157, 218] #purple :)
linecol = [123, 123, 123] #gray :)

#establish canvas
surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
ctx = cairo.Context(surface)

#line
lcol = convertrgb(linecol)
ctx.set_source_rgba(lcol[0], lcol[1], lcol[2])
ctx.set_line_width(1)
ctx.move_to(50,250)        #(x,y)
ctx.line_to(450,250)
ctx.stroke()

#rectangle
ctx.rectangle(150, 225, 150, 50)  # Rectangle(x0, y0, x1, y1)
rcol = convertrgb(rectcol)
ctx.set_source_rgb(rcol[0], rcol[1], rcol[2])
ctx.fill()

#saving png
surface.write_to_png("ooptest.png")  # Output to PNG