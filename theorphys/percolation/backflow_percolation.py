import math
#import matplotlib.pyplot as plt
#import matplotlib.patches as patches
from PIL import Image, ImageDraw
import numpy as np
import pprint


def check_percolation(array):
    percolation_map = np.zeros(array.size).reshape(array.shape)
    nx, ny = array.shape
    
    percolation_map[:,0] = array[:,0]
    percolation_exists = True

    stop = False
    counter = 0
#    print array
#    print percolation_map
    while not stop:# and counter < 2:
        stop = True 
#        counter += 1
        for i in xrange(nx):
            for j in xrange(ny):
                is_neighbor = ((i > 0    and percolation_map[i-1,j] == 1) or
                               (i < nx-1 and percolation_map[i+1,j] == 1) or
                               (j > 0    and percolation_map[i,j-1] == 1) or
                               (j < ny-1 and percolation_map[i,j+1] == 1))
#                print (i,j), is_neighbor
                if is_neighbor and array[i,j] == 1 and percolation_map[i,j] == 0: 
                    stop = False
                    percolation_map[i,j] = 1
#        print array
#        print percolation_map

    return percolation_exists, percolation_map 


def add_rectangle(draw, x, y, width, heigth, color):
    draw.rectangle([x, y, x + width, y + heigth], fill=color)

def plot_percolation(array, pmap):
    rect_size = 10
    im_size = rect_size*array.shape[0], rect_size*array.shape[1] 
    
    im = Image.new('RGBA', im_size, (0,0,0,0))
    draw = ImageDraw.Draw(im)

    nx, ny = array.shape

    for i in xrange(nx):
        for j in xrange(ny):
            x,y = i*rect_size, j*rect_size
            if pmap[i,j] == 1:
                color = 'blue'
            elif array[i,j] == 1:
                color = 'white'
            else:
                color = 'black'
            add_rectangle(draw, x, y, rect_size, rect_size, color)
    im.save('new_fig_1.png')
    


if __name__ == '__main__':
#    array = np.array([[1,1,1,0,0],
#                      [0,0,1,0,0],
#                      [0,1,1,0,0],
#                      [0,1,0,0,0],
#                      [0,1,1,1,1]])

    pcrit = (math.sqrt(5)-1)/2
    delta_p = -0.02

#    L = int(math.pow(1/abs(delta_p), 1/math.log(4*(pcrit - pcrit**3), 2)))
#    L = 10
#    nx = ny = L
    nx = 100
    ny = 100

#    print 'length: {}'.format(L)
    p = 0.5 #pcrit + delta_p
    array = np.random.binomial(1, pcrit + delta_p,size=nx*ny).reshape(nx,ny)
    percolation, pmap = check_percolation(array)
    plot_percolation(array, pmap)
