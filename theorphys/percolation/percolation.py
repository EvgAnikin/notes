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

    for line in xrange(1,ny):
        for col in xrange(nx):
            if percolation_map[col,line-1] == 1 and array[col,line] == 1:
                percolation_map[col,line] = 1 

        if not 1 in percolation_map[:,line]:
            percolation_exists = False
            break
            

        col = 0 
        while True:
            while col < nx and array[col, line] == 0:
                col += 1
            if col == nx: 
                break
            cluster_start = col
            cluster_percolates= False
            while col < nx and array[col, line] == 1:
                if percolation_map[col, line] == 1:
                    cluster_percolates = True
                col += 1
            cluster_end = col
            if cluster_percolates:
                percolation_map[cluster_start:cluster_end,line] = \
                                np.ones(cluster_end - cluster_start)
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
                color = 'black'
            elif array[i,j] == 1:
                color = 'grey'
            else:
                color = 'white'
            add_rectangle(draw, x, y, rect_size, rect_size, color)
    im.save('new_fig.png')
    


if __name__ == '__main__':
#    array = np.array([[1,1,1,0,0],
#                      [0,0,1,0,0],
#                      [0,1,1,0,0],
#                      [0,1,0,0,0],
#                      [0,1,1,1,1]])

    pcrit = (math.sqrt(5)-1)/2
    delta_p = 0

#    L = int(math.pow(1/abs(delta_p), 1/math.log(4*(pcrit - pcrit**3), 2)))
    L = 400
    nx = ny = L

    print 'length: {}'.format(L)
    array = np.random.binomial(1, pcrit + delta_p,size=nx*ny).reshape(nx,ny)
    percolation, pmap = check_percolation(array)
    plot_percolation(array, pmap)
