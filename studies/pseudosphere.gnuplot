#
#
set terminal jpeg size 500,500
set output "pseudosphere.jpg"

set xrange [-1:1]
set yrange [-1:1]
set zrange [0:3]
set isosample 60
set ticslevel 0
set view equal xyz
unset tics
unset key
unset border

#r(x,y) = sqrt(x**2 + y**2)
#set mapping spherical
#splot -sqrt(1 - r(x,y)**2) - log((1 - sqrt(1 - r(x,y)**2))/r(x,y)) title 'Pseudosphere'
#splot -sqrt(1 - r**2) - log((1 - sqrt(1 - r**2))/r) title 'Pseudosphere'

set multiplot
set parametric
set urange [0:pi/2-0.2]
set vrange [0:11*pi/6]

alpha1 = 0
alpha2 = 11*pi/6
line_width = 4

splot [t=alpha1:alpha2] cos(t),sin(t),0 lt rgb 'green' lw line_width notitle
splot sin(u)*cos(v), sin(u)*sin(v), -log(tan(u/2)) - cos(u) lt rgb 'dark-grey' notitle
splot [u=0:pi/2] sin(u)*cos(alpha1), sin(u)*sin(alpha1), -log(tan(u/2)) - cos(u) lt rgb 'blue' lw  line_width notitle
splot [u=0:pi/2] sin(u)*cos(alpha2), sin(u)*sin(alpha2), -log(tan(u/2)) - cos(u) lt rgb 'red' lw  line_width notitle

set nomultiplot
