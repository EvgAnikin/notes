xi = 0.3
m = 1
tt = -0.01
energy(px,py) = sqrt((xi + 1./m*(2 - cos(px) - cos(py)))**2 + 4*tt**2*(sin(px)**2 + sin(py)**2))

set xrange[-pi:pi]
set yrange[-pi:pi]
set isosample 50
set ticslevel 0
splot energy(x,y), -energy(x,y)
