xi = -0.01
m = 0.4
t = 1
set xrange[-pi:pi]
set yrange[-pi:pi]
set ticslevel 0
set isosample 50
splot sqrt((xi + 1./m*(2 - cos(x) - cos(y)))**2 + 4*t**2*(sin(x)**2 + sin(y)**2)) title "", \
     -sqrt((xi + 1./m*(2 - cos(x) - cos(y)))**2 + 4*t**2*(sin(x)**2 + sin(y)**2)) title ""
