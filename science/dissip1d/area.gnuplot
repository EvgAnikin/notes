set terminal jpeg
set output 'area.jpg'

tt = 1.

f(omega,epsilon2) = acos(1./2/tt*(omega - epsilon2/sqrt(4*tt**2 - omega**2)))
#set parametric
#set dummy t
#set xrange [0:pi]
#set yrange [-2*tt:2*tt]
#set trange [-2*tt:2*tt]
#plot f(t, 0.2), t title '0.01', f(t, -0.2), t title '-0.1'


set xrange [-2*tt:2*tt]
set samples 1000
plot f(x,0.2), f(x, -0.2)
