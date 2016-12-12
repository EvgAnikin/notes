set terminal jpeg
set output 'levels.jpg'

xi = 2.
tt = 2.
dE = -81
dt = -9
i = {0,1}

#omega(x) = -sqrt(xi**2 + 2*tt**2 + 2*tt**2*cosh(x))
#g(x) = -2*tt**2*sinh(x) + dE*(omega(x) + xi)
#h(x) = (dt*(1 + exp(-x)) - tt*sinh(x))**2 \
#			- 0.5*dE*sinh(x)*(omega(x) + xi) \
#			- dt**2*(1 + cosh(x))*(1 + exp(-x))

#f(x) is equal to 2*t**2*h(x)
#f(x) = tt*(dt*(1 + exp(-x)) - tt*sinh(x))* \
#	(2*tt*dt*(1 + exp(-x)) - 2*tt**2*sinh(x) + dE*(omega(x) + xi)) - \
#	dt*(omega(x) + xi)*(1 + exp(-x))*(dt*(omega(x) - xi) + dE*tt)

ch(x) = 1 + (xi**2 - x**2)/2./tt
sh(x) = sqrt(ch(x)**2 - 1)
ex(x) = ch(x) + sh(x)

r(x) = (dt*(1 - ex(-x)) + tt*sh(x))**2 \
			+ 0.5*dE*sh(x)*(x + xi) \
			- dt**2*(1 - ch(x))*(1 - ex(-x))
set xrange[-xi:xi]
plot r(x)
