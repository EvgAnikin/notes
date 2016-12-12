set terminal jpeg
set output "gamma_asymptotics.jpg"

set logscale y
set xrange [0.1:10]

stirling(x) = x**x*exp(-x)*sqrt(2*pi*x)

plot gamma(x+1), stirling(x), stirling(x)*(1 + 1./12/x)
