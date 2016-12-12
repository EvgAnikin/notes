set terminal epslatex color dashed
set output 'theorphys-gnuplottex-fig1.tex'
set xrange [0:3]
set yrange [0:1.3]
set xlabel '$\beta J$'
set ylabel '$-\frac{U}{J}$'
plot tanh(x) title "Одна цепочка",\
     1./3*(((4*sinh(2*x)*cosh(2*x)*cosh(x)**2) + \
2*sinh(x)*cosh(2*x)**2*cosh(x) - \
4*sinh(2*x)*cosh(2*x)) /  \
(2*sqrt(cosh(x)**2*cosh(2*x)**2 - sinh(2*x)**2)) \
+ sinh(x)*cosh(2*x)+2*sinh(2*x)*cosh(x))\
/(cosh(x)*cosh(2*x)+sqrt(cosh(x)**2*cosh(2*x)**2 - sinh(2*x)**2)) title "Две цепочки"
