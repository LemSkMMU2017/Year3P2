set term jpeg size 1024,720
set output "E.jpeg"
plot "dotsE.dat" using 1:2 title "E" with lines
