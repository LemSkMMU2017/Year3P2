set term jpeg size 1024,720
set output "U.jpeg"
plot "funcU.txt" using 1:2 title "u" with lines

