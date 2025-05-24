set terminal pngcairo enhanced font 'arial,10' fontscale 1.0 size 700, 500
set output 'DispPError.png'
set title 'SNR vs Bit Error Probability for QAM (p=1)'
set xlabel 'SNR (dB)'
set ylabel 'Bit Error Probability'
set logscale y
set grid
set key top right
plot 'qam_results(QAM16).dat' with linespoints title 'Simulated Bit Error', \
     3*erfc(sqrt(10**(x/10)))*(2**(1+1)-3*erfc(sqrt(10**(x/10)))/2)/(4**2) title 'Theoretical Upper Bound', \
     3*erfc(sqrt(10**(x/10)))*(2**(1+1)-3*erfc(sqrt(10**(x/10)))/2)/(4**2)/2 title 'Theoretical Lower Bound', \
