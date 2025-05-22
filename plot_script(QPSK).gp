set terminal pngcairo enhanced font 'arial,10' fontscale 1.0 size 700, 500
set output 'DispPError.png'
set title 'D vs Bit Error Probability for QAM (p=0)'
set xlabel 'D (dB)'
set ylabel 'Bit Error Probability'
set logscale y
set grid
set key top right
plot 'qam_results(QPSK).dat' with linespoints title 'Simulated Bit Error', \
     0.5*erfc(sqrt(0.5/(10**(x/10)))) title 'Theoretical BER (p=0)'
