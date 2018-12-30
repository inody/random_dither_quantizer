clear;
close all;

set(0,'defaultAxesFontSize',14)
set(0,'defaultTextFontSize',14)

Ts = 0.05;
T = 10;
N = 5000;
d = 2;

A1 = [-15 6;
      3 -15];
B1 = [0;
      1];
C1 = [1 0];
sys1 = ss(A1,B1,C1,0);
sd1 = c2d(sys1,Ts);
y = impulse(sd1,2)
stem(y)