clear;
ds = 0.01;
d = 2;
Ts = 0.01;
vs = 0.001;
vmax = 2;

A1 = [0 4;
     -3 2];
B1 = [0;
      1];
C1 = -[0.3 4];
sys1 = ss(A1,B1,C1,0);
sd1 = c2d(sys1,Ts);
Ad1 = sd1.a + sd1.b*C1;
Bd1 = sd1.b;
sd2 = ss(Ad1,Bd1,C1,0,Ts);
hl2 = norm(sd2);

v = 0:ds:d;
for i = 1: size(v,2)
    txi(i) = -v(i)^2+v(i)*d;
    for n = 1: 99
        txi(i+n*d/ds) = txi(i);
    end
end
w = -50*d:ds:50*d;
% plot(w,txi)

j=0;
for v = 0: vs: vmax;
    j=j+1
    for i = 1: size(w,2)
        fun(i) = txi(i)*exp(-w(i)^2/(2*v))/sqrt(2*pi*v);
    end
    i_v(j) = v;
    f_v(j) = trapz(w,fun);
end
plot(i_v,f_v)
% [y,sgmy] = min(abs(f_v - i_v/hl2^2));
% plot(abs(f_v - i_v/hl2^2))
% vy_approx = sgmy*vs