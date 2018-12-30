clear;
tau = 0.01;
ifinal = 100;
A1 = [0 4;
     -3 2];
B1 = [0;
      1];
C1 = -[0.3 4];
D1 = 0;
sys2 = ss(A1,B1,C1,D1); 

Ts(1) = tau;
for i = 1:ifinal
    sd2_1 = c2d(sys2, Ts(i));
    Ad2 = sd2_1.a + sd2_1.b*C1;
    Bd2 = sd2_1.b;
    sd2_2 = ss(Ad2, Bd2, C1, D1, Ts(i));
    
    H2norm(i) = norm(sd2_2);
    
    for k = 0:1000
        h(k+1,i) = C1*Ad2^k*Bd2;
    end
    hmax(i) = max(abs(h(:,i)));
    L2norm(i) = norm(h(:,i));
    
    hnorm(i) = H2norm(i)^2;
    epsilon(i) = hmax(i)/H2norm(i);
    
    Ts(i+1) = Ts(i)/2;
end
figure
plot(Ts(1:100), hnorm)
figure
plot(Ts(1:100), epsilon)
