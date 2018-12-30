clear;

delta = 0.01;

g = 9.81;
m = 0.63;           % cart weight (0.466+0.110)
M = 2.25;           % seesaw weight

kt = 0.023;         % torque constant
rm = 0.006;         % radious of motor gear
Ra = 8.3;           % resistance

kf= kt/(Ra*rm);     % volt 2 force

ls = 0.08595;       
lc = 0.12000;  
% ls = 0.069;
% lc = 0.098;

J = 0.21731;
Jbar = J + m*lc^2;

mbar = 0.9253;
ccbar= 2.0023;
cs   = 0;

a1=-m*m*lc*g/(Jbar*mbar-m*m*lc*lc);
a2=-Jbar*ccbar/(Jbar*mbar-m*m*lc*lc);
a3=(Jbar*m-m*lc*(m*lc+M*ls))/(Jbar*mbar-m*m*lc*lc)*g;
a4=m*lc*ccbar/(Jbar*mbar-m*m*lc*lc);
a5=mbar*m*g/(Jbar*mbar-m*m*lc*lc);
a6=m*lc*ccbar/(Jbar*mbar-m*m*lc*lc);
a7=-(m*m*lc-mbar*(m*lc+M*ls))/(Jbar*mbar-m*m*lc*lc)*g;
a8=-mbar*cs/(Jbar*mbar-m*m*lc*lc);

b1=Jbar*kf/(Jbar*mbar-m*m*lc*lc);
b2=-kf*m*lc/(Jbar*mbar-m*m*lc*lc);

Ac=[0 1 0 0; 
    a1 a2 a3 a4;
    0 0 0 1;
    a5 a6 a7 a8];

Bc=[0;
    b1;
    0;
    b2];
    
Cc=[1 0 0 0;
    0 0 1 0];

Dc=0;

rank(ctrb(Ac,Bc));
rank(obsv(Ac,Cc));

sys_seesaw = ss(Ac, Bc, Cc, Dc);
c2d_seesaw = c2d(sys_seesaw,delta);
Ad=c2d_seesaw.a;
Bd=c2d_seesaw.b;
Cd=c2d_seesaw.c;
Dd=c2d_seesaw.d;

Q = [10 0 0 0;
     0 1000 0 0;
     0 0 10000 0;
     0 0 0 1];
R = 1;
K = dlqr(Ad,Bd,Q,R)
% Pole_ctrl = 6*[0.1+0.1j, 0.1-0.1j, -0.1+0.1j, -0.1-0.1j];
% K = place(Ad, Bd, Pole_ctrl)
% ctrltest = eig(Ad-Bd*K)

Pole_obsv = 0.5*[1.01, 1.02, 1.03, 1.04];
% Pole_obsv = 0.1*[4, 3, 2, 1];
L = place(Ad', Cd', Pole_obsv)'
% obsvtest = eig(Ad-L*Cd)

Aob = Ad-Bd*K-L*Cd;
Bob = L;
Cob = -K;
Dob = zeros(1,2);

A_tilde = [Ad, zeros(size(Ad));
           Bob*Cd, Aob];
B_tilde = [Bd;
           zeros(size(Bd))];
C_tilde = [zeros(size(Cob)), Cob];

sys_tilde = ss(A_tilde+B_tilde*C_tilde,B_tilde,C_tilde,0,1);
norm_sys_tilde = norm(sys_tilde)^2;

K = 800;
N = 1000;
d = 4;
w = d*(rand(N,K)-0.5);

for n = 1:N
    n
%     x_tilde(:,1) = rand(8,1);
    x_tilde(:,1) = ones(8,1);

    for k = 1:K
        y_tilde(n,k) = C_tilde*x_tilde(:,k);
        u_tilde(k) = d*round((y_tilde(n,k) + w(n,k))/d) - w(n,k);
        x_tilde(:,k+1) = A_tilde*x_tilde(:,k) + B_tilde*u_tilde(k);
    end
end
for k = 1:K
    Vy(k) = var(y_tilde(:,k));
    Vy_bound(k) = norm_sys_tilde*d^2/4;
    Vy_approx(k) = 0.2415;
    Vy_subtractive(k) = norm_sys_tilde*d^2/12;
end


figure
hold on
plot(Vy,'Color', [0,0,0], 'LineWidth', 2)
plot(Vy_bound,'Color', [1,0,0], 'LineWidth', 2)
plot(Vy_approx, 'Color', [0,0,1], 'LineWidth', 2)
plot(Vy_subtractive, 'Color', [0,1,0], 'LineWidth', 2)
xlim([0 K]);
hold off


