clear;
close all;

set(0,'defaultAxesFontSize',18);
set(0,'defaultTextFontSize',18);
set(0,'defaultAxesLineWidth', 1.5);
set(0,'defaultLineLineWidth', 2);
clorder = ...
    [  0,   0, 255; % ê¬
       0, 128,   0; % óŒ
     255,   0,   0; % ê‘
       0, 255,   0;
     204,   8, 204; % éá
     222, 125,   0; % íÉ
       0,  51, 153; % çÆ (ê¬Ç∆ãÊï ÇµÇ√ÇÁÇ¢ÅCíçà”ÅI)
      64,  64,  64];% îZÇ¢äDêF  
set(0,'defaultAxesColorOrder',clorder/255);

Ts = 0.01;
T = 10;
N = 1000;
d = 2;

A1 = [0 4;
    -3 2];
B1 = [0;
      1];
C1 = -[0.3 4];
sys1 = ss(A1,B1,C1,0);
sd1 = c2d(sys1,Ts);
Ad1 = sd1.a;
Bd1 = sd1.b;
x(:,1) = 0.5*[-1;
           1];
x_2(:,1) = x(:,1);
w = d*(rand(N,T/Ts)-0.5);

for k = 1:T/Ts
    r(k) = 0.5*sin(k/20);
    r(k) = 0;
    yc(k) = C1*x(:,k) + r(k);  
    x(:,k+1) = Ad1*x(:,k) + Bd1*yc(k);
end

for k = 1:T/Ts
    yq(k) = C1*x(:,k);  
    uq(k) = d*round(yq(k)/d) + r(k);
    x(:,k+1) = Ad1*x(:,k) + Bd1*uq(k);
end

for n = 1:N
%     x(:,1) = [rand-0.5;
%               rand-0.5];
    countL = N-n
    for k = 1:T/Ts
        yd(k) = C1*x(:,k);
        ud(k) = d*round((yd(k)+w(n,k))/d) + r(k);
        x(:,k+1) = Ad1*x(:,k) + Bd1*ud(k);
    end
    y(:,n) = yd;
    u(:,n) = ud;
    eta(:,n) = ud - yd;
end

% for n = 1:N
% %     x(:,1) = [rand-0.5;
% %               rand-0.5];
%     countL = N-n
%     for k = 1:T/Ts
%         yd_2(k) = C1*x_2(:,k);
%         ud_2(k) = d*round((yd_2(k)+w(n,k))/d) - w(n,k) + r(k);
%         x_2(:,k+1) = Ad1*x_2(:,k) + Bd1*ud_2(k);
%     end
%     y_2(:,n) = yd_2;
%     u_2(:,n) = ud_2;
%     eta_2(:,n) = ud_2 - yd_2;
% end

vy(1) = 0;
vy_2(1) = 0;
vy_bound(1) = 0;
vy_subtractive(1) = 0;
vy_approx(1) = 0.0043;
for k = 1:T/Ts
    vy(k+1) = var(y(k,:)-yc(k));
%     vy_2(k+1) = var(y_2(k,:)-yc(k));
    vy_bound(k+1) = vy_bound(k) + (C1*(Ad1+Bd1*C1)^(k-1)*Bd1)^2 * d^2/4;
    vy_approx(k+1) = 0.0043;
    vy_subtractive(k+1) = vy_subtractive(k) + (C1*(Ad1+Bd1*C1)^(k-1)*Bd1)^2 * d^2/12;
end

% for k = 1:T/Ts
%     for l = 1:T/Ts
%         E_etaeta(k,l) = mean(times(eta(k,:),eta(l,:)));
%     end
% end

figure
hold on
plot(vy,'Color', [0,0,0], 'LineWidth', 2)
plot(vy_bound,'Color', [1,0,0], 'LineWidth', 2)
plot(vy_approx, 'Color', [0,0,1], 'LineWidth', 2)
% plot(vy_2, 'Color', [0,0,0], 'LineWidth', 2)
% plot(vy_subtractive, 'Color', [1,0,0], 'LineWidth', 2)
xlim([0,T/Ts]);
hold off

% k = 0:Ts:T-Ts;
% figure
% hold on
% plot(k,yc, 'Color', [0.5,0.5,0.5], 'LineWidth', 2);
% % ylim([-2.2,1.1]);
% % ylim([-2.2,1.1]);
% % stairs(k,yq, 'Color', [1,0,0], 'LineWidth', 1);
% % stairs(k,y(:,:),'LineWidth', 2);
% % stairs(k,y_2(:,:),'LineWidth', 2);
% hold off