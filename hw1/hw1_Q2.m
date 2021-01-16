clc;
clear all;

alpha = 9*pi()/2;
ubegin = 1;
uend = 2;


xexact_20 = linspace(0,1,21);
uexact_20 = cos(alpha.*xexact_20)+(2-2./(alpha^2)-cos(alpha))./(sin(alpha)).*sin(alpha.*xexact_20) + 2*xexact_20./(alpha^2);

xexact_40 = linspace(0,1,41);
uexact_40 = cos(alpha.*xexact_40)+(2-2./(alpha^2)-cos(alpha))./(sin(alpha)).*sin(alpha.*xexact_40) + 2*xexact_40./(alpha^2);

xexact_80 = linspace(0,1,81);
uexact_80 = cos(alpha.*xexact_80)+(2-2./(alpha^2)-cos(alpha))./(sin(alpha)).*sin(alpha.*xexact_80) + 2*xexact_80./(alpha^2);

xexact_160 = linspace(0,1,161);
uexact_160 = cos(alpha.*xexact_160)+(2-2./(alpha^2)-cos(alpha))./(sin(alpha)).*sin(alpha.*xexact_160) + 2*xexact_160./(alpha^2);

xexact = linspace(0,1,10001);
uexact = cos(alpha.*xexact)+(2-2./(alpha^2)-cos(alpha))./(sin(alpha)).*sin(alpha.*xexact) + 2*xexact./(alpha^2);
% 
[A_20,B_20,xu_20] = makeFD4C(20,alpha,ubegin,uend);
[A_40,B_40,xu_40] = makeFD4C(40,alpha,ubegin,uend);
[A_80,B_80,xu_80] = makeFD4C(80,alpha,ubegin,uend);
[A_160,B_160,xu_160] = makeFD4C(160,alpha,ubegin,uend);

U_20 = linsolve(A_20,B_20);
U_40 = linsolve(A_40,B_40);
U_80 = linsolve(A_80,B_80);
U_160 = linsolve(A_160,B_160);
% 
error4_20 = max(abs(U_20-transpose(uexact_20)));
error4_40 = max(abs(U_40-transpose(uexact_40)));
error4_80 = max(abs(U_80-transpose(uexact_80)));
error4_160 = max(abs(U_160-transpose(uexact_160)));

ratio4_4020 = error4_20/error4_40;
ratio4_8040 = error4_40/error4_80;
ratio4_16080 = error4_80/error4_160;

LineWidth=1.5;

figure(1)
plot(xexact,uexact,'r-','LineWidth',LineWidth);
hold on;
plot(xu_20,U_20,'b--','LineWidth',LineWidth);
legend('exact','IL=20','Location','best');
xlabel("x",'FontSize',20);
ylabel("u",'FontSize',20);
title('4th-order compact scheme (IL=20)','FontSize',20);

saveas(gcf,'Q2_20','epsc')
hold off;

figure(2)
plot(xexact,uexact,'r-','LineWidth',LineWidth);
hold on;
plot(xu_40,U_40,'b--','LineWidth',LineWidth);
legend('exact','IL=40','Location','best');
xlabel("x",'FontSize',20);
ylabel("u",'FontSize',20);
title('4th-order compact scheme (IL=40)','FontSize',20);

saveas(gcf,'Q2_40','epsc')
hold off;

figure(3)
plot(xexact,uexact,'r-','LineWidth',LineWidth);
hold on;
plot(xu_80,U_80,'b--','LineWidth',LineWidth);
legend('exact','IL=80','Location','best');
xlabel("x",'FontSize',20);
ylabel("u",'FontSize',20);
title('4th-order compact scheme (IL=80)','FontSize',20);

saveas(gcf,'Q2_80','epsc')
hold off;

 function [A,B,xu] = makeFD4C(IL,alpha,ubegin,uend)
 
    xu = linspace(0,1,IL+1);
    
    dx=1/IL;
    real_P=[1;10;1].*dx^(2)./12;
    real_Q=[1;-2;1];

    
    
    P=zeros(IL+1,IL+1);
    Q=zeros(IL+1,IL+1);
    q=zeros(IL+1,1);
    f=zeros(IL+1,1);
    

    for i = 2:IL
        P(i,:)=[zeros(i-2,1);real_P;zeros(IL-i,1)];
    end
    
    Q(1,1)=-1;
    Q(end,end)=-1;
    for i = 2:IL
        Q(i,:)=[zeros(i-2,1);real_Q;zeros(IL-i,1)];
    end
    
    q=transpose(linspace(0,1,IL+1).*2);
    
    f(1)=ubegin;
    f(end)=uend;
    
    A = Q + alpha^2.*P;
    B = P*q - f;
    
 end