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
[A_20,B_20,xu_20] = makeFD2(20,alpha,ubegin,uend);
[A_40,B_40,xu_40] = makeFD2(40,alpha,ubegin,uend);
[A_80,B_80,xu_80] = makeFD2(80,alpha,ubegin,uend);
[A_160,B_160,xu_160] = makeFD2(160,alpha,ubegin,uend);

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
title('2nd-order explicit scheme (IL=20)','FontSize',20);

saveas(gcf,'Q3_20','epsc')
hold off;

figure(2)
plot(xexact,uexact,'r-','LineWidth',LineWidth);
hold on;
plot(xu_40,U_40,'b--','LineWidth',LineWidth);
legend('exact','IL=40','Location','best');
xlabel("x",'FontSize',20);
ylabel("u",'FontSize',20);
title('2nd-order explicit scheme (IL=40)','FontSize',20);

saveas(gcf,'Q3_40','epsc')
hold off;

figure(3)
plot(xexact,uexact,'r-','LineWidth',LineWidth);
hold on;
plot(xu_80,U_80,'b--','LineWidth',LineWidth);
legend('exact','IL=80','Location','best');
xlabel("x",'FontSize',20);
ylabel("u",'FontSize',20);
title('2nd-order explicit scheme (IL=80)','FontSize',20);

saveas(gcf,'Q3_80','epsc')
hold off;

 function [A,B,xu] = makeFD2(IL,alpha,ubegin,uend)
 
    xu = linspace(0,1,IL+1);
    
    dx=1/IL;
    real_A=[1;-2+alpha^2.*dx^(2);1];

    
    
    A=zeros(IL+1,IL+1);

    B=zeros(IL+1,1);

    

    for i = 2:IL
        A(i,:)=[zeros(i-2,1);real_A;zeros(IL-i,1)];
    end
    
    A(1,1)=1;
    A(end,end)=1;

    
    B=transpose(linspace(0,1,IL+1).*2).*dx^(2);
    
    B(1)=ubegin;
    B(end)=uend;
    
    
    
    
    
 end