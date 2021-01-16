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

[AFD4_20,BFD4_20,xu_20] = makeAFD4(20,alpha,ubegin,uend);
[AFD4_40,BFD4_40,xu_40] = makeAFD4(40,alpha,ubegin,uend);
[AFD4_80,BFD4_80,xu_80] = makeAFD4(80,alpha,ubegin,uend);
[AFD4_160,BFD4_160,xu_160] = makeAFD4(160,alpha,ubegin,uend);

U_20 = linsolve(AFD4_20,BFD4_20);
U_40 = linsolve(AFD4_40,BFD4_40);
U_80 = linsolve(AFD4_80,BFD4_80);
U_160 = linsolve(AFD4_160,BFD4_160);

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
title('4th-order explicit scheme (IL=20)','FontSize',20);

saveas(gcf,'Q1_20','epsc')
hold off;

figure(2)
plot(xexact,uexact,'r-','LineWidth',LineWidth);
hold on;
plot(xu_40,U_40,'b--','LineWidth',LineWidth);
legend('exact','IL=40','Location','best');
xlabel("x",'FontSize',20);
ylabel("u",'FontSize',20);
title('4th-order explicit scheme (IL=40)','FontSize',20);

saveas(gcf,'Q1_40','epsc')
hold off;

figure(3)
plot(xexact,uexact,'r-','LineWidth',LineWidth);
hold on;
plot(xu_80,U_80,'b--','LineWidth',LineWidth);
legend('exact','IL=80','Location','best');
xlabel("x",'FontSize',20);
ylabel("u",'FontSize',20);
title('4th-order explicit scheme (IL=80)','FontSize',20);

saveas(gcf,'Q1_80','epsc')
hold off;

 function [AFD4,BFD4,xu] = makeAFD4(IL,alpha,ubegin,uend)
 
    xu = linspace(0,1,IL+1);
    
    A=sym([1 1 1 1 1;-1 0 1 2 3;1 0 1 4 9; -1 0 1 8 27;1 0 1 16 81]);
    B=sym([0 0 2 0 0]') ;
 
    FD4=[-1;16;-30;16;-1];
    fFD3=double([linsolve(A,B)].*12);
    bFD3=flip(fFD3);

	dx=1/IL;
    
    real_fFD3=fFD3+[0;12*alpha^(2)*dx^(2);0;0;0];
    real_bFD3=flip(real_fFD3);
    real_FD4=FD4+[0;0;12*alpha^(2)*dx^(2);0;0];
    
    AFD4=zeros(IL+1,IL+1);
    BFD4=zeros(IL+1,1);
    
    AFD4(1,1)=1;
    AFD4(end,end)=1;
    AFD4(2,:)=[real_fFD3;zeros(IL-4,1)];
    AFD4(end-1,:)=[zeros(IL-4,1);real_bFD3];
    
    for i = 3:IL-1
        AFD4(i,:)=[zeros(i-3,1);real_FD4;zeros(IL-i-1,1)];
    end
    
    center_BFD4=transpose(12*dx^(2)*2*[1:IL-1]./IL);
    BFD4=[ubegin;center_BFD4;uend];
 end
 
 
