clc;
clear all;

alpha = pi();
IL = 100;
c = 2;
ubegin = 0;
nu = 0.8;
% nu = 1.2; %1.2


% xexact_20 = linspace(0,1,21);
% uexact_20 = cos(alpha.*xexact_20)+(2-2./(alpha^2)-cos(alpha))./(sin(alpha)).*sin(alpha.*xexact_20) + 2*xexact_20./(alpha^2);
% 
% xexact_40 = linspace(0,1,41);
% uexact_40 = cos(alpha.*xexact_40)+(2-2./(alpha^2)-cos(alpha))./(sin(alpha)).*sin(alpha.*xexact_40) + 2*xexact_40./(alpha^2);
% 
% xexact_80 = linspace(0,1,81);
% uexact_80 = cos(alpha.*xexact_80)+(2-2./(alpha^2)-cos(alpha))./(sin(alpha)).*sin(alpha.*xexact_80) + 2*xexact_80./(alpha^2);
% 
% xexact_160 = linspace(0,1,161);
% uexact_160 = cos(alpha.*xexact_160)+(2-2./(alpha^2)-cos(alpha))./(sin(alpha)).*sin(alpha.*xexact_160) + 2*xexact_160./(alpha^2);
% 
% xexact = linspace(0,1,10001);
% uexact = cos(alpha.*xexact)+(2-2./(alpha^2)-cos(alpha))./(sin(alpha)).*sin(alpha.*xexact) + 2*xexact./(alpha^2);
% % 
% [A_20,B_20,xu_20] = makeFD2(20,alpha,ubegin,uend);
% [A_40,B_40,xu_40] = makeFD2(40,alpha,ubegin,uend);
% [A_80,B_80,xu_80] = makeFD2(80,alpha,ubegin,uend);
% [A_160,B_160,xu_160] = makeFD2(160,alpha,ubegin,uend);
% 
% U_20 = linsolve(A_20,B_20);
% U_40 = linsolve(A_40,B_40);
% U_80 = linsolve(A_80,B_80);
% U_160 = linsolve(A_160,B_160);
% % 
% error4_20 = max(abs(U_20-transpose(uexact_20)));
% error4_40 = max(abs(U_40-transpose(uexact_40)));
% error4_80 = max(abs(U_80-transpose(uexact_80)));
% error4_160 = max(abs(U_160-transpose(uexact_160)));
% 
% ratio4_4020 = error4_20/error4_40;
% ratio4_8040 = error4_40/error4_80;
% ratio4_16080 = error4_80/error4_160;

LineWidth=1.5;

% figure(1)
% plot(xexact,uexact,'r-','LineWidth',LineWidth);
% hold on;
% plot(xu_20,U_20,'b--','LineWidth',LineWidth);
% legend('exact','IL=40','Location','best');
% xlabel("x",'FontSize',20);
% ylabel("u",'FontSize',20);
% title('2nd-order explicit scheme (IL=20)','FontSize',20);
% 
% saveas(gcf,'Q3_20','epsc')
% hold off;
% 
% figure(2)
% plot(xexact,uexact,'r-','LineWidth',LineWidth);
% hold on;
% plot(xu_40,U_40,'b--','LineWidth',LineWidth);
% legend('exact','IL=40','Location','best');
% xlabel("x",'FontSize',20);
% ylabel("u",'FontSize',20);
% title('2nd-order explicit scheme (IL=40)','FontSize',20);
% 
% saveas(gcf,'Q3_40','epsc')
% hold off;
% 
% figure(3)
% plot(xexact,uexact,'r-','LineWidth',LineWidth);
% hold on;
% plot(xu_80,U_80,'b--','LineWidth',LineWidth);
% legend('exact','IL=80','Location','best');
% xlabel("x",'FontSize',20);
% ylabel("u",'FontSize',20);
% title('2nd-order explicit scheme (IL=80)','FontSize',20);
% 
% saveas(gcf,'Q3_80','epsc')
% hold off;

%  function [A,B,xu] = makeFD2(IL,alpha,ubegin,uend)
 
    xu = linspace(-1,1.5,IL+1)';
    u_initial = zeros(size(xu));
    for i = 1:length(xu)
        if xu(i) <= 0
            u_initial(i) = 5*sin(alpha*(xu(i)+1));
        end
    end
    
    u = u_initial(2:end);

    
    
    dx=1/(IL+1);
    
    real_A=[1+nu;-nu];

    
    
    A=zeros(IL,IL);

    B=zeros(IL,1);

    

    for i = 1:IL-1
        A(i,:)=[zeros(i-1,1);real_A;zeros(IL-i-1,1)];
    end
    
%     A(1,1)=1-nu;

%     A(end,end-1:end)=[nu;1-nu]; % 1st order
    A(end,end-3:end)=[-1*nu^2/2;-1*nu/2+4*nu^2/2;4*nu/2-5*nu^2/2;1-3*nu/2+2*nu^2/2]; % 2nd order
    
    
    NI = 40;
    u = (A^NI)*u;
    
    u_theo = zeros(size(u_initial));
    
    for i = 1:length(xu)
        if xu(i) <= nu && xu(i) >= -1+nu
            u_theo(i) = 5*sin(alpha*(xu(i)-nu+1));
        end
    end
    
    figure(1)
plot(xu(2:end),u_theo(2:end),'r-','LineWidth',1.5);
hold on;
plot(xu(2:end),u,'b-','LineWidth',1.5);
hold on;
legend('$u_{exact}$','$u$','Interpreter','latex','Location','best','FontSize',18);
xlabel("x",'Interpreter','latex','FontSize',20);
ylabel("u",'Interpreter','latex','FontSize',20);
title('Forward explicit $\nu$=0.8','Interpreter','latex','FontSize',20);
saveas(gcf,'Q2_forward_08','epsc')
% title('Forward explicit $\nu$=1.2','Interpreter','latex','FontSize',20);
% saveas(gcf,'Q2_forward_12','epsc')


%     max(abs(real(eig(A))))
%  end