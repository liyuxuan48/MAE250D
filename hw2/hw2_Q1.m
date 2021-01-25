clc;
clear all;

% alpha = 9*pi()/2;
IL = 100;
nulinspace = linspace(-1.2,1.2,200);
ubegin = 0;
% uend = 2;


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
% 
% [AFD4_20,BFD4_20,xu_20] = makeAFD4(20,alpha,ubegin,uend);
% [AFD4_40,BFD4_40,xu_40] = makeAFD4(40,alpha,ubegin,uend);
% [AFD4_80,BFD4_80,xu_80] = makeAFD4(80,alpha,ubegin,uend);
% [AFD4_160,BFD4_160,xu_160] = makeAFD4(160,alpha,ubegin,uend);
% 
% U_20 = linsolve(AFD4_20,BFD4_20);
% U_40 = linsolve(AFD4_40,BFD4_40);
% U_80 = linsolve(AFD4_80,BFD4_80);
% U_160 = linsolve(AFD4_160,BFD4_160);
% 
% error4_20 = max(abs(U_20-transpose(uexact_20)));
% error4_40 = max(abs(U_40-transpose(uexact_40)));
% error4_80 = max(abs(U_80-transpose(uexact_80)));
% error4_160 = max(abs(U_160-transpose(uexact_160)));
% 
% ratio4_4020 = error4_20/error4_40;
% ratio4_8040 = error4_40/error4_80;
% ratio4_16080 = error4_80/error4_160;
% 
% LineWidth=1.5;
% 
% figure(1)
% plot(xexact,uexact,'r-','LineWidth',LineWidth);
% hold on;
% plot(xu_20,U_20,'b--','LineWidth',LineWidth);
% legend('exact','IL=20','Location','best');
% xlabel("x",'FontSize',20);
% ylabel("u",'FontSize',20);
% title('4th-order explicit scheme (IL=20)','FontSize',20);
% 
% saveas(gcf,'Q1_20','epsc')
% hold off;
% 
% figure(2)
% plot(xexact,uexact,'r-','LineWidth',LineWidth);
% hold on;
% plot(xu_40,U_40,'b--','LineWidth',LineWidth);
% legend('exact','IL=40','Location','best');
% xlabel("x",'FontSize',20);
% ylabel("u",'FontSize',20);
% title('4th-order explicit scheme (IL=40)','FontSize',20);
% 
% saveas(gcf,'Q1_40','epsc')
% hold off;
% 
% figure(3)
% plot(xexact,uexact,'r-','LineWidth',LineWidth);
% hold on;
% plot(xu_80,U_80,'b--','LineWidth',LineWidth);
% legend('exact','IL=80','Location','best');
% xlabel("x",'FontSize',20);
% ylabel("u",'FontSize',20);
% title('4th-order explicit scheme (IL=80)','FontSize',20);
% 
% saveas(gcf,'Q1_80','epsc')
% hold off;
% 
%  function [AFD4,BFD4,xu] = makeAFD4(IL,alpha,ubegin,uend)
%  

eigenlinspaceA1 = zeros(size(nulinspace));
eigenlinspaceA2 = zeros(size(nulinspace));


for i = 1:size(nulinspace,2)
    
    nu = nulinspace(i);
    
    xu = linspace(0,1,IL+1);
	dx=1/IL;
    
    slice_A=[nu/2+(nu^2)/2;1-nu^2;-nu/2+(nu^2)/2];
    
    A=zeros(IL,IL);
    B=zeros(IL,1);
    
    A(1,1:2)=[1-nu^2;-nu/2+(nu^2)/2];
    
    for j = 2:IL-1
        A(j,:)=[zeros(j-2,1);slice_A;zeros(IL-j-1,1)];
    end
    
     A(end,end-1:end)=[nu;1-nu]; % 1st order
     
     A1 = A;
     
     A2 = A;
     A2(end,end-3:end)=[-1*nu^2/2;-1*nu/2+4*nu^2/2;4*nu/2-5*nu^2/2;1-3*nu/2+2*nu^2/2]; % 2nd order
    
%     eigenlinspaceA1(i) = max(abs(real(eig(A1))));
%     eigenlinspaceA2(i) = max(abs(real(eig(A2))));
     eigenlinspaceA1(i) = max(abs((eig(A1))));
     eigenlinspaceA2(i) = max(abs((eig(A2))));
end

figure(1)
plot(nulinspace,eigenlinspaceA1,'r-','LineWidth',1.5);
hold on;
plot(nulinspace,eigenlinspaceA2,'b-','LineWidth',1.5);
hold on;
plot([-1.5,1.5],[1,1],'k--','LineWidth',2);
legend('First-order one-sided','Second-order one-sided','Location','best','FontSize',10);
xlabel("$\nu$",'Interpreter','latex','FontSize',20);
ylabel("$|\lambda|_{max}$",'Interpreter','latex','FontSize',20);
title('matrix stability','FontSize',20);

saveas(gcf,'Q1_1','epsc')
hold off;
    
%     
%     center_BFD4=transpose(12*dx^(2)*2*[1:IL-1]./IL);
%     BFD4=[ubegin;center_BFD4;uend];
%  end

 
