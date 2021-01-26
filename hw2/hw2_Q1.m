% question 1 matrix stability 

clc;
clear all;

% alpha = 9*pi()/2;
IL = 100;
nulinspace = linspace(-1.2,1.2,201);
ubegin = 0;
% uend = 2;

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
   
 
