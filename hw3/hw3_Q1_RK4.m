%RK4 scheme,
clc;
clear all;

alpha = pi();
IL = 100;
c = 2;
ubegin = 0;
% nu = 0.8;
nu = 1.2; %1.2


    xu = linspace(-1,1.5,IL+1)';
    u_initial = zeros(size(xu));
    for i = 1:length(xu)
        if xu(i) <= 0
            u_initial(i) = 5*sin(alpha*(xu(i)+1));
        end
    end
    
    
    dx=(xu(end)-xu(1))/(IL);
    dt = nu*dx/c;

    LineWidth=1.5;
    
    u = u_initial(2:end);

    
    
    real_A=[-1;0;1];
    A=zeros(IL,IL);

    
    for i = 2:IL-1
        A(i,:)=[zeros(i-2,1);real_A;zeros(IL-i-1,1)];
    end
    
    A(1,1:2)=[0;1];

    A(end,end-1:end)=[-2;2]; % 1st order right end
    
    NI = 40;
    
    for j = 1:NI
    
    k1 = -c/2/dx*A*u;
    
    k2 = -c/2/dx*A*(u + dt/2*k1);
    
    k3 = -c/2/dx*A*(u + dt/2*k2);
    
    k4 = -c/2/dx*A*(u + dt*k3);
    
    u = u + dt/6*(k1+2*k2+2*k3+k4);
    
    end

%     u = (A^NI)*u;
    
    u_theo = zeros(size(u_initial));
    
    for i = 1:length(xu)
        if xu(i) <= nu*NI/40 && xu(i) >= -1+nu*NI/40
            u_theo(i) = 5*sin(alpha*(xu(i)-nu*NI/40+1));
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
% title('4th-order RK $\nu$=0.8','Interpreter','latex','FontSize',20);
% saveas(gcf,'4th_order_RK_08','epsc')
title('4th-order RK $\nu$=1.2','Interpreter','latex','FontSize',20);
saveas(gcf,'4th_order_RK_12','epsc')
% 
% %  end