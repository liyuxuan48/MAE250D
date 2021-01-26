%central explicit scheme,
clc;
clear all;

alpha = pi();
IL = 100;
c = 2;
ubegin = 0;
nu = 0.8;
% nu = 1.2; %1.2


LineWidth=1.5;


    xu = linspace(-1,1.5,IL+1)';
    u_initial = zeros(size(xu));
    for i = 1:length(xu)
        if xu(i) <= 0
            u_initial(i) = 5*sin(alpha*(xu(i)+1));
        end
    end
    
    u = u_initial(2:end);

    
    
    dx=1/(IL+1);
    
    real_A=[nu/2;1;-nu/2];

    
    
    A=zeros(IL,IL);

    B=zeros(IL,1);

    

    for i = 2:IL-1
        A(i,:)=[zeros(i-2,1);real_A;zeros(IL-i-1,1)];
    end
    
    A(1,1:2)=[1;-nu/2];

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
title('Central explicit $\nu$=0.8','Interpreter','latex','FontSize',20);
saveas(gcf,'Q2_central_08','epsc')
% title('Central explicit $\nu$=1.2','Interpreter','latex','FontSize',20);
% saveas(gcf,'Q2_central_12','epsc')

%  end