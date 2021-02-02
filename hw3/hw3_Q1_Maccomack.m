% MacCormack scheme,
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
    
    %% predictor
    real_A_predict=[1+nu;-nu];
    A_predict=zeros(IL,IL);
    for i = 1:IL-1
        A_predict(i,:)=[zeros(i-1,1);real_A_predict;zeros(IL-i-1,1)];
    end
%     A(1,1:2)=[1;nu/2];
%     A_predict(end,end-1:end)=[nu;1-nu]; % 1st order
%     A(end,end-3:end)=[-1*nu^2/2;-1*nu/2+4*nu^2/2;4*nu/2-5*nu^2/2;1-3*nu/2+2*nu^2/2]; % 2nd order
%     A(end,end-1:end)=[-nu;1+nu]; % 1st order for implicit
%      A(end,end-3:end)=[1*nu^2/2;+1*nu/2-4*nu^2/2;-4*nu/2+5*nu^2/2;1+3*nu/2-2*nu^2/2]; % 2nd order for implicit

    %% corrector
    real_A_correct=[nu/2;1/2-nu/2];
    A_correct=zeros(IL,IL);
    A2_correct=A_correct;
    Abc=A_correct;
    for i = 2:IL-1
        A_correct(i,:)=[zeros(i-2,1);real_A_correct;zeros(IL-i,1)];
        
        A2_correct(i,i)=1/2;
    end
    A_correct(1,1)=[1/2-nu/2];
%     A_correct(end,end-1:end)=[nu;1-nu]; % 1st order
%     A(end,end-3:end)=[-1*nu^2/2;-1*nu/2+4*nu^2/2;4*nu/2-5*nu^2/2;1-3*nu/2+2*nu^2/2]; % 2nd order
%     A(end,end-1:end)=[-nu;1+nu]; % 1st order for implicit
%      A(end,end-3:end)=[1*nu^2/2;+1*nu/2-4*nu^2/2;-4*nu/2+5*nu^2/2;1+3*nu/2-2*nu^2/2]; % 2nd order for implicit
    A2_correct(1,1)=1/2;
    Abc(end,end-1:end)=[nu;1-nu]; % 1st order
    
    
    NI = 40;
    
    A_final = (A_correct*A_predict+A2_correct+Abc);
    
    u = (A_final^NI)*u;
    
    
    
    u_theo = zeros(size(u_initial));
    
    for i = 1:length(xu)
        if xu(i) <= nu*100/IL && xu(i) >= -1+nu*100/IL
            u_theo(i) = 5*sin(alpha*(xu(i)-nu*100/IL+1));
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
title('MacCormack scheme $\nu$=0.8','Interpreter','latex','FontSize',20);
saveas(gcf,'Q2_MacCormack_08','epsc')
% title('MacCormack scheme $\nu$=1.2','Interpreter','latex','FontSize',20);
% saveas(gcf,'Q2_MacCormack_12','epsc')

%  end