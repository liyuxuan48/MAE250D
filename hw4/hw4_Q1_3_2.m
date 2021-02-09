clc;
clear all;
%% duct, stretching stretching
%% grid construction

IL=61;
IT=41;
IE=21;
JL=41;


xend=1;
yend=2;
dx=xend/(IL-1);
dy=yend/(JL-1);


%% get alpha

alpha_number=3.13771;


% %% Construct airfoil 
% t=0.12;
% xET=linspace(0,1,(IT-IE)+1);
% yET=abs((t/0.2)*(0.2969*sqrt(xET) - 0.1260*xET - 0.3537*xET.^2 + 0.2843*xET.^3 - 0.1015*xET.^4));
% 
% 
% yL=zeros(IL,1);
% yL(IE:IT)=yET;

%% Construct airfoil 

xIL=linspace(0,xend,IL);
yT=(7+cos(4*pi()*xIL))/4;
yL=(1-cos(2*pi()*xIL))/4;


%% Construct the grid and boundaries of xy domain
xgrid=zeros(JL,IL);
ygrid=zeros(size(xgrid));

ygrid(1,:)=yL;
ygrid(end,:)=yT;
% ygrid(end,:)=yend;


xgrid(:,1)=0;
xgrid(:,end)=xend;


%% Initualize xygrid by interpolation

for i = 1:IL
        x_temp = linspace(0,xend,IL);
        xgrid(:,i)=x_temp(i);
end

for i = 1:IL
        yinter = (linspace(ygrid(1,i),ygrid(end,i),JL)-ygrid(1,i))./(ygrid(end,i)-ygrid(1,i));
        y_temp = (ygrid(end,i)-ygrid(1,i))*(exp(yinter.*alpha_number)-1)./(exp(alpha_number)-1)+ygrid(1,i);
        ygrid(:,i)=y_temp;
        
%         yinter = (linspace(ygrid(1,i),ygrid(end,i),JL)-ygrid(1,i))./(ygrid(end,i)-ygrid(1,i));
%         y_temp = (ygrid(end,i)-ygrid(1,i))*(yinter.*alpha_number-1)./(alpha_number-1)+ygrid(1,i);
%         ygrid(:,i)=y_temp;
end

% %% test initial grids
% 
% zgrid = zeros(size(xgrid))
% mesh(xgrid,ygrid,zgrid)


%% start iteration
error = 1;
counter = 1;
while error>1e-5

%% Get partial derivatives for x
x_s=zeros(size(xgrid));
x_n=zeros(size(xgrid));
x_ss=zeros(size(xgrid));
x_sn=zeros(size(xgrid));
x_nn=zeros(size(xgrid));
for i = 2:IL-1
    for j = 2:JL-1
        x_s(j,i)=(xgrid(j,i+1)-xgrid(j,i-1))/2;
        x_n(j,i)=(xgrid(j+1,i)-xgrid(j-1,i))/2;
%         x_ss(j,i)=(xgrid(j,i+1)-2*xgrid(j,i)+xgrid(j,i-1))/1;
%         x_sn(j,i)=(xgrid(j+1,i+1)-xgrid(j+1,i-1))/4 - (xgrid(j-1,i+1)-xgrid(j-1,i-1))/4;
%         x_nn(j,i)=(xgrid(j+1,i)-2*xgrid(j,i)+xgrid(j-1,i))/1;
    end
end

%% Get partial derivatives for y
y_s=zeros(size(ygrid));
y_n=zeros(size(ygrid));
y_ss=zeros(size(ygrid));
y_sn=zeros(size(ygrid));
y_nn=zeros(size(ygrid));
for i = 2:IL-1
    for j = 2:JL-1
        y_s(j,i)=(ygrid(j,i+1)-ygrid(j,i-1))/2;
        y_n(j,i)=(ygrid(j+1,i)-ygrid(j-1,i))/2;
%         y_ss(j,i)=(ygrid(j,i+1)-2*ygrid(j,i)+ygrid(j,i-1))/1;
%         y_sn(j,i)=(ygrid(j+1,i+1)-ygrid(j+1,i-1))/4 - (ygrid(j-1,i+1)-ygrid(j-1,i-1))/4;
%         y_nn(j,i)=(ygrid(j+1,i)-2*ygrid(j,i)+ygrid(j-1,i))/1;
    end
end

alpha=(x_n.^2)+(y_n.^2);
beta=x_n.*x_s+y_n.*y_s;
gamma=(x_s.^2)+(y_s.^2);


new_xgrid=xgrid;
new_ygrid=ygrid;
for i = 2:IL-1
    for j = 2:JL-1
        new_xgrid(j,i)=(0-alpha(j,i).*(xgrid(j,i+1)+xgrid(j,i-1)) + 2.*beta(j,i).*(xgrid(j+1,i+1)-xgrid(j-1,i+1)-xgrid(j+1,i-1)+xgrid(j-1,i-1))/4 - gamma(j,i).*(xgrid(j+1,i)+xgrid(j-1,i)))/(-2.*(alpha(j,i)+gamma(j,i)));
        new_ygrid(j,i)=(0-alpha(j,i).*(ygrid(j,i+1)+ygrid(j,i-1)) + 2.*beta(j,i).*(ygrid(j+1,i+1)-ygrid(j-1,i+1)-ygrid(j+1,i-1)+ygrid(j-1,i-1))/4 - gamma(j,i).*(ygrid(j+1,i)+ygrid(j-1,i)))/(-2.*(alpha(j,i)+gamma(j,i)));
    end
end


error_x_grid=new_xgrid-xgrid;
error_y_grid=new_ygrid-ygrid;

error=max(max(error_x_grid(:)),max(error_y_grid(:)));

xgrid=new_xgrid;
ygrid=new_ygrid;

counter= counter + 1;

end

zgrid = zeros(size(xgrid));

mesh(xgrid,ygrid,zgrid,'edgecolor', 'k');
view(2);
xlabel('x');
ylabel('y');
title('exponential grid stretching, duct','FontSize',20);
saveas(gcf,'exp_duct','epsc')

