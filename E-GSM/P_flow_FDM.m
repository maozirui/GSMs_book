% This code is writen to solve the 2D Poiseuelle flow 
% using Finite Difference Scheme.
% Author: Zirui Mao (maozirui0827@gmail.com)
% Last Updated: Oct., 2020
clear all; clc; clf;
%%%%%%%%%% controlling parameters %%%%%%%%
N=101; % amount of nodes on one edge of the square domain
ds=1/(N-1); % grid size
P=1; % -dpdx pressure gradient term
mu=1; nu = mu/1; % mu <=1 is dynamic viscosity of fluid; rho = 1
dt=ds^2*0.5; % time interval
writeInterval=1000; % time-step interval for plotting results 
criterion_u=1.0e-7; % convergence criterion of velocity calculation
%%%%%%%%%% define domain grid %%%%%%%%%
x=zeros(N,N);   y=zeros(N,N);
for i=1:N
    for j=1:N
         x(i,j)=(j-1)/(N-1); y(i,j)=(i-1)/(N-1)-0.5; 
    end
end
%%%%%%%%%%%%%%  ICs and BCs   %%%%%%%%%%
u=zeros(N,N); % initial condition settings
u(1,:)=0.0; u(N,:)=0.0; % no-slip boundary condition
time=0; s=0;  t=0;  iter=0;    Res_u=1.0;
%%%%%%%%%%%%% velocity update %%%%%%%%%%%
while (Res_u>criterion_u)
    uold=u;   iter=iter+1;   t=t+dt;
    u_yy=zeros(N,N); % second derivative of u w.r.t. y
    u_n=u(1:N-2,:); u_s=u(3:N,:); % the upper and lower grid
    u_yy=(u_n+u_s-2.0*u(2:N-1,:))/ds/ds; % computate u_xx with central FD scheme
    u(2:N-1,:)=u(2:N-1,:)+dt*(P+u_yy*nu); % update velocity based governing equation
    u(1,:)=0.0; u(N,:)=0.0; % no-slip boundary condition on top and bottom
    if iter/writeInterval == round(iter/writeInterval)
        s=s+1;  times(s)=t; 
        res_u(s)=max(max(abs(u-uold))); % maximum residual of increament of velocity
        Res_u=res_u(s);
        fprintf('Calculating ... time = %0.2f | Res_u-7=%0.2f \n', t, log10(Res_u))    
        hh=suptitle(sprintf('Res_u-7=%0.2f | Time =%0.3f sec ',log10(Res_u),t));
        subplot(1,2,1) %%%%%%%%% plot velocity contour
        h=surf(x,y,u, 'facecolor', 'interp');
        view(2); axis equal; axis([0 1 -0.5 0.5])
        title('U'); xlabel('x');  ylabel('y')
        colorbar;  colormap jet; caxis([0 P/2/mu*0.5^2]);
        set(gca, 'Fontname', 'Times New Roman','FontSize',15);
        subplot(1,2,2)  %%%%%%%%% compare numerical solution to exact solution
        plot(u(:,1),y(:,1),'o',P/2/mu*(0.5^2-y(:,1).^2),y(:,1),'-','linewidth',1.5);
        legend('FDM solution','Exact solution'); xlabel('U'); ylabel('y');
        drawnow;
    end
end