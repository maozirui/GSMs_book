% This code is writen to solve the 2D Poiseuelle flow 
% using Gradient Smoothing Method.
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
X=[0:ds:1]; X=[X X]; Y=[zeros(1,N) ones(1,N)]; % generate top and bottom no-slip nodes
X=[X zeros(1,N-2) ones(1,N-2)]; % generate x coordinates of periodical boundary nodes
Y=[Y [ds:ds:1-ds] [ds:ds:1-ds]]; % generate y coordinates of periodical boundary nodes
x=zeros(N-2,N-2);   y=zeros(N-2,N-2); % irregularly distributed interior nodes
dx=0.4*(rand(N-2)-0.5)*2*ds; dy=0.4*(rand(N-2)-0.5)*2*ds; % noise assigned to uniform nodes
for i=2:N-1
    for j=2:N-1
         x(i-1,j-1)=(j-1)/(N-1); y(i-1,j-1)=(i-1)/(N-1); % uniform interior nodes 
    end
end
x=x+dx; y=y+dy; % generate irregularly distributed interior nodes
X=[X reshape(x,[1,(N-2)*(N-2)])]; Y=[Y reshape(y,[1,(N-2)*(N-2)])]; Y=Y-0.5; % reshape
index_L=[1 2*N+1:3*N-2 N+1]; % index array of left boundary nodes
index_R=[N 3*N-1:4*N-4 2*N]; % index array of right boundary nodes
index_L_virt=[N*N+1:N*N+N]; % index array of left virtual nodes
index_R_virt=[N*N+N+1:N*N+2*N]; % index array of right virtual nodes
X=[X X(index_L)+1.0+ds 1.0-X(index_R)-ds]; % assign virtual boundary nodes 
Y=[Y Y(index_L) Y(index_R)]; % assign virtual boundary nodes 
T=delaunay(X,Y); %%%%%%%%% generate unstructured mesh via Delaunay algorithm
%%%%%%%%%%%%%%  ICs and BCs   %%%%%%%%%%
u=zeros(length(X),1); % initial condition of velocity u 
u(1:2*N)=0.0; % no-slip BCs for top and bottom boundary
u(index_L_virt)=u(index_L); u(index_R_virt)=u(index_R); % Periodic left and right BCs
%%%%%%%%%%%%% velocity update %%%%%%%%%%%
time=0; s=0;  t=0;  iter=0;    Res_u=1.0;   
while (Res_u>criterion_u)
    uold=u;   iter=iter+1;   t=t+dt;
    [u_x,u_y]=GSM_gradient(T,X,Y,length(X),u); % calculate velocity gradient
    u_y(index_L_virt)=u_y(index_L); u_y(index_R_virt)=u_y(index_R); % Periodic BCs
    [u_xy,u_yy]=GSM_gradient(T,X,Y,length(X),u_y); % calculate velocity u_yy
    u_yy(index_L_virt)=u_yy(index_L); u_yy(index_R_virt)=u_yy(index_R); % Periodic BCs
    u(2*N+1:N*N)=u(2*N+1:N*N)+dt*(P+u_yy(2*N+1:N*N)*nu); % update u of interior domain
    u(1:2*N)=0.0; % no-slip BCs for top and bottom boundary
    u(index_L_virt)=u(index_L); u(index_R_virt)=u(index_R); % Periodic left and right BCs
    if iter/writeInterval == round(iter/writeInterval)
        s=s+1;  times(s)=t; 
        res_u(s)=max(max(abs(u-uold)));  % maximum residual of increament of velocity
        Res_u=res_u(s);
        fprintf('Calculating ... time = %0.2f | Res_u-7=%0.2f \n', t, log10(Res_u))    
        hh=suptitle(sprintf('Res_u-7=%0.2f | Time =%0.3f sec ',log10(Res_u),t));
        subplot(1,2,1) %%%%%%%%% plot velocity contour
        h=trimesh(T,X,Y,u,'FaceColor','interp','EdgeColor','k');
        view(2); axis equal; axis([0 1 -0.5 0.5])
        title('U'); xlabel('x');  ylabel('y')
        colorbar;  colormap jet; caxis([0 P/2/mu*0.5^2]);
        set(gca, 'Fontname', 'Times New Roman','FontSize',15);
        subplot(1,2,2) %%%%%%%%% compare numerical solution to exact solution
        plot(u(index_L),Y(index_L),'o',P/2/mu*(0.5^2-Y(index_L).^2),Y(index_L),'-','linewidth',1.5);
        legend('GSM solution','Exact solution'); xlabel('U'); ylabel('y');
        drawnow;
    end
end