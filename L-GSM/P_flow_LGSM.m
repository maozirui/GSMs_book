% This code is writen to solve the 2D Couette flow 
% using Lagrangian Gradient Smoothing Method. For detailed informationl,
% please refer to "Zirui Mao, G.R. Liu, Int J Numer Methods Eng. 2018;113:858–890.".
clear all; clc; clf; 
%%%%%%%%%% controlling parameters %%%%%%%%
Lx=0.5e-3; Ly=1e-3; % computational domain
Nx=21; Ny=41; % amount of nodes on x, y edge of the domain
ds=Lx/(Nx-1); % grid size
rho0=1e3; % initial density
mu=1e-6; % kinetic viscosity
U0 = 2.5e-5; % Horizontal velocity of top lid
c = 1e-4; % sound speed in fluid 
TimeStep=20000; % total time steps to run
dt=1e-4; % critical time step
writeInterval=200; % time-step interval for results plot
%%%%%%%%%% generate intial particles %%%%%%%%%
X=[0:ds:Lx]; X=[X X]; Y=[zeros(1,Nx) Ly*ones(1,Nx)]; % generate no-slip boundary particles
Nb=length(X); % amount of no-slip boundary particles
x=zeros(Ny-2,Nx);   y=zeros(Ny-2,Nx); % interior nodes
for i=1:Ny-2
    for j=1:Nx
         x(i,j)=(j-1)*ds; y(i,j)=(i)*ds; % initially uniform interior nodes 
    end
end
X=[X reshape(x,[1,(Ny-2)*(Nx)])]; X=X-Lx/2; Y=[Y reshape(y,[1,(Ny-2)*(Nx)])]; % reshape
%%%%%%%%%%%%%%  ICs  %%%%%%%%%%
N=length(X); vx=zeros(1,N); vx(Nb/2+1:Nb)=U0; vy=zeros(1,N); rho=rho0*ones(1,N); pres=c*c*rho; 
drho=zeros(1,N); dvelocity_x=zeros(1,N); dvelocity_y=zeros(1,N);
rho_min=rho; vx_min=vx; vy_min=vy;  % store the mid-step variable by using leaf-frog time scheme
%%%%%%%%%%%%% update of field variable and locations %%%%%%%%%%%
time=0; s=0;  t=0;  
for itimestep=1:TimeStep
    t=t+dt;
    if itimestep~=1
        rho_min=rho; vx_min=vx; vy_min=vy; 
        rho=rho+dt/2*drho; pres=c*c*rho; vx=vx+dt/2*dvelocity_x; vy=vy+dt/2*dvelocity_y;
        vx(1:Nb/2)=0.0; vy(1:Nb/2)=0.0; vx(Nb/2+1:Nb)=U0; vy(Nb/2+1:Nb)=0; 
    end
    XX=[X X(1,1:2*Nx)]; YY=[Y -ds*ones(1,Nx) (Ly+ds)*ones(1,Nx)]; % add a layer of virtual particles for no-slip BCs
    vvx=[vx vx(1,1:2*Nx)]; vvy=[vy vy(1,1:2*Nx)]; rrho=[rho rho0*ones(1,Nx*2)]; 
    index_virt_L=find(XX>Lx/2-1.5*ds); index_virt_R=find(XX<1.5*ds-Lx/2); % find particles applying Peridical BCs
    XX=[XX XX(index_virt_L)-Lx-ds XX(index_virt_R)+Lx+ds]; % Periodical BC assignment
    YY=[YY YY(index_virt_L) YY(index_virt_R)];
    vvx=[vvx vvx(index_virt_L) vvx(index_virt_R)];  vvy=[vvy vvy(index_virt_L) vvy(index_virt_R)];
    rrho=[rrho rrho(index_virt_L) rrho(index_virt_R)];   ppres=c*c*rrho; % equation of state
    [drrho]=rho_partial (XX,YY,length(XX),ds, vvx, vvy, rrho); % calculate partial rho in continuity equation
    [dvelocity_px,dvelocity_py]=velocity_partial (XX,YY,length(XX),ds,ppres,rrho); % in momentum equation
    [avx, avy] = art_visc (XX, YY, length(XX), ds, vvx, vvy, rrho, mu); % calculate viscosity-acceleration term
    drho=drrho(1:N);% = right hand side of continuity equation
    dvelocity_x=-(dvelocity_px(1:N)+avx(1:N)); % = right hand side of momentum equation in x
    dvelocity_y=-(dvelocity_py(1:N)+avy(1:N)); % = right hand side of momentum equation in x
    if itimestep==1  % update field variables and positions
        rho=rho+dt/2*drho; pres=c*c*rho; vx=vx+dt/2*dvelocity_x;  vy=vy+dt/2*dvelocity_y;
        vx(1:Nb/2)=0.0; vy(1:Nb/2)=0.0; vx(Nb/2+1:Nb)=U0; vy(Nb/2+1:Nb)=0; 
        X=X+dt*vx; Y=Y+dt*vy; 
        X(X>Lx/2+ds/2)=X(X>Lx/2+ds/2)-Lx-ds; % periodical BC treatment
        X(X<-Lx/2-ds/2)=X(X<-Lx/2-ds/2)+Lx+ds; % periodical BC treatment
    else % update field variables and positions of particles
        rho=rho_min+dt*drho; pres=c*c*rho; vx=vx_min+dt*dvelocity_x;  vy=vy_min+dt*dvelocity_y;
        vx(1:Nb/2)=0.0; vy(1:Nb/2)=0.0; vx(Nb/2+1:Nb)=U0; vy(Nb/2+1:Nb)=0; 
        X=X+dt*vx; Y=Y+dt*vy; 
        X(X>Lx/2+ds/2)=X(X>Lx/2+ds/2)-Lx-ds; % periodical BC treatment
        X(X<-Lx/2-ds/2)=X(X<-Lx/2-ds/2)+Lx+ds; % periodical BC treatment
    end
    if itimestep/writeInterval == round(itimestep/writeInterval) % plot results
        s=s+1;  times(s)=t; 
        fprintf('Calculating ... time = %0.2f | percentage = %0.2f \n', t, itimestep/TimeStep*100.0)    
        subplot(1,2,1) %%%%%%%%% plot velocity contour
        scatter(X,Y,[],vx,'filled'); hold on;
        quiver(X,Y,vx,vy); hold off;
        view(2); axis equal; colorbar;  colormap jet;
        axis([-3.0e-4 3.0e-4 -0.5e-4 1.05e-3]);
        title('U'); xlabel('x');  ylabel('y')
        set(gca, 'Fontname', 'Times New Roman','FontSize',18);
        subplot(1,2,2) %%%%%%%%% compare numerical solution to exact solution
        vx_plot=[vx(1) vx(Nb+1:Nb+Ny-2) vx(Nx+1)];
        y_plot=[Y(1) Y(Nb+1:Nb+Ny-2) Y(Nx+1)];
        vx_theory_infinity=U0/Ly*y_plot; vx_theory=vx_theory_infinity;
        for n=1:10
            vx_theory=vx_theory+2*U0/n/pi*(-1)^n*sin(n*pi*y_plot/Ly)*exp(-mu*n^2*pi^2*t/Ly^2);
        end
        plot(vx_plot,y_plot,'o',vx_theory,y_plot,'-',vx_theory_infinity,y_plot,'k','linewidth',1.5);
        legend('L-GSM solution','Exact instant solution', 'Ux(t=infinity) in theory','location','southeast'); xlabel('U'); ylabel('y');
        set(gca, 'Fontname', 'Times New Roman','FontSize',18);
        drawnow;
        hh=suptitle(sprintf('Time =%0.3f sec ',t));
    end
end