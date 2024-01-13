function [avx,avy] = art_visc (X, Y, NN, ds, vx, vy, rho, mu)
avx=zeros(1,NN); avy=zeros(1,NN);
for k=1:NN % for all real particles
    %%%%%%% find out the neighboring particles %%%%%%%%%%%
    distance=sqrt((X-X(k)).^2+(Y-Y(k)).^2); % distance from other particles
    L=find(distance<=1.5*ds); % supporting particles list
    L(L==k)=[];
    %%%%%%% ordering supporting particles anti-clockwisely %%%%%%%%%
    N=length(L);
    P_local=[X(L)-X(k); Y(L)-Y(k)];
    angle=zeros(1,N);
    for i=1:N
        if P_local(2,i)>=0
            angle(i)=acosd(P_local(1,i)/distance(L(i)));
        else
            angle(i)=360-acos(P_local(1,i)/distance(L(i)));
        end
    end
    [angle,I]=sort(angle);
    L=L(I);   
    x=X(L); y=Y(L);
    %%%%%% calculate area of n-GSD %%%%%%%%%%
    area=(x(N).*y(1)-x(1).*y(N))/6;
    o=[X(k); Y(k)];
    DIST=(o(1)-x).^2+(o(2)-y).^2;
    for j=1:N-1
        area=area+(x(j).*y(j+1)-x(j+1).*y(j))/6; % calculate the area V of n-GSD
    end
    %%%%%%%%%%%%%% gradient approximation %%%%%%%%%%%%%%%%
    for j=1:N
        jm1=j-1;jp1=j+1;
        if j==1
            jm1=N;
        elseif j==N
            jp1=1;
        end
        x_up=(x(jp1)+x(j)+o(1))/3; y_up=(y(jp1)+y(j)+o(2))/3;
        dSx=y_up-(o(2)+y(j))/2; dSy=x_up-(o(1)+x(j))/2;
        x_low=(x(jm1)+x(j)+o(1))/3; y_low=(y(jm1)+y(j)+o(2))/3;
        dSx=dSx-y_low+(o(2)+y(j))/2; dSy=dSy-x_low+(o(1)+x(j))/2;
        dfx=(o(1)-x(j))*dSx/area;
        dfy=-(o(2)-y(j))*dSy/area;
        %%%%%%%% calculate the viscosity-caused acceleration %%%%%%%%%%%
        avx(k)=avx(k)-mu*(rho(k)+rho(L(j)))/2/rho(k)/DIST(j)*(dfx+dfy)*(vx(k)-vx(L(j)));
        avy(k)=avy(k)-mu*(rho(k)+rho(L(j)))/2/rho(k)/DIST(j)*(dfx+dfy)*(vy(k)-vy(L(j)));
    end
end
end