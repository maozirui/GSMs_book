function [dfx] = GSM_standard( N,x,y,o,angle,f,f_ref)
% This function is written to approximate gradient of scalar function f
% with the standard GSM gradient operator based on n-GSD. 
dfx=0;
%%%%%%%%%%[start] area calculation [start] %%%%%%%%%%%%%%
area=(x(N).*y(1)-x(1).*y(N))/6;
for j=1:N-1
    area=area+(x(j).*y(j+1)-x(j+1).*y(j))/6; % calculate the area V of n-GSD
end
k=0;
%%% evaluate the target particle is on boundary or not
for j=1:N
    jp1=j+1;
    if j==N
        jp1=1;
    end
    angle_up=angle(jp1)-angle(j);
    if angle_up<0
        angle_up=angle_up+360;
    end
    if angle_up>180
        k=j;
    end
end
if k~=0 % meaning the target particle locates on boundary
    X=x; Y=y;
    if k~=N
        j=k;
        for i=1:(N-j)
            X(i)=x(j+i); Y(i)=y(j+i);
        end
        for i=1:j
            X(N-j+i)=x(i); Y(N-j+i)=y(i);
        end
    end
    X=[X; o(1)]; Y=[Y; o(2)];
    area=(X(N+1).*Y(1)-X(1).*Y(N+1))/6; % recalcuate the area of n-GSD for boundary particle
    for j=1:N
        area=area+(X(j).*Y(j+1)-X(j+1).*Y(j))/6;
    end
end
%%%%%%%%%%%%%% gradient approximation %%%%%%%%%%%%%%%%
for j=1:N
    jm1=j-1;jp1=j+1;
    if j==1
        jm1=N;
    elseif j==N
        jp1=1;
    end
    angle_up=angle(jp1)-angle(j); % Left part
    if angle_up<0
        angle_up=angle_up+360;
    end
    angle_low=angle(j)-angle(jm1); % right part
    if angle_low<0
        angle_low=angle_low+360;
    end
    if angle_up<180
        x_up=(x(jp1)+x(j)+o(1))/3; y_up=(y(jp1)+y(j)+o(2))/3;
        dSx=y_up-(o(2)+y(j))/2;
        dfx=dfx+(f(j)+f_ref)/2.*dSx/area;
    end
    if angle_low<180
        x_low=(x(jm1)+x(j)+o(1))/3; y_low=(y(jm1)+y(j)+o(2))/3;
        dSx=-y_low+(o(2)+y(j))/2;
        dfx=dfx+(f(j)+f_ref)/2.*dSx/area;
    end
end
end

