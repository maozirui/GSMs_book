function [dfx] = GSM_normalized( N,x,y,o,angle,f,f_ref)
% This function is written to approximate gradient of scalar function f
% with the normalized GSM gradient operator. For detailed information, 
% please refer to "Zirui Mao., et al.,Int J Numer Methods Eng. 2018;113:858–890."
sigma1=0; sigma2=0; sigma3=0; sigma4=0; A=0; B=0;  
for j=1:N
    jm1=j-1;jp1=j+1; % jm1 = last supporting particle; jp1 = next particle
    if j==1
        jm1=N;
    elseif j==N
        jp1=1;
    end
    angle_up=angle(jp1)-angle(j); % angle from next supporting particle
    if angle_up<0
        angle_up=angle_up+360;
    end
    angle_low=angle(j)-angle(jm1); % angle from last supporting particle
    if angle_low<0
        angle_low=angle_low+360;
    end
    if angle_up>180 % boundary particle
        x_up=(o(1)+x(j))/2; y_up=(o(2)+y(j))/2;
    else
        x_up=(x(jp1)+x(j)+o(1))/3;  y_up=(y(jp1)+y(j)+o(2))/3;
    end
    if angle_low>180 % boundary particle
        x_low=(o(1)+x(j))/2 ;  y_low=(o(2)+y(j))/2;
    else
        x_low=(x(jm1)+x(j)+o(1))/3;  y_low=(y(jm1)+y(j)+o(2))/3;
    end
    dSx=y_up-y_low; dSy=x_up-x_low;
    sigma1=sigma1+dSx.*(x(j)-o(1));     sigma2=sigma2+dSx.*(y(j)-o(2));
    sigma3=sigma3+dSy.*(x(j)-o(1));     sigma4=sigma4+dSy.*(y(j)-o(2));
    A=A+dSx.*(f(j)-f_ref);     B=B+dSy.*(f(j)-f_ref);
end
dfx=(sigma4.*A-sigma2.*B)./(sigma1.*sigma4-sigma2.*sigma3); % normalized GSM operator
end

