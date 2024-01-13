function [dFx,dFy] = GSM_gradient (T, x, y, N, F)
% This function is written to approximate gradient of scalar function F
% by GSM gradient operator based on n-GSD.
NN=length(T); % total number of triangles representing the whole domain
dSx_F=zeros(N,1); dSy_F=zeros(N,1); % DeltaS*F in the GSM gradient operator
V=zeros(N,1); % volume or area of GSD in the GSM gradient operator
for t=1:NN
    i1=T(t,1); i2=T(t,2); i3=T(t,3); % indicies are ordered in anticlockwise sequence
    area=abs((x(i2)-x(i1))*(y(i3)-y(i1)) - (x(i3)-x(i1))*(y(i2)-y(i1)))/6;
    V(i1)=V(i1)+area; V(i2)=V(i2)+area; V(i3)=V(i3)+area; % area of GSD
    dSx_F(i1)=dSx_F(i1)+(y(i3)/3.0-y(i1)/6.0-y(i2)/6.0)*(F(i2)+F(i1))/2.0;
    dSy_F(i1)=dSy_F(i1)+(x(i3)/3.0-x(i1)/6.0-x(i2)/6.0)*(F(i2)+F(i1))/2.0;
    dSx_F(i1)=dSx_F(i1)+(y(i1)/6.0+y(i3)/6.0-y(i2)/3.0)*(F(i3)+F(i1))/2.0;
    dSy_F(i1)=dSy_F(i1)+(x(i1)/6.0+x(i3)/6.0-x(i2)/3.0)*(F(i3)+F(i1))/2.0;
    dSx_F(i2)=dSx_F(i2)+(y(i1)/3.0-y(i2)/6.0-y(i3)/6.0)*(F(i3)+F(i2))/2.0;
    dSy_F(i2)=dSy_F(i2)+(x(i1)/3.0-x(i2)/6.0-x(i3)/6.0)*(F(i3)+F(i2))/2.0;
    dSx_F(i2)=dSx_F(i2)+(y(i1)/6.0+y(i2)/6.0-y(i3)/3.0)*(F(i1)+F(i2))/2.0;
    dSy_F(i2)=dSy_F(i2)+(x(i1)/6.0+x(i2)/6.0-x(i3)/3.0)*(F(i1)+F(i2))/2.0;
    dSx_F(i3)=dSx_F(i3)+(y(i2)/3.0-y(i1)/6.0-y(i3)/6.0)*(F(i1)+F(i3))/2.0;
    dSy_F(i3)=dSy_F(i3)+(x(i2)/3.0-x(i1)/6.0-x(i3)/6.0)*(F(i1)+F(i3))/2.0;
    dSx_F(i3)=dSx_F(i3)+(y(i3)/6.0+y(i2)/6.0-y(i1)/3.0)*(F(i2)+F(i3))/2.0;
    dSy_F(i3)=dSy_F(i3)+(x(i3)/6.0+x(i2)/6.0-x(i1)/3.0)*(F(i2)+F(i3))/2.0;
end
dFx=dSx_F./V; dFy=-dSy_F./V; % GSM gradient operator for gradient components of F
end