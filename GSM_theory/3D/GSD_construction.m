function [dsx,dsy,dsz] = GSD_construction (xx, yy, zz)
% This function is writen to construct 3D n-GSD structure for 3D L-GSM.
% For more detail, please refer to ""Zirui Mao, et al., IInt J Numer Methods Eng. 
% 2020. vol. 121, no. 6, pp. 1268–1296.""
X=xx; Y=yy; Z=zz; 
N=length(X); % total number of neighbors
x0=X; y0=Y; z0=Z;
%%%%%%%%%%%%%% step 1: normalized the distances %%%%%%%%%%%%%%%%%
for i=1:N
    l(i)=sqrt(X(i)^2+Y(i)^2+Z(i)^2);
    x(i)=X(i)/l(i); y(i)=Y(i)/l(i);z(i)=Z(i)/l(i);
end
xt=0; yt=0; zt=0; % coordinate of target particle
NS=0; NE=0;
%%%%%%%% step 2:  find the first boundary triangle %%%%%%%%%%%%%%%%%%
i1=1; min_dis=1000;
for i2=2:N
    distance=sqrt((x(i1)-x(i2))^2+(y(i1)-y(i2))^2+(z(i1)-z(i2))^2);
    if distance<min_dis
        min_dis=distance;
        ii=i2;
    end
end
i2=ii;
V=0; surf_area=zeros(N,3);
if N>3
    ind=0;
    for i=2:N-1
        if i>=ind
            sig=1;
            A0=[0 0 0 1; x(i1) y(i1) z(i1) 1; x(i2) y(i2) z(i2) 1; x(i) y(i) z(i) 1];
            V0=det(A0)/6;
            for j=i+1:N
                A1=[x(j) y(j) z(j) 1; x(i1) y(i1) z(i1) 1; x(i2) y(i2) z(i2) 1; x(i) y(i) z(i) 1];
                V1=det(A1)/6;
                if V0*V1<0
                    sig=-1; ind=j;
                    break
                end
            end
            if sig>0&&abs(V0)>0
                break;
            end
        end
        if i==N-1
            i=N;
        end
    end
else
    i1=1;i2=2;i=3;
    A0=[0 0 0 1; x(i1) y(i1) z(i1) 1; x(i2) y(i2) z(i2) 1; x(i) y(i) z(i) 1];
    V0=det(A0)/6;
end
if abs(V0)>0 % is triangle
    NS=1; % number of surface
    surface(NS, 1:3)=[i1 i2 i];
    NE=3;
    A=[0 0 0 1; x0(i1) y0(i1) z0(i1) 1; x0(i2) y0(i2) z0(i2) 1; x0(i) y0(i) z0(i) 1];
    V=V+abs(det(A))/24;
    x11=x0(i1)/2; y11=y0(i1)/2; z11=z0(i1)/2;
    x12=(x0(i1)+x0(i2))/3; y12=(y0(i1)+y0(i2))/3; z12=(z0(i1)+z0(i2))/3;
    x13=(x0(i1)+x0(i)+x0(i2))/4; y13=(y0(i1)+y0(i)+y0(i2))/4; z13=(z0(i1)+z0(i)+z0(i2))/4;
    x21=x0(i2)/2; y21=y0(i2)/2; z21=z0(i2)/2;
    x22=(x0(i2)+x0(i))/3; y22=(y0(i2)+y0(i))/3; z22=(z0(i2)+z0(i))/3;
    x23=(x0(i2)+x0(i1)+x0(i))/4; y23=(y0(i2)+y0(i1)+y0(i))/4; z23=(z0(i2)+z0(i1)+z0(i))/4;
    x31=x0(i)/2; y31=y0(i)/2; z31=z0(i)/2;
    x32=(x0(i)+x0(i1))/3; y32=(y0(i)+y0(i1))/3; z32=(z0(i)+z0(i1))/3;
    x33=(x0(i)+x0(i2)+x0(i1))/4; y33=(y0(i)+y0(i2)+y0(i1))/4; z33=(z0(i)+z0(i2)+z0(i1))/4;   
    a13=sqrt((x11-x12)^2+(y11-y12)^2+(z11-z12)^2);
    a11=sqrt((x12-x13)^2+(y13-y12)^2+(z13-z12)^2);
    a12=sqrt((x11-x13)^2+(y11-y13)^2+(z11-z13)^2);
    a23=sqrt((x21-x22)^2+(y21-y22)^2+(z21-z22)^2);
    a21=sqrt((x22-x23)^2+(y23-y22)^2+(z23-z22)^2);
    a22=sqrt((x21-x23)^2+(y21-y23)^2+(z21-z23)^2);
    a33=sqrt((x31-x32)^2+(y31-y32)^2+(z31-z32)^2);
    a31=sqrt((x32-x33)^2+(y33-y32)^2+(z33-z32)^2);
    a32=sqrt((x31-x33)^2+(y31-y33)^2+(z31-z33)^2);
    p3=a31+a32+a33; p3=p3/2;  p1=(a11+a12+a13)/2;  p2=(a21+a22+a23)/2;
    area1=sqrt(p1*(p1-a11)*(p1-a12)*(p1-a13))*2;
    area2=sqrt(p2*(p2-a21)*(p2-a22)*(p2-a23))*2;
    area3=sqrt(p3*(p3-a31)*(p3-a32)*(p3-a33))*2;
    l1(1)=(y12-y11)*(z13-z11)-(y13-y11)*(z12-z11);
    l1(2)=(z12-z11)*(x13-x11)-(z13-z11)*(x12-x11);
    l1(3)=(x12-x11)*(y13-y11)-(x13-x11)*(y12-y11);
    l2(1)=(y22-y21)*(z23-z21)-(y23-y21)*(z22-z21);
    l2(2)=(z22-z21)*(x23-x21)-(z23-z21)*(x22-x21);
    l2(3)=(x22-x21)*(y23-y21)-(x23-x21)*(y22-y21);
    l3(1)=(y32-y31)*(z33-z31)-(y33-y31)*(z32-z31);
    l3(2)=(z32-z31)*(x33-x31)-(z33-z31)*(x32-x31);
    l3(3)=(x32-x31)*(y33-y31)-(x33-x31)*(y32-y31);   
    ll1=sqrt(l1(1)^2+l1(2)^2+l1(3)^2);
    ll2=sqrt(l2(1)^2+l2(2)^2+l2(3)^2);
    ll3=sqrt(l3(1)^2+l3(2)^2+l3(3)^2);
    l1=l1/ll1; l2=l2/ll2; l3=l3/ll3;    
    if V0<0
        edge(1,1:2)=[i1 i2]; edge(2,1:2)=[i2 i]; edge(3,1:2)=[i i1];
        surf_area(i1,1:3)=surf_area(i1,1:3)+area1*l1;
        surf_area(i2,1:3)=surf_area(i2,1:3)+area2*l2;
        surf_area(i,1:3)=surf_area(i,1:3)+area3*l3;
    else
        edge(1,1:2)=[i2 i1]; edge(2,1:2)=[i i2]; edge(3,1:2)=[i1 i];
        surf_area(i1,1:3)=surf_area(i1,1:3)-area1*l1;
        surf_area(i2,1:3)=surf_area(i2,1:3)-area2*l2;
        surf_area(i,1:3)=surf_area(i,1:3)-area3*l3;
    end
    edge_particle(1,1)=i; edge_particle(2,1)=i1; edge_particle(3,1)=i2;
end
Edge=edge; edge_used=0;
%%%%%%%%%%%%%%% step 3: find the subsequent boundary triangle %%%%%%%%%%%%%%%%%
while ~isempty(edge)
    i1=edge(1,1); i2=edge(1,2);
    min_i=0;
    for i=1:N
        if i~=i1&&i~=i2&&i~=min_i&&i~=edge_particle(1)
            belong_used=0;
            for j=1:length(edge_used(:,1))
                if i==edge_used(j,1)&&i1==edge_used(j,2)
                    belong_used=1;
                    break
                elseif i1==edge_used(j,1)&&i==edge_used(j,2)
                    belong_used=1;
                    break
                elseif i==edge_used(j,1)&&i2==edge_used(j,2)
                    belong_used=1;
                    break
                elseif i2==edge_used(j,1)&&i==edge_used(j,2)
                    belong_used=1;
                    break
                end
            end
            if belong_used==0
                if min_i==0
                    min_i=i;
                else
                    A0=[x(i) y(i) z(i) 1; x(i1) y(i1) z(i1) 1; x(i2) y(i2) z(i2) 1; x(min_i) y(min_i) z(min_i) 1] ;
                    if det(A0)<0
                        min_i=i;
                    end
                end
            end
        end
    end
    if min_i==0
        edge=0;
        break
    end
    A0=[0 0 0 1; x(i1) y(i1) z(i1) 1; x(i2) y(i2) z(i2) 1; x(min_i) y(min_i) z(min_i) 1];
    V0=det(A0)/6;
    if V0>0  % inner
        NS=NS+1;
        A=[0 0 0 1; x0(i1) y0(i1) z0(i1) 1; x0(i2) y0(i2) z0(i2) 1; x0(min_i) y0(min_i) z0(min_i) 1];
        V=V+abs(det(A))/24;
        x11=x0(i1)/2; y11=y0(i1)/2; z11=z0(i1)/2;
        x12=(x0(i1)+x0(i2))/3; y12=(y0(i1)+y0(i2))/3; z12=(z0(i1)+z0(i2))/3;
        x13=(x0(i1)+x0(min_i)+x0(i2))/4; y13=(y0(i1)+y0(min_i)+y0(i2))/4; z13=(z0(i1)+z0(min_i)+z0(i2))/4;
        x21=x0(i2)/2; y21=y0(i2)/2; z21=z0(i2)/2;
        x22=(x0(i2)+x0(min_i))/3; y22=(y0(i2)+y0(min_i))/3; z22=(z0(i2)+z0(min_i))/3;
        x23=(x0(i2)+x0(i1)+x0(min_i))/4; y23=(y0(i2)+y0(i1)+y0(min_i))/4; z23=(z0(i2)+z0(i1)+z0(min_i))/4;
        x31=x0(min_i)/2; y31=y0(min_i)/2; z31=z0(min_i)/2;
        x32=(x0(min_i)+x0(i1))/3; y32=(y0(min_i)+y0(i1))/3; z32=(z0(min_i)+z0(i1))/3;
        x33=(x0(min_i)+x0(i2)+x0(i1))/4; y33=(y0(min_i)+y0(i2)+y0(i1))/4; z33=(z0(min_i)+z0(i2)+z0(i1))/4;       
        a13=sqrt((x11-x12)^2+(y11-y12)^2+(z11-z12)^2);
        a11=sqrt((x12-x13)^2+(y13-y12)^2+(z13-z12)^2);
        a12=sqrt((x11-x13)^2+(y11-y13)^2+(z11-z13)^2);
        a23=sqrt((x21-x22)^2+(y21-y22)^2+(z21-z22)^2);
        a21=sqrt((x22-x23)^2+(y23-y22)^2+(z23-z22)^2);
        a22=sqrt((x21-x23)^2+(y21-y23)^2+(z21-z23)^2);
        a33=sqrt((x31-x32)^2+(y31-y32)^2+(z31-z32)^2);
        a31=sqrt((x32-x33)^2+(y33-y32)^2+(z33-z32)^2);
        a32=sqrt((x31-x33)^2+(y31-y33)^2+(z31-z33)^2);
        p3=a31+a32+a33; p3=p3/2; p1=(a11+a12+a13)/2; p2=(a21+a22+a23)/2;
        area1=sqrt(p1*(p1-a11)*(p1-a12)*(p1-a13))*2;
        area2=sqrt(p2*(p2-a21)*(p2-a22)*(p2-a23))*2;
        area3=sqrt(p3*(p3-a31)*(p3-a32)*(p3-a33))*2;        
        l1(1)=(y12-y11)*(z13-z11)-(y13-y11)*(z12-z11);
        l1(2)=(z12-z11)*(x13-x11)-(z13-z11)*(x12-x11);
        l1(3)=(x12-x11)*(y13-y11)-(x13-x11)*(y12-y11);
        l2(1)=(y22-y21)*(z23-z21)-(y23-y21)*(z22-z21);
        l2(2)=(z22-z21)*(x23-x21)-(z23-z21)*(x22-x21);
        l2(3)=(x22-x21)*(y23-y21)-(x23-x21)*(y22-y21);
        l3(1)=(y32-y31)*(z33-z31)-(y33-y31)*(z32-z31);
        l3(2)=(z32-z31)*(x33-x31)-(z33-z31)*(x32-x31);
        l3(3)=(x32-x31)*(y33-y31)-(x33-x31)*(y32-y31);        
        ll1=sqrt(l1(1)^2+l1(2)^2+l1(3)^2);
        ll2=sqrt(l2(1)^2+l2(2)^2+l2(3)^2);
        ll3=sqrt(l3(1)^2+l3(2)^2+l3(3)^2);
        l1=l1/ll1; l2=l2/ll2; l3=l3/ll3;        
        surf_area(i1,1:3)=surf_area(i1,1:3)-area1*l1;
        surf_area(i2,1:3)=surf_area(i2,1:3)-area2*l2;
        surf_area(min_i,1:3)=surf_area(min_i,1:3)-area3*l3;
        surface0=sort([i1,i2,min_i]);
        surface(NS, 1:3)=surface0;
        if edge_used~=0
            edge_used=[edge_used; edge(1,:)];
        else
            edge_used=edge(1,:);
        end
        edge(1,:)=[]; edge_particle(1)=[];
        belong_i1=0; belong_i2=0;
        for j=1:length(Edge(:,1)) % examine whether [i1 min] [min i2] in edge?
            if (Edge(j,1)==i1&&Edge(j,2)==min_i)||(Edge(j,1)==min_i&&Edge(j,2)==i1)
                belong_i1=1;
            end
            if (Edge(j,1)==min_i&&Edge(j,2)==i2)||(Edge(j,2)==min_i&&Edge(j,1)==i2)
                belong_i2=1;
            end
        end
        if belong_i1==0
            edge=[edge; i1 min_i];  edge_particle=[edge_particle; i2];
        else
            for i=1:length(edge(:,1))
                if (edge(i,1)==i1&&edge(i,2)==min_i)||(edge(i,2)==i1&&edge(i,1)==min_i)
                    if edge_used~=0
                        edge_used=[edge_used; edge(i,:)];
                    else
                        edge_used=edge(i,:);
                    end
                    edge(i,:)=[]; edge_particle(i)=[];
                    break
                end
            end
        end
        if belong_i2==0
            edge=[edge; min_i i2];   edge_particle=[edge_particle; i1];
        else
            for i=1:length(edge(:,1))
                if (edge(i,1)==i2&&edge(i,2)==min_i)||(edge(i,2)==i2&&edge(i,1)==min_i)
                    if edge_used~=0
                        edge_used=[edge_used; edge(i,:)];
                    else
                        edge_used=edge(i,:);
                    end
                    edge(i,:)=[]; edge_particle(i)=[];
                    break
                end
            end
        end
        Edge=[Edge; i1 min_i; min_i i2];
    else
        if edge_used~=0
            edge_used=[edge_used; edge(1,:)];
        else
            edge_used=edge(1,:);
        end
        edge(1,:)=[];  edge_particle(1)=[];
        for i=1:length(edge(:,1))
            if (edge(i,1)==i1&&edge(i,2)==min_i)||(edge(i,2)==i1&&edge(i,1)==min_i)
                if edge_used~=0
                    edge_used=[edge_used; edge(i,:)];
                else
                    edge_used=edge(i,:);
                end
                edge(i,:)=[];   edge_particle(i)=[];
                break
            end
        end
        for i=1:length(edge(:,1))
            if (edge(i,1)==i2&&edge(i,2)==min_i)||(edge(i,2)==i2&&edge(i,1)==min_i)
                if edge_used~=0
                    edge_used=[edge_used; edge(i,:)];
                else
                    edge_used=edge(i,:);
                end
                edge(i,:)=[];   edge_particle(i)=[];
                break
            end
        end    
    end
end
x=x0; y=y0; z=z0;
E11=0; E12=0; E13=0; E21=0; E22=0; E23=0; E31=0; E32=0; E33=0;
for j=1:N
    E11=E11+x(j)*surf_area(j,1); E12=E12+y(j)*surf_area(j,1); E13=E13+z(j)*surf_area(j,1);
    E21=E21+x(j)*surf_area(j,2); E22=E22+y(j)*surf_area(j,2); E23=E23+z(j)*surf_area(j,2);
    E31=E31+x(j)*surf_area(j,3); E32=E32+y(j)*surf_area(j,3); E33=E33+z(j)*surf_area(j,3);
end
CCz=(E32-E12*E31/E11)/(E22-E12*E21/E11); CCy=(E33-E13*E31/E11)/(E23-E13*E21/E11); CCx=(E13-E23*E12/E22)/(E33-E23*E32/E22);
Cz=E33-E13*E31/E11-(E23-E13*E21/E11)*CCz; Cy=E32-E12*E31/E11-(E22-E12*E21/E11)*CCy; Cx=E11-E21*E12/E22-(E31-E21*E32/E22)*CCx;
dsz=((surf_area(:,3)-surf_area(:,1)*E31/E11)-(surf_area(:,2)-surf_area(:,1)*E21/E11)*CCz)/Cz;
dsy=((surf_area(:,3)-surf_area(:,1)*E31/E11)-(surf_area(:,2)-surf_area(:,1)*E21/E11)*CCy)/Cy;
dsx=((surf_area(:,1)-surf_area(:,2)*E12/E22)-(surf_area(:,3)-surf_area(:,2)*E32/E22)*CCx)/Cx;
%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%
ox=[xt;x(1)]; oy=[yt;y(1)]; oz=[zt; z(1)];
plot3(x,y,z,'bo','MarkerFaceColor','b','markersize',10); hold on;
plot3(xt,yt,zt,'ro','MarkerFaceColor','r','markersize',10); hold on;
NN=length(surface(:,1));
set(gca,'FontName','times new Roman','FontSize',18);
axis equal; grid on;
xlabel('x'); ylabel('y'); zlabel('z'); legend('Supporting particles','target particle')
end

