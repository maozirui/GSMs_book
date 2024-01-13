function [] = plot_particles( N,x,y,o,angle)
X=x; Y=y;
XX=X;YY=Y;
k=0;
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

if k==0
    X=[x; x(1)]; Y=[y; y(1)];
else
    if k~=N
        j=k;
        for i=1:(N-j)
            X(i)=XX(j+i); Y(i)=YY(j+i);
        end
        for i=1:j
            X(N-j+i)=XX(i); Y(N-j+i)=YY(i);
        end
    end
end
        
Ox=[o(1); x(1)]; Oy=[o(2); y(1)];
for i=2:N
    Ox=[Ox;o(1);x(i)]; Oy=[Oy;o(2);y(i)];
end

H=zeros(1,4);
set(0,'defaultfigurecolor','w')    % set the background as white.
H(1:2)=plot(X,Y,'b-',Ox,Oy,'b-','linewidth',1.5);
hold on;
H(3)=plot(x,y,'bo','MarkerFaceColor','b','markersize',10);
hold on;
H(4)=plot(o(1),o(2),'ro','MarkerFaceColor','r','markersize',10);
legend(H([3 4]),'Neighboring particles','Target particles');
axis equal
xlabel('x'); ylabel('y'); 
set(gca,'FontName','TimesRoman','FontSize',14);