% This code is writen to solve the 3D gradient approximation of scalar F 
% using 3D normalized GSM gradient operator. For detailed information, please refer to
% "Zirui Mao, et al., IInt J Numer Methods Eng. 2020. vol. 121, no. 6, pp. 1268–1296."
% Author: Zirui Mao (maozirui0827@gmail.com)
% Last Updated: Sept., 2020
clc; clear all; clf;
%%%%%%%% controlling parameters %%%%%%%%%
o1=1; o2=1; o3=1; % define 3D coordinates of target particle
type_support=1; % type switch of supporting particles
                % = 1 regular box
                % = 2 mixed convex-concave domain
                % = 3 corner boundary particle
r=0.01; % define the spacing size
f=o1^2/2+o2^3/3+o3^4/4+10*(o1+o2+o3); % testing scalar function F
dfx_theory=o1+10; dfy_theory=o2^2+10; dfz_theory=o3^3+10; % gradient of F in theory
%%%%%%% generate positions of neighbors relative to  target particle %%%%%%%%
if type_support == 1 % regular box
    x=[0; 1; 1.; 0; -1; -1.; -1; 0; 1; 1; 1; 0; -1.; -1; -1; 0; 1; 0; 1; 1; 0; -1; -1.; -1; 0; 1];
    y=[0; 0; 1.; 1; 1; 0; -1; -1; -1; 0; 1; 1.; 1; 0; -1; -1; -1; 0; 0; 1; 1; 1.; 0; -1; -1.; -1];
    z=[-1; -1.;  -1;  -1; -1; -1; -1.; -1; -1; 0; 0; 0; 0; 0; 0; 0; 0; 1.; 1; 1; 1.; 1; 1; 1; 1; 1];
elseif type_support == 2 % mixed convex-concave domain
    x=[0; 1.5; 0.9; 0; -1; -1.; -1; 0; 1; 0.8; 1; 0; -1.; -0.9; -1; -0.2; 1; 0; 1; 1; 0; -1; -1.5; -1; 0; 1];
    y=[0; 0; 1.; 1; 1.5; 0; -1; -0.8; -1; 0; 1; 1.2; 1; 0; -1; -1; -1; 0; 0; 1.2; 1; 0.7; 0; -1.5; -1.; -1];
    z=[-1.5; -1.1; -1; -1; -0.75; -1; -1.; -1.2; -1; 0; 0; 0.3; 0; 0.2; 0; 0; -0.3; 0.8; 1; 1; 0.75; 1; 1.2; 1; 1; 1.5];
elseif type_support == 3 % corner
    x=[0;   1;   1.25;  0;  -0.75];
    y=[0;   0;   1.25;  1;  0.75];
    z=[-1; -1.;  -1.25;  -1; -0.75];
end
x=x*r; y=y*r; z=z*r; % generate the relative position of supporting particles
%%%%%%%%%%%%%% 3D n-GSD construction %%%%%%%%%%%%%%
[dsx,dsy,dsz] = GSD_construction (x, y, z); % output dS_x, dS_y, and dS_z in GSM operator
% %%%% Approximation of gradient with the 3D normalized GSM operation %%%%%%%%%%%%
dfx=0; dfy=0; dfz=0; % gradient of f in x, y, and z directions, respectively
for k=1:length(x)
    f_neighbor=(o1+x(k))^2/2+(o2+y(k))^3/3+(o3+z(k))^4/4+10*(o1+x(k)+o2+y(k)+o3+z(k));
    dfx=dfx+(f_neighbor-f)*dsx(k); % x componenet of gradient of F
    dfy=dfy+(f_neighbor-f)*dsy(k); % y componenet of gradient of F
    dfz=dfz+(f_neighbor-f)*dsz(k); % z componenet of gradient of F
end
error=abs(dfz-dfz_theory)/dfz_theory*100;
fprintf('dF_x: GSM solution = %0.2f | Exact solution = %0.2f \n', dfx,dfx_theory);
fprintf('dF_y: GSM solution = %0.2f | Exact solution = %0.2f \n', dfy,dfy_theory);
fprintf('dF_z: GSM solution = %0.2f | Exact solution = %0.2f | Error percentage of GSM solution = %0.2e \n', dfz,dfz_theory, error);