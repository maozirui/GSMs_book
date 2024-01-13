% This code is writen to approximate gradient of scalar function F 
% using GSM gradient operator.
% Author: Zirui Mao (maozirui0827@gmail.com)
% Last Updated: Sept., 2020
clc; clear all; clf;
%%%%%%%% Controlling parameters %%%%%%%%%%
o=[1.0 1.0]; % 2D coordinate of target node/particle
type_support=2; % type switch of supporting particles
                % = 1 uniformly distributed supporting particles
                % = 2 irregular distribution with small distance variation
                % = 3 irregular distribution with large distance variation
                % = 4 boundary particle with particle deficiency 
r=0.001; % define the spacing size
f_ref=o(1)^3; % define the testing scalar function F
dfx_ref=3*o(1)^2; % give the exact solution of F's gradient 
%%%%%%%% generate supporting particles %%%%%%%%%%%%%%%
if type_support==1 % uniformly distributed supporting particles
    N=6; % amount of supporting particles
    factor_length=[1 1 1 1 1 1]; % normalized length of edges
    factor_angle=[1 1 1 1 1 1]; % normzlied angles between each pair of adjacent supporting particles
elseif type_support==2 % irregular distribution with small distance variation
    N=6;
    factor_length=[1 1.05 0.83 1.0 0.85 1.1];
    factor_angle=[1 1.1 0.78 1.2 0.83 1];
elseif type_support==3 %irregular distribution with large distance variation
    N=6;
    factor_length=[1 5 5.5 1.0 0.85 1.1];
    factor_angle=[1 1.1 0.78 1.2 0.83 1];
elseif type_support==4 % boundary particle with particle deficiency 
    N=3;
    factor_length=[1 1.1 0.95 1.0 0.85 1.1];
    factor_angle=[1 1.1 0.78 1.2 0.83 1];
end
x=zeros(N,1); y=zeros(N,1); Length=r*factor_length; angle=zeros(N,1);
f=zeros(N,1); % scalar function F
for i=1:N
    angle(i)=360/sum(factor_angle)*(sum(factor_angle(1:i)));
    x(i)=o(1)+Length(i)*cosd(angle(i)); % specify the x coordinates of supporting particles
    y(i)=o(2)+Length(i)*sind(angle(i)); % specify the y coordinates of supporting particles
    f(i)=x(i).^3; % define the function values of supporting particles 
end
plot_particles( N,x,y,o,angle);%%%%%%%%%%% plot the particles configuration
%%%%%%%%% Gradient approximation with GSM gradient operator %%%%%%%%%%%%%%%
[dfx_1] = GSM_standard( N,x,y,o,angle,f,f_ref); % calculate gradient with standard GSM operator
[dfx_2] = GSM_normalized( N,x,y,o,angle,f,f_ref); % calcuate gradient with nomalized GSM operator
error1=abs(dfx_1-dfx_ref)/dfx_ref*100; error2=abs(dfx_2-dfx_ref)/dfx_ref*100;
%%%%%%%%% print out the solutions and compare to the exact solution %%%%%%%
fprintf('Standard GSM operator = %0.2f | Exact solution = %0.2f | Error percentage of Standard GSM operator = %0.2e \n', dfx_1,dfx_ref,error1);
fprintf('Normalized GSM operator = %0.2f | Exact solution = %0.2f | Error percentage of Normalized GSM operator = %0.2e \n', dfx_2,dfx_ref, error2);