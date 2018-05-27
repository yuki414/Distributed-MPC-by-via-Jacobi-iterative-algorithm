clc
clear all
close all
solver = 'gurobi';  % gurobi or cplex
controller = 'dmpc'; % cmpc, dmpc or free
terminal = 'constraint' % constraint or cost
%% parameters
k1 = 0.4;   k2 = 0.3;
fs = 0.4;
Ts = 0.05;
RUNSTEP = 4/Ts;
m = 1;
M = 5; % a number of agent
N = 3; % predictive horizon
Qi = [100 0; 0 0]; % weighting matrix for state
Ri = 10; % for input
%% local dynamics
Aij = zeros(2,2);
Aim1 = [0,      0
        Ts*k2,  0];
Aii = [ 1,                  Ts
        Ts*(-k1-2*k2),   1-Ts*fs]; % nominal parameter is k1-2*k2
Aip1= Aim1;
Bij = zeros(2,1);
Bii = [0, Ts]';
%% global dynamics
n = size(Aii);  % nuber of subsystem
A = zeros(n*M);
B = zeros(n*M,M);
n_x = size(A,2); % size of state
n_u = size(B,2); % size of input

for i = 1:M %
    A(n*i-1:n*i,n*i-1:n*i) = Aii; % diagonal component
    if i < M
        A(n*(i+1)-1:n*(i+1),n*i-1:n*i) = Aim1; % 
        A(n*i-1:n*i,n*(i+1)-1:n*(i+1)) = Aip1; % 
    end
end
for i = 1:M
    B(n*i-1:n*i,i) = Bii;
end
C = eye(n*M);
D = zeros(n*M,M);
%% weighting matrices
Q = zeros(n*M);
R = zeros(M);
for i = 1:M
    Q(n*i-1:n*i,n*i-1:n*i) = Qi;
    R(i,i) = Ri;
end
%% lqr
[K,S,E] = dlqr(A,B,Q,R);
rand('state',1);
x0 = 2*(rand(n_x,1)-0.5);   % x_i=[position, velocity]';
% x0 = zeros(n_x,1);
% x0(1,1) = 1;
%% Centralized MPC
if(strcmpi(controller,'cmpc'))
mpc_pre_centralized
x_c = x0;
data_c = [];    data_cx=[];
draw_calc_time(0)
for step_c = 1:RUNSTEP
    tic
    mpc_calc_centralized
%     u_c=zeros(n_u,1);
    calc_time(step_c) = toc;
    x_c = A*x_c + B*u_c;
    data_c = [data_c; step_c-1, x_c'*Q*x_c + u_c'*R*u_c,norm(x_c),norm(u_c)];
    data_cx = [data_cx; x_c'];
    calc_time(step_c) = toc;
    draw_calc_time(1)
end
draw_calc_time(3)
end
%% Distributed MPC
if(strcmpi(controller,'dmpc'))
%% unique parameter for dmpc
x = x0;
omega = 1;
omega_i = omega/M; % set a mean
pmax = 2;   % a number of iteration
r = 1;  % neighborhood set
r = r*2 + 1;    % e.g. r=1 is {i-1,i,i+1}
AB = [];
for j = 1:N
    AB = blkdiag(AB,A,B);
end
%%
mpc_pre_distributed
% dmpc preparation for each subsystems
% especial case is on i=1,M.
% design each unique system to H, Aeq, Aineq and so on...
% but beq is constant value changed by step or iteration
% thereby beq is programmed later

z_p=zeros((n_x+n_u)*N,1); 
% X = x0(2*i-1:2*i,1);
% data_i = [];    data_d=[]; data_u=[];
% data_x=[];
% RUNSTEP = 1;
draw_calc_time(0)
for step_d = 1:RUNSTEP
    calc_time(step_i) = toc;
%         tic
    opt=zeros((n_x_i+n_u_i)*N,1);
    for p = 0:pmax % loop for iteration
        
        data_t=[];
        for i = 1:M
%             should review, because of z=[x', u']' but zi = [xi', ui']'
            z_mi = z_p;
            if i == 1
                z_mi(1:(n_x_i+n_u_i)*N*2,1) = 0;
            elseif i == M
                z_mi((n_x_i+n_u_i)*N*(M-1)+1:(n_x_i+n_u_i)*N*M) = 0;
            else
                z_mi((n_x_i+n_u_i)*N*(i-2)+1:(n_x_i+n_u_i)*N*(i+1)) = 0;
            end
            x_mi = AB*z_mi;
            mpc_calc_distributed
            draw_calc_time(2)
%             eval('z_', num2str(i), '_p = z_i_p;');
            % this is a case when optimal sequence of only i-th phiscal system
            % is extended to global sequence
            % [zeros((n_x_i+n_u_i)*N*(i-1),1)', z_i_p', zeros((n_x_i+n_u_i)*N*(M-i),1)']
            
%             z_i_p = zeros((n_x+n_u)*N); % reset by subsystem
            if i == 1
                z_i_p = z_mi + [z_i_ast', zeros((n_x_i+n_u_i)*N*(M-2),1)'];
            elseif i == M
                z_i_p = z_mi + [zeros((n_x_i+n_u_i)*N*(M-2),1)', z_i_ast'];
            else
                z_i_p = z_mi + [zeros((n_x_i+n_u_i)*N*(i-2),1)', z_i_ast', zeros((n_x_i+n_u_i)*N*(M-i-1),1)'];
            end
            z_pp1 = z_pp1 + omega_i*z_i_p;  % convex combination
        end
        z_p = z_pp1;
%         opt = z_p; % save z_p
%         opt1= opt;
        data_i = [data_i, z_p]; %  sequence of subsystem
        data_u = [data_u; ([zeros(n_u,n_u*k) eye(n_u) zeros(n_u,n_u*(N-1-k))]*[zeros(n_u*N,n_x*N), eye(n_u*N)]*opt)'];
    end
    k=0;
    u = [zeros(n_u,n_u*k) eye(n_u) zeros(n_u,n_u*(N-1-k))]*[zeros(n_u*N,n_x*N), eye(n_u*N)]*z_p;
    x = A*x + B*u;
    data_d = [data_d; step_i-1, x'*Q*x, u'*R*u, norm(x), norm(u)];
    data_x = [data_x; step_i, x'];
%     end
%     data_d = [data_d, data_i];
% draw_calc_time(2)
end
draw_calc_time(3)
end
%% calclation time
% mean_calc_time=mean(calc_time(2:RUNSTEP));
% [max_calc_time, I]=max(calc_time(2:RUNSTEP));
% X1 = sprintf('mean : %8.2f [ms] ',1000*mean_calc_time);
% X2 = sprintf('max  : %8.2f [ms] (%4.0f [step])',1000*max_calc_time, I);
% disp('calc_time--')
% disp(X1); disp(X2);
%% plot
figure(1)
hold on; grid on
plot(data_c(:,1)*Ts,data_c(:,2),'k')
plot(data_d(:,1)*Ts,data_d(:,2),'r:')

figure(2)
hold on; grid on
plot(data_c(:,1)*Ts,data_c(:,2)+data_c(:,3),'k')
plot(data_d(:,1)*Ts,data_d(:,2)+data_d(:,3),'r:')

figure(3)
hold on; grid on
for i=1:M
    plot(data_x(:,1)*Ts,data_x(:,2*i+1))
end

figure(4)
hold on; grid on
for i=1:M
    plot(data_c(:,1)*Ts,data_c(:,2*i+1))
end