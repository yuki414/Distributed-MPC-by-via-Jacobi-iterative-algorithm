clc
clear all
close all
solver = 'gurobi';  % gurobi or cplex
%% parameters
k1 = -0.4;   k2 = 0.3;
fs = 0.4;
Ts = 0.05;
RUNSTEP = 4/Ts;
m = 1;
M = 5; % a number of agent
N = 10; % predictive horizon
Qi = [100 0; 0 0]; % weighting matrix for state
Ri = 10; % for input
%% local dynamics
Aij = zeros(2,2);
Aim1 = [0,      0
        Ts*k2,  0];
Aii = [ 1,              Ts
        Ts*(k1-2*k2),   1-Ts*fs];
Aip1= Aim1;
Bij = zeros(2,1);
Bii = [0, Ts]';
%% global dynamics
n = size(Aii);  %�T�u�V�X�e���̏�Ԑ�
A = zeros(n*M);
B = zeros(n*M,M);
n_x = size(A,2); %��Ԃ̎���
n_u = size(B,2); %���͂̎���

for i = 1:M %
    A(n*i-1:n*i,n*i-1:n*i) = Aii; % �Ίp����
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
if 1
%     RUNSTEP = 10;
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
    data_c = [data_c; step_c-1, x_c'*Q*x_c , u_c'*R*u_c,norm(x_c),norm(u_c)];
    data_cx = [data_cx; x_c'];
    calc_time(step_c) = toc;
    draw_calc_time(1)
end
draw_calc_time(3)
end
%% Distributed MPC
if 0
omega = 1;
omega_i = omega/M; % �ʌ����F����
pmax = 2;   % a number of iteration
z_p=[]; 
x = x0;
X = x0(2*i-1:2*i,1);
data_i = []; data_d=[]; data_u=[]; data_x=[];
% RUNSTEP = 1;
% draw_calc_time(0)
for step_i = 1:RUNSTEP
%     calc_time(step_i) = toc;
%         tic
    opt=zeros((n_x+n_u)*N,1); % �œK���̓��ꕨ
    for p = 0:pmax % loop for iteration
        z_p=0; % iteration���Ƃɏ�����
        data_t=[];
        for i = 1:M
            switch i
                case 1
                    mpc_pre_distributed_1 % ������
                case M 
                    mpc_pre_distributed_M
                otherwise
                    mpc_pre_distributed
            end
            mpc_calc_distributed_iteration
            draw_calc_time(2)
            opt1 = opt;
            z_p = z_p + omega_i*z_i_p; % �ʌ���
        end
%         norm(z_p)
        opt = z_p; % z_p�̕ۑ�
%         opt1= opt;
        data_i = [data_i, opt]; % �e�T�u�V�X�e���̉��n��
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