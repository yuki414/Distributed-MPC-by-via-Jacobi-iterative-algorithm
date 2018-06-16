clc
clear all
close all
solver = 'gurobi';  % gurobi or cplex
controller = 'dmpc'; % cmpc, dmpc or free
terminal = 'constraint'; % constraint or cost
%% parameters
k1 = 0.4;   k2 = 0.3;
fs = 0.4;
Ts = 0.05;
RUNSTEP = 100;
m = 1;
M = 5; % a number of agent
N = 10; % predictive horizon
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
c_fig = 0;
d_fig = 0;
% x0 = zeros(n_x,1);
% x0(1,1) = 1;
%% Centralized MPC
if(strcmpi(controller,'cmpc'))
c_fig = 0;
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
save data_of_candidate Usol % CMPC‚ð‰ñ‚µ‚Äx~,u~‚ÌŒó•â‚ð“¾‚é
end
%% Distributed MPC
if(strcmpi(controller,'dmpc'))
d_fig = 0;
%% unique parameter for dmpc
x = x0;
omega = 1;
omega_i = omega/M; % set a mean
pmax = 5;   % a number of iteration
r = 1;  % neighborhood set
r = r*2 + 1;    % e.g. r=1 is {i-1,i,i+1}
AB = [];
% ˜_•¶’†‚ÌƒJƒŠƒOƒ‰ƒtƒVƒXƒeƒ€s—ñ
Ak = kron(eye(N),A);    % diag(A,A,...)
Bk = kron(eye(N),B);    % diag(B,B,...)
%%
mpc_pre_distributed
% dmpc preparation for each subsystems
% especial case is on i=1,M
% design each unique system to H, Aeq, Aineq and so on...
% but beq is constant value changed by step or iteration
% thereby beq is programmed later
load data_of_candidate Usol % CMPC‚Ì‰ðŒn—ñ‚ð“¾‚é
z_p = Usol;
data_i = []; data_d=[]; data_u=[]; data_x=[];
% RUNSTEP = 1;
% draw_calc_time(0)
for step_d = 1:RUNSTEP
%     calc_time(step_i) = toc;
%         tic
    z_mi = z_p; % ˆê‚Â‘O‚Ì‘S‘Ì‚Ì‰ðŒn—ñ
    for p = 0:pmax % starts iteration
        X = z_p(1:n_x*N,1);     % ‰ðŒn—ñ‚©‚çX‚ð”²‚«o‚·
        U = z_p(n_x*N+1:end,1); % ‰ðŒn—ñ‚©‚çU‚ð”²‚«o‚·
        X_A = Ak*X; % §–ñ—p
        X_i = zeros(n_x_i*N,M); % {i-2,..,i+2}‚Ìx‚Ì‰ðŒn—ñ‚ð•À‚×‚é‚½‚ß‚Ì‚¢‚ê‚à‚Ì
        X_A_i = zeros(n_x_i*N,M); 
        U_i = zeros(n_u_i*N,M); % {i-2,..,i+2}‚Ìu‚Ì‰ðŒn—ñ‚ð•À‚×‚é‚½‚ß‚Ì‚¢‚ê‚à‚Ì
        % i—ñ–Ú‚ªi”Ô–Ú‚ÌƒTƒuƒVƒXƒeƒ€‚Ì‰ðŒn—ñ‚É‚È‚é‚æ‚¤‚É¬Œ`
        for i = 1:M
            % ‰ðŒn—ñ‚ð“ü‚ê‚Ä‚¢‚é
            for k = 1:N
                X_i(n_x_i*k-1:n_x_i*k,i) = X(n_x_i*(k-1)*M+(i-1)*n_x_i+1:n_x_i*(k-1)*M+(i-1)*n_x_i+2,1);
                X_A_i(n_x_i*k-1:n_x_i*k,i) = X_A(n_x_i*(k-1)*M+(i-1)*n_x_i+1:n_x_i*(k-1)*M+(i-1)*n_x_i+2,1);
                U_i(n_u_i*k,i) = U(n_u_i*(k-1)*M+(i-1)*n_u_i+1,1);
            end
        end
        data_t=[];
        z_pp1 = 0;
        for i = 1:M
            mpc_calc_distributed % generate an optimal input
%           z_i_p‚Íi-1,i,i+1‚Ì‰ðŒn—ñ
%           ‚Ü‚¸ó‘Ô‚Æ“ü—Í‚É•ª‚¯‚é
%           ‚»‚Ì‚ ‚ÆXCU‚Ì‰ðŒn—ñ‚Ì‘Î‰ž‚·‚éêŠ‚É‘ã“ü
            if (i == 1)||(i == M) % case of i=1orM
                b = -round(i/M);   % i=1‚È‚çb=0 i=M‚È‚çb=-1
                % ¡‰ñ“¾‚½‰ðŒn—ñ‚Ì•ª‰ð
                % ‚ ‚Æ‚Å‘S‘Ì‚Ì‰ðŒn—ñ‚É“ü‚ê‚â‚·‚¢‚æ‚¤‚É‚·‚é‚½‚ß‚Éì¬
                Xs_i_p = z_i_p(1:n_x_i*N*2,1);
                Us_i_p = z_i_p(n_x_i*N*2+1:end);
                % ‘S‘Ì‚Ì‰ðŒn—ñ‚É‚¢‚ê‚é
                for j = 1:2
                    for k = 1:N
                        X_i(n_x_i*k-1:n_x_i*k,i+j+b-1) = Xs_i_p(n_x_i*(k-1)*2+(j-1)*n_x_i+1:n_x_i*(k-1)*2+(j-1)*n_x_i+2,1);
                        U_i(n_u_i*k,i+j+b-1) = Us_i_p(n_u_i*(k-1)*2+(j-1)*n_u_i+1,1);
                    end
                end
            else % case of i=2~M-1
                Xs_i_p = z_i_p(1:n_x_i*N*3,1);
                Us_i_p = z_i_p(n_x_i*N*3+1:end);
                for j = 1:3
                    for k = 1:N
                        X_i(n_x_i*k-1:n_x_i*k,i+j-2) = Xs_i_p(n_x_i*(k-1)*3+(j-1)*n_x_i+1:n_x_i*(k-1)*3+(j-1)*n_x_i+2,1);
                        U_i(n_u_i*k,i+j-2) = Us_i_p(n_u_i*(k-1)*3+(j-1)*n_u_i+1,1);
                    end
                end
            end
            % ‘S‘Ì‚Ì‰ðŒn—ñ‚É–ß‚· “‡
            for l = 1:M
                for k = 1:N
                    X_i_pp1(n_x_i*(k-1)*M+(l-1)*n_x_i+1:n_x_i*(k-1)*M+(l-1)*n_x_i+2,1) = X_i(n_x_i*k-1:n_x_i*k,l);
                    U_i_pp1(n_u_i*k+(l-1)*N*n_u_i,1) = U_i(n_u_i*k,l);
                end
            end
            z_i_pp1 = [X_i_pp1', U_i_pp1']';
%             draw_calc_time(2)
            z_pp1 = z_pp1 + omega_i*z_i_pp1;  % convex combination
        end
        z_p = z_pp1;
        data_i = [data_i, z_p]; % sequence of subsystem
    end
    k=0;
    u = [zeros(n_u,n_u*k) eye(n_u) zeros(n_u,n_u*(N-1-k))]*[zeros(n_u*N,n_x*N), eye(n_u*N)]*z_p; % apply input to plant
    x = A*x + B*u;
    data_d = [data_d; step_d-1, x'*Q*x + u'*R*u, norm(x), norm(u)];
    data_x = [data_x; step_d, x'];
end
% draw_calc_time(3)
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
if c_fig
plot(data_c(:,1)*Ts,data_c(:,2),'k')
end
if d_fig
plot(data_d(:,1)*Ts,data_d(:,2),'r:')
end

figure(3)
hold on; grid on
for i=1:M
    plot(data_x(:,1)*Ts,data_x(:,2*i+1))
end

% figure(4)
% hold on; grid on
% for i=1:M
%     plot(data_c(:,1)*Ts,data_c(:,2*i+1))
% end