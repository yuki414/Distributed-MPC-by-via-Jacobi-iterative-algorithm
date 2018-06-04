%% parameter
n_x_i = 2;  % size of subsystem's state
n_u_i = 1;  % input
xmax_ii = [2, inf]';
xmax_i = repmat(xmax_ii,r,1);    % constraints for state
xmin_i = -xmax_i;
xmax_1 = repmat(xmax_ii,r-1,1);  xmin_1 = -xmax_1;
xmax_M = xmax_1;    xmin_M = xmin_1;

Ak_1 = [Aii,    Aim1;
        Aim1,   Aii];
Bk_1 = [Bii, Bij;
        Bij, Bii];
Ak_M = Ak_1;
Bk_M = Bk_1;
Qk_1 = blkdiag(Qi,Qi);  Qk_M = Qk_1;
Rk_1 = blkdiag(Ri,Ri);  Rk_M = Rk_1;
%% local system matrix
Ak_i = zeros(n_x_i*r);
Ak_mi = zeros(n_x_i*(M-r));
Bk_i = zeros(n_x_i*r,n_u_i*r);
Qk_i = zeros(n_x_i*r);
Rk_i = zeros(n_u_i*r);
for i = 1:r %
    Ak_i(n*i-1:n*i,n*i-1:n*i) = Aii;
    if i < r
        Ak_i(n*(i+1)-1:n*(i+1),n*i-1:n*i) = Aim1; % 
        Ak_i(n*i-1:n*i,n*(i+1)-1:n*(i+1)) = Aip1; % 
    end
    Bk_i(n*i-1:n*i,i) = Bii;
    Qk_i(n*i-1:n*i,n*i-1:n*i) = Qi;
    Rk_i(i,i) = Ri;
end
%% setting
H_i = zeros(r*(n_x_i+n_u_i)*N);   
% objective function
% but i designed for all variable, H=zeros(n_x+n_u)*N
% hence consider only the variable for subsystem
f_i = zeros(r*(n_x_i+n_u_i)*N,1);

H_1 = zeros((r-1)*(n_x_i+n_u_i)*N);     H_M = H_1;
f_1 = zeros((r-1)*(n_x_i+n_u_i)*N,1);   f_M = f_1;

AineqC_i = [];  bineqC_i = [];
AineqC_1 = [];  bineqC_1 = [];
AineqC_M = [];  bineqC_M = [];
AeqC_i = [];    beqC_i = [];
AeqC_1 = [];    beqC_1 = [];
AeqC_M = [];    beqC_M = [];

% preparation for z^i
Tr_X_i = [eye(r*n_x_i*N), zeros(r*n_x_i*N,r*n_u_i*N)];
Tr_U_i = [zeros(r*n_u_i*N,r*n_x_i*N), eye(r*n_u_i*N)];
% for i=1
Tr_X_1 = [eye((r-1)*n_x_i*N), zeros((r-1)*n_x_i*N,(r-1)*n_u_i*N)];
Tr_U_1 = [zeros((r-1)*n_u_i*N,(r-1)*n_x_i*N), eye((r-1)*n_u_i*N)];
% for i=M
Tr_X_M = [eye((r-1)*n_x_i*N), zeros((r-1)*n_x_i*N,(r-1)*n_u_i*N)];
Tr_U_M = [zeros((r-1)*n_u_i*N,(r-1)*n_x_i*N), eye((r-1)*n_u_i*N)];
%% for loop
for k = 0:N-1
%% preparation for i
    Tr_x_i_kp1 = [zeros(r*n_x_i,r*n_x_i*k) eye(r*n_x_i) zeros(r*n_x_i,r*n_x_i*(N-1-k))]*Tr_X_i; % distributed x
    Tr_u_i_k = [zeros(r*n_u_i,r*n_u_i*k) eye(r*n_u_i) zeros(r*n_u_i,r*n_u_i*(N-1-k))]*Tr_U_i; % disributed u
    Tr_x_1_kp1 = [zeros((r-1)*n_x_i,(r-1)*n_x_i*k) eye((r-1)*n_x_i) zeros((r-1)*n_x_i,(r-1)*n_x_i*(N-1-k))]*Tr_X_1;
    Tr_u_1_k = [zeros((r-1)*n_u_i,(r-1)*n_u_i*k) eye((r-1)*n_u_i) zeros((r-1)*n_u_i,(r-1)*n_u_i*(N-1-k))]*Tr_U_1;
    Tr_x_M_kp1 = [zeros((r-1)*n_x_i,(r-1)*n_x_i*k) eye((r-1)*n_x_i) zeros((r-1)*n_x_i,(r-1)*n_x_i*(N-1-k))]*Tr_X_M;
    Tr_u_M_k = [zeros((r-1)*n_u_i,(r-1)*n_u_i*k) eye((r-1)*n_u_i) zeros((r-1)*n_u_i,(r-1)*n_u_i*(N-1-k))]*Tr_U_M;

    if k > 0
        Tr_x_i_k = [zeros(r*n_x_i,r*n_x_i*(k-1)) eye(r*n_x_i) zeros(r*n_x_i,r*n_x_i*(N-k))]*Tr_X_i;
        Tr_x_1_k = [zeros((r-1)*n_x_i,(r-1)*n_x_i*(k-1)) eye((r-1)*n_x_i) zeros((r-1)*n_x_i,(r-1)*n_x_i*(N-k))]*Tr_X_1;
        Tr_x_M_k = [zeros((r-1)*n_x_i,(r-1)*n_x_i*(k-1)) eye((r-1)*n_x_i) zeros((r-1)*n_x_i,(r-1)*n_x_i*(N-k))]*Tr_X_M;
    end
%% Objective function
    H_i = H_i + Tr_x_i_kp1'*Qk_i*Tr_x_i_kp1 + Tr_u_i_k'*Rk_i*Tr_u_i_k;
    H_1 = H_1 + Tr_x_1_kp1'*Qk_1*Tr_x_1_kp1 + Tr_u_1_k'*Rk_1*Tr_u_1_k;
    H_M = H_M + Tr_x_M_kp1'*Qk_M*Tr_x_M_kp1 + Tr_u_M_k'*Rk_M*Tr_u_M_k;
%% Constraints
%% x(k+1|t)=Ax(k|t)+Bu(k|t)
    if (k > 0)&&(k < N-1)
        % for local variable z^i
        AeqC_i = [AeqC_i; 
        Tr_x_i_kp1 - Ak_i*Tr_x_i_k - Bk_i*Tr_u_i_k
        ];
        % for z^1
        AeqC_1 = [AeqC_1; 
        Tr_x_1_kp1 - Ak_1*Tr_x_1_k - Bk_1*Tr_u_1_k
        ];
        % for z^M
        AeqC_M = [AeqC_M; 
        Tr_x_M_kp1 - Ak_M*Tr_x_M_k - Bk_M*Tr_u_M_k
        ];
    end
%% input constraint
%% state constraint
% first
    % for i
%     AineqC_i = [AineqC_i;
%     Tr_x_i_kp1;
%     -Tr_x_i_kp1;
%     ];
%     bineqC_i = [bineqC_i;
%     xmax_i;
%     -xmax_i;
%     ];
%     % for 1
%     AineqC_1 = [AineqC_1;
%     Tr_x_1_kp1;
%     -Tr_x_1_kp1;
%     ];
%     bineqC_1 = [bineqC_1;
%     xmax_1;
%     -xmax_1;
%     ];
%     % for M
%     AineqC_M = [AineqC_M;
%     Tr_x_M_kp1;
%     -Tr_x_M_kp1;
%     ];
%     bineqC_M = [bineqC_M;
%     xmax_M;
%     -xmax_M;
%     ];
% % second
    if k>0
%         for i
        AineqC_i = [AineqC_i;
        Tr_x_i_kp1 - Tr_x_i_k;
        -Tr_x_i_kp1 + Tr_x_i_k;
        ];
        bineqC_i = [bineqC_i;
        xmax_i;
        -xmin_i;
        ];
%         for 1
        AineqC_1 = [AineqC_1;
        Tr_x_1_kp1 - Tr_x_1_k;
        -Tr_x_1_kp1 + Tr_x_1_k;
        ];
        bineqC_1 = [bineqC_1;
        xmax_1;
        -xmin_1;
        ];
%         for M
        AineqC_M = [AineqC_M;
        Tr_x_M_kp1 - Tr_x_M_k;
        -Tr_x_M_kp1 + Tr_x_M_k;
        ];
        bineqC_M = [bineqC_M;
        xmax_M;
        -xmin_M;
        ];
    end
%% terminal constraint
    if k == N-1
        if(strcmpi(terminal,'constraint'))
            % for i
            AeqC_i = [AeqC_i;
            Tr_x_i_kp1
            ];
            % for 1
            AeqC_1 = [AeqC_1;
            Tr_x_1_kp1
            ];
            % for M
            AeqC_M = [AeqC_M;
            Tr_x_M_kp1
            ];
        end
    elseif(strcmpi(terminal,'cost'))
        H_i = H_i - Tr_x_i_kp1'*Qk_i*Tr_x_i_kp1 + Tr_x_i_kp1'*Sk_i*Tr_x_i_kp1;
    end
end
%%
k = 0;
Tr_u_i_0 = [zeros(r*n_u_i,r*n_u_i*k) eye(r*n_u_i) zeros(r*n_u_i,r*n_u_i*(N-1-k))]*Tr_U_i;
Tr_x_i_1 = [zeros(r*n_x_i,r*n_x_i*k) eye(r*n_x_i) zeros(r*n_x_i,r*n_x_i*(N-1-k))]*Tr_X_i;

Tr_u_1_0 =  [zeros((r-1)*n_u_i,(r-1)*n_u_i*k) eye((r-1)*n_u_i) zeros((r-1)*n_u_i,(r-1)*n_u_i*(N-1-k))]*Tr_U_1;
Tr_x_1_1 =  [zeros((r-1)*n_x_i,(r-1)*n_x_i*k) eye((r-1)*n_x_i) zeros((r-1)*n_x_i,(r-1)*n_x_i*(N-1-k))]*Tr_X_1;

Tr_u_M_0 = [zeros((r-1)*n_u_i,(r-1)*n_u_i*k) eye((r-1)*n_u_i) zeros((r-1)*n_u_i,(r-1)*n_u_i*(N-1-k))]*Tr_U_M;
Tr_x_M_1 = [zeros((r-1)*n_x_i,(r-1)*n_x_i*k) eye((r-1)*n_x_i) zeros((r-1)*n_x_i,(r-1)*n_x_i*(N-1-k))]*Tr_X_M;