%% parameter
%xmax_i = [2 0]'; xmin_i = -xmax;
xmax = repmat([2 inf]',M,1);
xmin = -xmax;
% xxmax = 2; xxmin = -xxmax;
%%
H = zeros((n_x+n_u)*N);
f = zeros((n_x+n_u)*N,1);

AineqC = [];
bineqC = [];
AeqC = [];
beqC = [];

Tr_X = [eye(n_x*N), zeros(n_x*N,n_u*N)];
Tr_U = [zeros(n_u*N,n_x*N), eye(n_u*N)];
% ctype_X = repmat('C', [1, n_x*N]);
% ctype_U = repmat('C', [1, n_u*N]);
% ctype = [ctype_X, ctype_U];
%% for loop
for k = 0:N-1
%% preparation
    Tr_u_k = [zeros(n_u,n_u*k) eye(n_u) zeros(n_u,n_u*(N-1-k))]*Tr_U;
    Tr_x_kp1 = [zeros(n_x,n_x*k) eye(n_x) zeros(n_x,n_x*(N-1-k))]*Tr_X;
    if k > 0
        Tr_x_k = [zeros(n_x,n_x*(k-1)) eye(n_x) zeros(n_x,n_x*(N-k))]*Tr_X;
    end
%% Objective function
    H = H + Tr_x_kp1'*Q*Tr_x_kp1;
    H = H + Tr_u_k'*R*Tr_u_k;
%% Constraints
%% x(k+1|t)=Ax(k|t)+Bu(k|t)
    if k>0
        AeqC = [AeqC; 
        Tr_x_kp1 - A*Tr_x_k - B*Tr_u_k
        ];
        beqC = [beqC;
        zeros(n_x,1);
        ];
    end
%% input constraint
%% state constraint
% first
    AineqC = [AineqC;
    Tr_x_kp1
    ];
    bineqC = [bineqC;
    xmax
    ];
    AineqC = [AineqC;
    -Tr_x_kp1
    ];
    bineqC = [bineqC;
    -xmin
    ];
% second
    if k>0
        AineqC = [AineqC;
        Tr_x_kp1 - Tr_x_k
        ];
        bineqC = [bineqC;
        xmax
        ];
        AineqC = [AineqC;
        -Tr_x_kp1 + Tr_x_k
        ];
        bineqC = [bineqC;
        -xmin
        ];
    end
%% terminal constraint
    if k==N-1
        AeqC = [AeqC;
        Tr_x_kp1
        ];
        beqC = [beqC;
        zeros(n_x,1)
        ];
    end
    %% 終端コスト
%     if k == N-1
%         H = H - Tr_x_kp1'*Q*Tr_x_kp1 + Tr_x_kp1'*S*Tr_x_kp1;
%     end
end
%%
k = 0;
Tr_u_0 = [zeros(n_u,n_u*k) eye(n_u) zeros(n_u,n_u*(N-1-k))]*Tr_U;
Tr_x_1 = [zeros(n_x,n_x*k) eye(n_x) zeros(n_x,n_x*(N-1-k))]*Tr_X; 