%% parameter
n_x_i = 2;  % i番目の状態のサイズ
n_u_i = 1;
% xmax_i = [zeros(n_x_i*(i-1),1); 2; inf; zeros(n_x_i*(M-i),1)];
xmax_i = [2, inf]';
xmin_i = -xmax_i;
%% setting
% Ai = [Aim1; Aii; Aip1]; % i番目に対応したカリグラフA
% Bi = [zeros(n_u_i*(i-1),n_u_i); Bii; zeros(n_u_i*(M-i),n_u_i)]; % カリグラフB

H_i = zeros((n_x+n_u)*N);
f_i = zeros((n_x+n_u)*N,1);

AineqC = [];
bineqC = [];
AeqC = [];
beqC = [];

Tr_X = [eye(n_x*N), zeros(n_x*N,n_u*N)];
Tr_U = [zeros(n_u*N,n_x*N), eye(n_u*N)];
Tr_X_i = [zeros(n_x_i*N,n_x_i*(i-1)*N), eye(n_x_i*N), zeros(n_x_i*N,n_x_i*(M-i)*N), zeros(n_x_i*N,n_u*N)];
Tr_U_i = [zeros(n_u_i*N,n_x*N), zeros(n_u_i*N,n_u_i*(i-1)*N), eye(n_u_i*N), zeros(n_u_i*N,n_u_i*(M-i)*N)];
if i > 1
    Tr_X_im1 = [zeros(n_x_i*N,n_x_i*(i-2)*N), eye(n_x_i*N), zeros(n_x_i*N,n_x_i*(M-i+1)*N), zeros(n_x_i*N,n_u*N)];
end
if i < M
    Tr_X_ip1 = [zeros(n_x_i*N,n_x_i*(i)*N), eye(n_x_i*N), zeros(n_x_i*N,n_x_i*(M-i-1)*N), zeros(n_x_i*N,n_u*N)];
end
%% for loop
for k = 0:N-1
%% preparation for i
    Tr_x_kp1 = [zeros(n_x,n_x*k) eye(n_x) zeros(n_x,n_x*(N-1-k))]*Tr_X;
    Tr_u_k = [zeros(n_u,n_u*k) eye(n_u) zeros(n_u,n_u*(N-1-k))]*Tr_U;
    Tr_u_i_k = [zeros(n_u_i,n_u_i*k) eye(n_u_i) zeros(n_u_i,n_u_i*(N-1-k))]*Tr_U_i;
    Tr_x_i_kp1 = [zeros(n_x_i,n_x_i*k) eye(n_x_i) zeros(n_x_i,n_x_i*(N-1-k))]*Tr_X_i;
    if k > 0
        Tr_x_i_k = [zeros(n_x_i,n_x_i*(k-1)) eye(n_x_i) zeros(n_x_i,n_x_i*(N-k))]*Tr_X_i;
    end
%% preparation for i-1 and i+1
    if i > 1
%     Tr_u_i_k = [zeros(n_u_i,n_u_i*k) eye(n_u_i) zeros(n_u_i,n_u_i*(N-1-k))]*Tr_U_i;
%         Tr_x_im1_kp1 = [zeros(n_x_i,n_x_i*k) eye(n_x_i) zeros(n_x_i,n_x_i*(N-1-k))]*Tr_X_im1;
        if k > 0
            Tr_x_im1_k = [zeros(n_x_i,n_x_i*(k-1)) eye(n_x_i) zeros(n_x_i,n_x_i*(N-k))]*Tr_X_im1;
        end
    end
    if i < M
%     Tr_u_i_k = [zeros(n_u_i,n_u_i*k) eye(n_u_i) zeros(n_u_i,n_u_i*(N-1-k))]*Tr_U_i;
%         Tr_x_ip1_kp1 = [zeros(n_x_i,n_x_i*k) eye(n_x_i) zeros(n_x_i,n_x_i*(N-1-k))]*Tr_X_ip1;
        if k > 0
            Tr_x_ip1_k = [zeros(n_x_i,n_x_i*(k-1)) eye(n_x_i) zeros(n_x_i,n_x_i*(N-k))]*Tr_X_ip1;
        end
    end    
%% Objective function
    H_i = H_i + Tr_x_kp1'*Q*Tr_x_kp1;
    H_i = H_i + Tr_u_k'*R*Tr_u_k;
%     H_i = H_i + Tr_x_i_kp1'*Qi*Tr_x_i_kp1;
%     H_i = H_i + Tr_u_i_k'*Ri*Tr_u_i_k;
%% Constraints
%% x(k+1|t)=Ax(k|t)+Bu(k|t)
    switch i
        case 1
            if k > 0
                AeqC = [AeqC; 
                Tr_x_i_kp1 - Aii*Tr_x_i_k - Aip1*Tr_x_ip1_k - Bii*Tr_u_i_k
                ];
                beqC = [beqC;
                zeros(n_x_i,1);
                ];
            end
        case M
            if k > 0
                AeqC = [AeqC; 
                Tr_x_i_kp1 - Aii*Tr_x_i_k - Aim1*Tr_x_im1_k  - Bii*Tr_u_i_k
                ];
                beqC = [beqC;
                zeros(n_x_i,1);
                ];
            end
        otherwise
            if k > 0
                AeqC = [AeqC; 
                Tr_x_i_kp1 - Aii*Tr_x_i_k - Aim1*Tr_x_im1_k - Aip1*Tr_x_ip1_k - Bii*Tr_u_i_k
                ];
                beqC = [beqC;
                zeros(n_x_i,1);
                ];
            end
    end
%% input constraint
%% state constraint
% first
    AineqC = [AineqC;
    Tr_x_i_kp1
    ];
    bineqC = [bineqC;
    xmax_i
    ];
    AineqC = [AineqC;
    -Tr_x_i_kp1
    ];
    bineqC = [bineqC;
    -xmin_i
    ];
% second
    if k>0
        AineqC = [AineqC;
        Tr_x_i_kp1 - Tr_x_i_k
        ];
        bineqC = [bineqC;
        xmax_i
        ];
        AineqC = [AineqC;
        -Tr_x_i_kp1 + Tr_x_i_k
        ];
        bineqC = [bineqC;
        -xmin_i
        ];
    end
%% terminal constraint
    if k==N-1
        AeqC = [AeqC;
        Tr_x_i_kp1
        ];
        beqC = [beqC;
        zeros(n_x_i,1)
        ];
    end
    %% 終端コスト
%     if k == N-1
%         H = H - Tr_x_i_kp1'*Q*Tr_x_i_kp1 + Tr_x_i_kp1'*S*Tr_x_i_kp1;
%     end
end
%%
k = 0;
Tr_u_i_0 = [zeros(n_u_i,n_u_i*k) eye(n_u_i) zeros(n_u_i,n_u_i*(N-1-k))]*Tr_U_i;
Tr_x_i_1 = [zeros(n_x_i,n_x_i*k) eye(n_x_i) zeros(n_x_i,n_x_i*(N-1-k))]*Tr_X_i;
% if i > 1
%     Tr_x_im1_1 = [zeros(n_x_i,n_x_i*k) eye(n_x_i) zeros(n_x_i,n_x_i*(N-1-k))]*Tr_X_im1;
% end
% if i < M
%     Tr_x_ip1_1 = [zeros(n_x_i,n_x_i*k) eye(n_x_i) zeros(n_x_i,n_x_i*(N-1-k))]*Tr_X_ip1;
% end
%% 多次元配列に格納
% dataD = [];
H_d(:,:,i) = H_i;
AeqC_d(:,:,i) = AeqC;
beqC_d(:,:,i) = beqC;
AineqC_d(:,:,i) = AineqC;
bineqC_d(:,:,i) = bineqC;
Tr_u_i_0_d(:,:,i) = Tr_u_i_0;
Tr_x_i_1_d(:,:,i) = Tr_x_i_1;
%% i番目のシステムに対して名前付け
% eval(['H_' num2str(i) ' = H_i;'])
% eval(['AeqC_' num2str(i) ' = AeqC;'])
% eval(['beqC_' num2str(i) ' = beqC;'])
% eval(['Aineq_' num2str(i) ' = AineqC;'])
% eval(['bineq_' num2str(i) ' = bineqC;'])
% eval(['Aineq_' num2str(i) ' = H;']);
% eval(['Tr_u_' num2str(i) '_0 = Tr_u_i_0;'])
% eval(['Tr_x_' num2str(i) '_1 = Tr_x_i_1;'])
% if i > 1
%     eval(['Tr_x_' num2str(i) '_m1 = Tr_x_im1_1;']);
% end
% if i < M
%     eval(['Tr_x_' num2str(i) '_p1 = Tr_x_ip1_1;']);
% end