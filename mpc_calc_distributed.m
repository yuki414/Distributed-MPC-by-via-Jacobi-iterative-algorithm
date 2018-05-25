%% Optimization
H_i = sparse((H_i + H_i')/2);
H_1 = sparse((H_1 + H_1')/2);
H_M = sparse((H_M + H_M')/2);
%% constraints
% AineqV = [];   bineqV = [];
AeqV_i = [];    beqV_i = [];
AeqV_1 = [];    beqV_1 = [];
AeqV_M = [];    beqV_M = [];
%% initial condition
x_i = x(2*i-3:2*i+2,1);
x_1 = x(1:4,1);
x_M = x(M-3:M,1);
% if i > 1
%     x_im1 = x(2*(i-1)-1:2*(i-1),1);
%     x_im2 = zeros(2,1);
% end
% if i < M
%     x_ip1 = x(2*(i+1)-1:2*(i+1),1);
%     x_ip2 = zeros(2,1);
% end
% if i > 2
%     x_im2 = x(2*(i-2)-1:2*(i-2),1);
% end
% if i < M - 1
%     x_ip2 = x(2*(i+2)-1:2*(i+2),1);
% end
%% x(k+1|t)=A*x(k|t)+B*u(k|t)
switch i
%% for i
    case 1
        AeqV_1 =[AeqV_1;
        Tr_x_1_1 - Bk_1*Tr_u_1_0
        ];
        beqV_1 = [beqV_1;
        Ak_1*x_1
        ];
%% for M
    case M
        AeqV_M =[AeqV_M;
        Tr_x_M_1 - Bk_M*Tr_u_M_0
        ];
        beqV_M = [beqV_M;
        Ak_M*x_M
        ];
%% for i
    otherwise
        AeqV_i =[AeqV_i;
        Tr_x_i_1 - Bk_i*Tr_u_i_0
        ];
        beqV_i = [beqV_i;
        Ak_i*x_i
        ];
        for k = 0:N-1
            beqC_i = [beqC_i;
            % constant value corresponded to A^-i
            
            ];
        end
end
%% constraint matrix
switch i
    case M
        H = H_M;
        f = f_M;
        Aeq = [AeqV_M ; AeqC_M];
        beq = [beqV_M ; beqC_M];
        Aineq = [AineqV_M ; AineqC_M];
        bineq = [bineqV_M ; bineqC_M];
        lb = ones(1,(r-1)*(n_u_i+n_x_i)*N);
        ub = ones(1,(r-1)*(n_u_i+n_x_i)*N);
    case 1
        H = H_1;
        f = f_1;
        Aeq = [AeqV_1 ; AeqC_1];
        beq = [beqV_1 ; beqC_1];
        Aineq = [AineqV_1 ; AineqC_1];
        bineq = [bineqV_1 ; bineqC_1];
        lb = ones(1,(r-1)*(n_u_i+n_x_i)*N);
        ub = ones(1,(r-1)*(n_u_i+n_x_i)*N);
    otherwise
        H = H_i;
        f = f_i;
        Aeq = [AeqV_i ; AeqC_i];
        beq = [beqV_i ; beqC_i];
        Aineq = [AineqV_i ; AineqC_i];
        bineq = [bineqV_i ; bineqC_i];
        lb = ones(1,r*(n_u_i+n_x_i)*N);
        ub = ones(1,r*(n_u_i+n_x_i)*N);
end
%% cplex
if(strcmpi(solver,'cplex'))
    options = cplexoptimset; options.Diagnostics = 'off'; options.Display = 'off';
    [z_i_p, fval, exitflag, output1] = cplexqp(H, f, Aineq, bineq, Aeq, beq, [], [], [], options);
    if isempty(z_i_p)~=0
        disp(['no optimal answers : exitflag = ', num2str(exitflag)])
        break;
    end
end
%% gurobi
if(strcmpi(solver,'gurobi'))
    params.outputflag = 0;
    model.modelsense = 'min'; % Å¬‰»‚©Å‘å‰»‚©
    model.Q = H;
    model.obj = f; % 
%     model.vtype =  ctype;
    model.lb = -inf*lb;
    model.ub = inf*ub;
    model.sense = [repmat('<',[1,length(bineq)]) repmat('=',[1,length(beq)])];
    model.A = sparse([Aineq;Aeq]);
    model.rhs = [bineq;beq];

    z_i_opt = gurobi(model, params); % gurobi ‚ÅŒvŽZ
    if strcmpi(z_i_opt.status,'optimal')
        z_i_p = z_i_opt.x;
    else
        disp(['no optimal answers : Uopt.status = ', z_i_opt.status])
        t = t - 1;
        break;
    end
end
%%
% x_i_1 = Tr_x_i_1_d(:,:,i)*z_i_p;
% u_i_0 = Tr_u_i_0_d(:,:,i)*z_i_p;