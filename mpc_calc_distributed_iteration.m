%% Optimization
H_i = sparse((H_i + H_i')/2);
%% constraints
AineqV = [];   bineqV = [];
AeqV = [];      beqV = [];
%% initial condition
x_i = x(2*i-1:2*i,1);
if i > 1
    x_im1 = x(2*(i-1)-1:2*(i-1),1);
    x_im2 = zeros(2,1);
end
if i < M
    x_ip1 = x(2*(i+1)-1:2*(i+1),1);
    x_ip2 = zeros(2,1);
end
if i > 2
    x_im2 = x(2*(i-2)-1:2*(i-2),1);
end
if i < M - 1
    x_ip2 = x(2*(i+2)-1:2*(i+2),1);
end
%% x(k+1|t)=A*x(k|t)+B*u(k|t)
switch i
    case 1
        AeqV =[AeqV;
        Tr_x_i_1_d(:,:,i) - Bii*Tr_u_i_0_d(:,:,i)
        ];
        beqV = [beqV;
        Aii*x_i + Aip1*x_ip1
        ];
        AeqV = [AeqV;
        Tr_x_ip1_1_d(:,:,i) - Bii*Tr_u_ip1_0_d(:,:,i) 
        ];
        beqV = [beqV;
        Aii*x_ip1 + Aim1*x_i 
        ];
    case M
        AeqV =[AeqV;
        Tr_x_i_1_d(:,:,i) - Bii*Tr_u_i_0_d(:,:,i)
        ];
        beqV = [beqV;
        Aii*x_i + Aim1*x_im1;
        ];
        AeqV =[AeqV;
        Tr_x_im1_1_d(:,:,i) - Bii*Tr_u_im1_0_d(:,:,i)
        ];
        beqV = [beqV;
        Aii*x_im1 + Aip1*x_i;
        ];
    otherwise
        AeqV =[AeqV;
        Tr_x_i_1_d(:,:,i) - Bii*Tr_u_i_0_d(:,:,i)
        ];
        beqV = [beqV;
        Aii*x_i + Aim1*x_im1 + Aip1*x_ip1;
        ];
        AeqV = [AeqV;
        Tr_x_ip1_1_d(:,:,i) - Bii*Tr_u_ip1_0_d(:,:,i) 
        ];
        beqV = [beqV;
        Aii*x_ip1 + Aim1*x_i + Aim1*x_im2
        ];
        AeqV =[AeqV;
        Tr_x_im1_1_d(:,:,i) - Bii*Tr_u_im1_0_d(:,:,i)
        ];
        beqV = [beqV;
        Aii*x_im1 + Aip1*x_i + Aip1*x_ip2;
        ];
        
end
%% constraint matrix
switch i
    case M
        Aeq = [AeqV ; AeqC_d_M(:,:,i)];
        beq = [beqV ; beqC_d_M(:,:,i)];
        Aineq = [AineqV ; AineqC_d_M(:,:,i)];
        bineq = [bineqV ; bineqC_d_M(:,:,i)];
    case 1
        Aeq = [AeqV ; AeqC_d_1(:,:,i)];
        beq = [beqV ; beqC_d_1(:,:,i)];
        Aineq = [AineqV ; AineqC_d_1(:,:,i)];
        bineq = [bineqV ; bineqC_d_1(:,:,i)];
    otherwise
        Aeq = [AeqV ; AeqC_d(:,:,i)];
        beq = [beqV ; beqC_d(:,:,i)];
        Aineq = [AineqV ; AineqC_d(:,:,i)];
        bineq = [bineqV ; bineqC_d(:,:,i)];
end

%% cplex
if(strcmpi(solver,'cplex'))
    options = cplexoptimset; options.Diagnostics = 'off'; options.Display = 'off';
    [z_i_p, fval, exitflag, output1] = cplexqp(H_d(:,:,i), f_i, Aineq, bineq, Aeq, beq, [], [], [], options);
    if isempty(z_i_p)~=0
        disp(['no optimal answers : exitflag = ', num2str(exitflag)])
        break;
    end
end
%% gurobi
if(strcmpi(solver,'gurobi'))
    params.outputflag = 0;
    model.modelsense = 'min'; % 最小化か最大化か
    model.Q = H_i;
    model.obj = f_i; % 
%     model.vtype =  ctype;
    model.lb = -inf*ones(1,(n_u+n_x)*N);
    model.ub = inf*ones(1,(n_u+n_x)*N);
    model.sense = [repmat('<',[1,length(bineq)]) repmat('=',[1,length(beq)])];
    model.A = sparse([Aineq;Aeq]);
    model.rhs = [bineq;beq];

    z_i_opt = gurobi(model, params); % gurobi で計算
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