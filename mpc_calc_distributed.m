%% Optimization
H_i = sparse((H_i + H_i')/2);
%% constraints
AineqV = [];   bineqV = [];
AeqV = [];      beqV = [];
%% initial condition
x_i = x0(2*i-1:2*i,1);
if i > 1
    x_im1 = x0(2*(i-1)-1:2*(i-1),1);
end
if i < M
    x_ip1 = x0(2*(i+1)-1:2*(i+1),1);
end
%% x(k+1|t)=A*x(k|t)+B*u(k|t)
switch i
    case 1
        AeqV =[AeqV;
        Tr_x_i_1_d(:,:,i) - Bii*Tr_u_i_0_d(:,:,i)
        ];
        beqV = [beqV;
        Aii*x_i + Aip1*x_ip1;
        ];
    case M
        AeqV =[AeqV;
        Tr_x_i_1_d(:,:,i) - Bii*Tr_u_i_0_d(:,:,i)
        ];
        beqV = [beqV;
        Aii*x_i + Aim1*x_im1;
        ];
    otherwise
        AeqV =[AeqV;
        Tr_x_i_1_d(:,:,i) - Bii*Tr_u_i_0_d(:,:,i)
        ];
        beqV = [beqV;
        Aii*x_i + Aim1*x_im1 + Aip1*x_ip1;
        ];
end
%% constraint matrix
Aineq = [AineqV ; AineqC_d(:,:,i)];
bineq = [bineqV ; bineqC_d(:,:,i)];
Aeq = [AeqV ; AeqC_d(:,:,i)];
beq = [beqV ; beqC_d(:,:,i)];
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
    model.modelsense = 'min'; % Å¬‰»‚©Å‘å‰»‚©
    model.Q = H_i;
    model.obj = f_i; % 
%     model.vtype =  ctype;
    model.lb = -inf*ones(1,(n_u+n_x)*N);
    model.ub = inf*ones(1,(n_u+n_x)*N);
    model.sense = [repmat('<',[1,length(bineq)]) repmat('=',[1,length(beq)])];
    model.A = sparse([Aineq;Aeq]);
    model.rhs = [bineq;beq];

    z_i_opt = gurobi(model, params); % gurobi ‚ÅŒvZ
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