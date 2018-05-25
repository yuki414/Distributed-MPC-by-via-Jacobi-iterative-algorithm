%% Optimization
H = sparse((H + H')/2);
%% constraints
AineqV = [];   bineqV = [];
AeqV = [];      beqV = [];
%% x(k+1|t)=A*x(k|t)+B*u(k|t)
AeqV =[AeqV;
    Tr_x_1 - B*Tr_u_0
    ];
beqV = [beqV;
     A*x_c;
    ];
%% constraint matrix
Aineq = [AineqV ; AineqC];
bineq = [bineqV ; bineqC];
Aeq = [AeqV ; AeqC];
beq = [beqV ; beqC];
%% cplex
if(strcmpi(solver,'cplex'))
    options = cplexoptimset; options.Diagnostics = 'off'; options.Display = 'off';
    [Usol, fval, exitflag, output1] = cplexqp(H, f, Aineq, bineq, Aeq, beq, [], [], [], options);
    if isempty(Usol)~=0
        disp(['no optimal answers : exitflag = ', num2str(exitflag)])
        break;
    end
end
%% gurobi
if(strcmpi(solver,'gurobi'))
    params.outputflag = 0;
    model.modelsense = 'min'; % ç≈è¨âªÇ©ç≈ëÂâªÇ©
    model.Q = H;
    model.obj = f; % 
%     model.vtype =  ctype;
    model.lb = -inf*ones(1,(n_u+n_x)*N);
    model.ub = inf*ones(1,(n_u+n_x)*N);
    model.sense = [repmat('<',[1,length(bineq)]) repmat('=',[1,length(beq)])];
    model.A = sparse([Aineq;Aeq]);
    model.rhs = [bineq;beq];

    Uopt = gurobi(model, params); % gurobi Ç≈åvéZ
    if strcmpi(Uopt.status,'optimal')
        Usol = Uopt.x;
    else
        disp(['no optimal answers : Uopt.status = ', Uopt.status])
        t = t - 1;
        break;
    end
end
%%
u_c = Tr_u_0*Usol;