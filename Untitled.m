N=3;
M=5;
X=[1:N*M*2]';
U=[1:N*M]';
X_i = zeros(2*N,M);
U_i = zeros(N,M);
for i = 1:M
    for k = 1:N
    X_i(n_x_i*k-1:n_x_i*k,i) = X(n_x_i*(k-1)*M+(i-1)*n_x_i+1:n_x_i*(k-1)*M+(i-1)*n_x_i+2,1);
    U_i(n_u_i*k,i) = U(n_u_i*k+(i-1)*N*n_u_i,1);
    end
end

xx = [];
uu = [];
for i = 1:M
    for k = 1:N
%         X_i(n_x_i*k-1:n_x_i*k,i)
        xx(n_x_i*(k-1)*M+(i-1)*n_x_i+1:n_x_i*(k-1)*M+(i-1)*n_x_i+2,1)=X_i(n_x_i*k-1:n_x_i*k,i);
%         U_i(n_u_i*k,i)
        uu(n_u_i*k+(i-1)*N*n_u_i,1)=U_i(n_u_i*k,i);
    end
end