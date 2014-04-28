function [W B Sigma] = est_sg(spike, glia, para);

%spike = Y; glia = Zd;
lambda = 1;

[dim_n, T] = size(spike); [dim_g, T] = size(glia);

 h_nn = para(1,1); h_ng = para(1,2);
 h_gn = para(2,1); h_gg = para(2,2);
 lambda_nn = para(3,1); lambda_ng = para(3,2);
 lambda_gn = para(4,1); lambda_gg = para(4,2);


% h_nn = 10; h_ng = 10;
% h_gn = 10; h_gg = 10;
% lambda_nn = 10; lambda_ng = 10;
% lambda_gn = 100; lambda_gg = 10;

theta = zeros(dim_n,1); theta_g = zeros(dim_g,1);
A_nn = zeros(dim_n, dim_n, h_nn); A_ng = zeros(dim_n, dim_g, h_ng);
A_gn = zeros(dim_g, dim_n, h_gn); A_gg = zeros(dim_g, dim_g, h_gg);

X = matrixX(spike, glia, [dim_n dim_g h_nn h_ng T]);

if h_nn == 0
    S = 0;
else
    S = matrixR([lambda_nn lambda_ng dim_n dim_g h_nn h_ng T]);
end

p = dim_n*h_nn+ dim_g*h_ng +1;
W = zeros(p,dim_n);

parfor i = 1:dim_n
    w = zeros(p,1); y = spike(i,:)';
    f =@(w)var_ridge(w,y,X,lambda,S);
    options = []; options.display = 'off';
    w = minFunc(f, w, options);
    W(:,i) = w;

    A_nn(i,:,:) = reshape(w(1:dim_n*h_nn), dim_n, h_nn);
    A_ng(i,:,:) = reshape(w(dim_n*h_nn+1:dim_n*h_nn+dim_g*h_ng), dim_g, h_ng);
    theta(i) = w(p);
end

X = matrixX(spike, glia, [dim_n dim_g h_gn h_gg T]);
if h_gg == 0
    S = 0;
else
    S = matrixR([lambda_gn lambda_gg dim_n dim_g h_gn h_gg T]);
end

q = dim_n*h_gn+ dim_g*h_gg +1;
D = (X*X' + lambda*ones(q,q) + S'*S);
B = D \ (X*glia'); B = B';

temp = glia' - X'*B';
Sigma = temp'*temp/T;
%size(S)

for i = 1:dim_g    
    A_gn(i, :, :) = reshape(B(i,1:dim_n*h_gn), dim_n, h_gn);
    A_gg(i, :, :) = reshape(B(i,dim_n*h_gn+1:dim_n*h_gn+dim_g*h_gg), dim_g, h_gg);
    theta_g(i) = B(i,q);
end

% W = cell(2,3);
% W{1,1} = A_nn; W{1,2} = A_ng; W{1,3} = theta;
% W{2,1} = A_gn; W{2,2} = A_gg; W{2,3} = theta_g;

end