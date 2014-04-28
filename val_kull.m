function [v_spike v_glia] = val_kull(spike, glia, W, B, Sigma, para);

[dim_n, T] = size(spike); [dim_g, T] = size(glia);

 h_nn = para(1,1); h_ng = para(1,2);
 h_gn = para(2,1); h_gg = para(2,2);
 lambda_nn = para(3,1); lambda_ng = para(3,2);
 lambda_gn = para(4,1); lambda_gg = para(4,2);

X = matrixX(spike, glia, [dim_n dim_g h_nn h_ng T]);
p = dim_n*h_nn+ dim_g*h_ng +1;

v_spike = zeros(dim_n,1);

for i = 1:dim_n
    w = W(:,i); y = spike(i,:)';
    v_spike(i) = -1*logistic_loss(w, y, X);
end

X = matrixX(spike, glia, [dim_n dim_g h_gn h_gg T]);
v_glia = zeros(dim_g,1);

for i = 1:dim_g
    w = B(i,:)'; y = glia(i,:)';
    v_glia(i) = -0.5*T*log(Sigma(i,i)) - (0.5/Sigma(i,i))*norm(y-X'*w,2)^2;
end

end