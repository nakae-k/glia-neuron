function [f,g] = var_ridge_lam(w,y,X,lambda,S,p_n)

    [f, g] = logistic_loss(w,y,X);
    
    lambda_n = lambda(1); lambda_g = lambda(2);
    
    tempw = [w(1:end-1) ; 0];    
    temp_n = zeros(size(tempw));
    temp_n(1:p_n) = w(1:p_n);
    temp_g = zeros(size(tempw));
    temp_g(p_n+1:end-1) = w(p_n+1:end-1);
    
    f = f + 0.5*lambda_n*norm(temp_n, 2)^2 + 0.5*lambda_g*norm(temp_g, 2)^2 + 0.5*norm(S*w, 2)^2;
    g = g + lambda_n*temp_n + lambda_g*temp_g + S'*(S*w);
end
