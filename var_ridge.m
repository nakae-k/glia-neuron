function [f,g] = var_ridge(w,y,X,lambda,S)

    [f, g] = logistic_loss(w,y,X);
    
    tempw = [w(1:end-1) ; 0];   
    
    f = f + 0.5*lambda*norm(tempw, 2)^2 + 0.5*norm(S*w, 2)^2;
    g = g + lambda*tempw + S'*(S*w);
end
