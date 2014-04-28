function [f, g] = logistic_loss(w,y,X)

z = X'*w; p = exp(z)./(1+exp(z));

f = -y'*z + sum(log(1+exp(z)));
g = X*(-y+p);

end