function [X] = matrixX(spike, glia, temp)
%UNTITLED この関数の概要をここに記述
%   詳細説明をここに記述

dim_n = temp(1); dim_g = temp(2); h_nn = temp(3); h_ng = temp(4); T = temp(5);
p = dim_n*h_nn+ dim_g*h_ng +1;
X = sparse(p, T);

for s = 1:h_nn
    pspike = [sparse(zeros(dim_n, s)) spike];
    pspike(:,T+1:T+s) = [];
    X(dim_n*(s-1)+1:dim_n*s, :) = pspike;
end
for s = 1:h_ng
    pglia = [sparse(zeros(dim_g, s)) glia];    
    pglia(:,T+1:T+s) = [];
    X(dim_n*h_nn + dim_g*(s-1)+1:dim_n*h_nn + dim_g*s, :) = pglia;
end
X(p,:) = ones(1,T);

end

