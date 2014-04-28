function [S] = matrixR(temp)
%UNTITLED2 この関数の概要をここに記述
%   詳細説明をここに記述
lambda_nn = temp(1);
lambda_ng = temp(2);
dim_n = temp(3);
dim_g = temp(4);
h_nn = temp(5);
h_ng = temp(6);
T = temp(7);

if h_nn == 0
Q_g = spdiag(sqrt(lambda_ng)*ones(dim_g*(h_ng-1),1)); R_g = [Q_g sparse(dim_g*(h_ng-1), dim_g)] + [sparse(dim_g*(h_ng-1), dim_g) -Q_g];
S = [R_g sparse(size(R_g,1), 1);
    sparse(1,size(R_g,2)+1)];
elseif h_ng == 0
Q = spdiag(sqrt(lambda_nn)*ones(dim_n*(h_nn-1),1));   R = [Q sparse(dim_n*(h_nn-1), dim_n)] + [sparse(dim_n*(h_nn-1), dim_n) -Q];
S = [R sparse(size(R,1), 1);
    sparse(1,size(R,2)+1)];
else
Q = spdiag(sqrt(lambda_nn)*ones(dim_n*(h_nn-1),1));   R = [Q sparse(dim_n*(h_nn-1), dim_n)] + [sparse(dim_n*(h_nn-1), dim_n) -Q];
Q_g = spdiag(sqrt(lambda_ng)*ones(dim_g*(h_ng-1),1)); R_g = [Q_g sparse(dim_g*(h_ng-1), dim_g)] + [sparse(dim_g*(h_ng-1), dim_g) -Q_g];
S =  [R                            sparse(size(R,1), size(R_g,2)+1);
     sparse(size(R_g,1),size(R,2)) R_g   sparse(size(R_g,1),1);
     sparse(1,size(R,2)+size(R_g,2)+1)];
end

end

