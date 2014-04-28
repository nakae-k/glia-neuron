function [spike_mat glia_mat] = cross(spike, glia, para,n)
%n = 10;

[dim_n, T] = size(spike); [dim_g, T] = size(glia);
T_unit = fix(T./n);
spike(:,(T_unit*n)+1:T) = []; spike = sparse(spike);
glia(:,(T_unit*n)+1:T) = [];

est_spike = cell(1,n); test_spike = cell(1,n);
est_glia = cell(1,n);  test_glia = cell(1,n);

for i = 1:n
    test_spike{i} = sparse( spike(:, ((i-1)*T_unit+1) : (i*T_unit) ) );
    temp_spike = spike; temp_spike(:, ((i-1)*T_unit+1) : (i*T_unit) ) = [];
    est_spike{i}  = sparse(temp_spike);

    test_glia{i} = sparse( glia(:, ((i-1)*T_unit+1) : (i*T_unit) ) );
    temp_glia = glia; temp_glia(:, ((i-1)*T_unit+1) : (i*T_unit) ) = [];
    est_glia{i}  = sparse(temp_glia);
end

spike_mat = zeros(dim_n,n); glia_mat = zeros(dim_g,n);

for i = 1:n
    i
    [W B S] = est_sg(est_spike{i},est_glia{i}, para);
    [v_spike v_glia] = val_kull(test_spike{i}, test_glia{i}, W, B, S, para);
    spike_mat(:,i) = v_spike; glia_mat(:,i) = v_glia;

end

end