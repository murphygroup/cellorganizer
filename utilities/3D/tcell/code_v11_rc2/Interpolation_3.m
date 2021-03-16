function Need = Interpolation_3(Need)
if isempty(Need)
    return;
end
% delete those with few time points in the begining
time_index = Need(Need(:, 1)~=0, 5) - Need(1,5) + 1;
if sum(ismember(1:5, time_index)) < 3 || sum(ismember(1:10, time_index)) < 6
    Need(2:end, :) = [];
    return;
end

% remove last frames without predicted coordinates.
Need(Need(:, 5) - Need(1,5) + 1 > time_index(end), :) = [];

Temp_m=find(Need(:,1)==0);
if isempty(Temp_m)
    return;
end

for T_i=1:length(Temp_m)-3
    if Temp_m(T_i) + 3 == Temp_m(T_i + 3)
        Need(Temp_m(T_i):end, :) = [];
        Temp_m(T_i:end) = [];
        break;
    end
end

% linear interporlation for each coordinate
if size(Need, 1) < 2 || numel(Temp_m) < 1
    return;
end
xq = Need(:, 5);
Need_v = Need;
Need_v(Temp_m, :) = [];
x = Need_v(:, 5);
% Need = [];
for i = 1 : 4
    v = Need_v(:, i);
    vq = interp1(x, v, xq);
    Need(:, i) = vq;   
end

end

