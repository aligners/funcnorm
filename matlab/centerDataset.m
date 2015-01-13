function T = centerDataset(T)

[tlen, N] = size(T);

try
    % Sometimes this produces an out-of-memory error
    T = T - repmat(sum(T)/tlen, tlen, 1);
catch
    for n = 1:N
        T(:, n) = T(:, n) - sum(T(:, n))/tlen;
    end
end