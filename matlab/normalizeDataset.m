function T = normalizeDataset(T)
% FUNCTION T = normalizeDataset(T)
%
% *** INPUT ARGUMNETS ***
%	T:  time-series dataset matrix (time-series as columns)
%
% *** OUTPUT ARGUMNETS ***
%	T:  time-series dataset matrix with columns normalized to unit norm (and centered)
%
%
% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

[tlen, N] = size(T);

% First center the dataset
try
	T = T - repmat(sum(T)/tlen, tlen, 1);
catch
	for n = 1:N
		T(:, n) = T(:, n) - sum(T(:,n))/tlen;
	end
end

for n = 1:N
	magT = sqrt(sum(T(:, n).^2));
	if magT > 1e-8
		T(:, n) = T(:, n)/magT;
	else
		T(:, n) = 0;
	end
end
