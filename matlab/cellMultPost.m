function C = cellMultPost(A, B, keepAsCell)
if nargin == 2
	keepAsCell = 0;
end

C = cellfun(@funMultPost, A, B, 'UniformOutput', false);
if keepAsCell == 0
	C = cell2mat(C);
end

function C = funMultPost(A, B)
C = A*B;
