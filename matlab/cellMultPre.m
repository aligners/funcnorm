function C = cellMultPre(A, B, keepAsCell)
if nargin == 2
	keepAsCell = 0;
end

C = cellfun(@funMultPre, A, B, 'UniformOutput', false);
if keepAsCell == 0
	C = cell2mat(C);
end

function C = funMultPre(A, B)
C = A*B;
