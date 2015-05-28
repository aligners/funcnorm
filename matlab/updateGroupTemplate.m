function [VG, sg, hasErrored]=updateGroupTemplate(VG, sg, ng, VA, sa, VS, ss)
% FUNCTION [V,s]=updateGroupTemplate(VG, sg, ng, VA, sa [, VS, ss])
% (VG, sg) represent the eigenvectors/sqrt(eigenvalues) of the previous group
%   template
% ng is the number of subjects averaged in the previous 
%    group template (VG, sg)
% ncomponents is the number of eigenvector/sqrt(eigenvalue) components to compute
% (VA, sa) represent the eigenvectors/sqrt(eigenvalues) of the subject to be
%   added into the template
% (VS, ss) represent the eigenvectors/sqrt(eigenvalues) of the subject to be
%   removed from the template (optional)

hasErrored = 0;
tic;
ticSum = 0;

if nargin == 5
    % Then we are only adding a subject in
    mode = 'ADD';

    sg = (ng/(ng+1))*sg.^2;
    sa = (1/(ng+1))*sa.^2;
elseif nargin == 7
    % Then we are adding one in, taking one out
    mode = 'ADD_SUBTRACT';

    % Leave sg alone!
    sg = sg.^2;
    sa = (1/ng)*sa.^2;
    ss = (1/ng)*ss.^2;
else
    error('updateGroupTemplate:  See prototype by typing ''help updateGroupTempate''');
end
clear ng;

N = size(VG, 1);

if strcmp(mode, 'ADD')
    % Compute projection matrices for subject to be added
    PA = VG'*VA;
    [XA, RA] = orth_qr(VA - VG*PA);

    PS = 0;
    ss = 0;
elseif strcmp(mode, 'ADD_SUBTRACT')
    PA = VG'*VA;
    [XA, RA] = orth_qr(VA - VG*PA);
    
    % Compute projection matrices for subject to be removed
    PS = VG'*VS;
end

% Compute block matrix components
A11 = diag(sg) + diagMultPost(PA, sa)*PA' - diagMultPost(PS, ss)*PS';
A21 = diagMultPost(RA, sa)*PA';
A12 = A21';
A22 = diagMultPost(RA, sa)*RA';
A = [A11 A12;
     A21 A22];


ticSum = ticSum + toc;
disp(['Completed setup phase']);
disp(['Time:  ', num2str(ticSum), ' seconds']);
pause(0.05);

tic;
% Changing sysev to eig to prevent errors
%[V,D] = syev(A, 'mrrr');
[V,D] = eig(A);
s = sqrt(diag(D));
[s, I] = sort(s, 1, 'descend');
V = V(:, I);
s = real(s); % Sometimes there are negative values in D (that become imaginary upon sqrt)
clear D;
ticSum = ticSum + toc;
disp(['Completed syev']);
disp(['Time:  ', num2str(ticSum), ' seconds']);
pause(0.05);

clear A A11 A12 A21 A22 PA PS RA VA sa;
tic;

try
	VG = [VG XA]*V;
catch exception
	disp('Couldn''t do it the first time');
	% Sometimes the above doesn't work because of out of memory problems
	% So we can overwrite VG one at a time
	try
		VG(:, (size(VG,2)+1):(size(VG,2)+size(XA,2))) = XA;
	catch
		disp('Must write out to file because there is not enough memory...');
		% We're going to have to do it the hard way...
		tmpFile = tempname;%'uGT_tmp.bin';
		fp = fopen(tmpFile, 'wb');
		fwrite(fp, VG, 'single');
		fwrite(fp, XA, 'single');
		fclose(fp);
		vg1 = size(VG,1);
		vg2 = size(VG,2) + size(XA,2);
		clear VG XA;
		fp = fopen(tmpFile, 'rb');
		try
			VG = fread(fp, [vg1, vg2], 'single=>single');
		catch
			fclose(fp);
			system(['rm -f ', tmpFile]);
			hasErrored = 1;
			whos
			return;
		end
		fclose(fp);
		system(['rm -f ', tmpFile]);
	end
	
	chunkSplit = 1000;
	for j = 1:chunkSplit:size(VG,1)
		jEnd = min(size(VG,1), j+chunkSplit-1);
		VG(j:jEnd,:) = VG(j:jEnd,:)*V;
	end
end
sg = s;

ticSum = ticSum + toc;
disp(['Completed everything']);
disp(['Time:  ', num2str(ticSum), ' seconds']);
