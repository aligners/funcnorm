function funcnorm_prepare_align(subjects, hems, saveSubjDir, outDir, mmRes, fs_surf)

nsubj = length(subjects);
nhem = length(hems);

for s=1:nsubj
        sub=subjects{s};
        for h=1:nhem
                hem=hems{h};

                disp(['subject ',sub,' ',hem,' hemisphere...']);
                
                orig_prefix=[saveSubjDir,'/',sub,'/SUMA/',hem,'.'];
                orig_sphere=[orig_prefix,fs_surf,'.asc'];
                orig_pial=[orig_prefix,'pial.asc'];
                orig_smoothwm=[orig_prefix,'smoothwm.asc'];
                warpDir=[outDir,'/warps/',hem];
                warp_sphere=[warpDir,'/standard',num2str(mmRes),'mm_',sub,'_',hem,'.funcnorm.asc'];

                % read in warped, low-res mesh and original, high-res mesh
                [warpCoords,warpTris] = read_fs_ascii_mesh(warp_sphere);
                origCoords = read_fs_ascii_mesh(orig_sphere);
                % Find nearest node in orig_sphere from warp_sphere
                N = size(warpCoords,2);
                if ~exist('nearest_nodes','var')
                    nearest_nodes = zeros(1, N, 'single');
                end
                disp('finding nearest nodes...');
                for n = 1:N
                    corrs = warpCoords(:, n)'*origCoords;
                    % Find node with max correlation
                    [garbage, nearest_nodes(n)] = max(corrs);
                end
          
                % use nearest nodes to get rows of X Y Z from orig_pial and orig_smoothwm
                % and use those X Y Z coordinates with the triangles info from warp_sphere
                origCoords = read_fs_ascii_mesh(orig_pial);
                warpPrefix = [warpDir,'/',sub,'_',hem];
                make_warped_nonsphere(origCoords,nearest_nodes,warpTris,[warpPrefix,'.pial.funcnorm.asc']);
                origCoords = read_fs_ascii_mesh(orig_smoothwm);
                make_warped_nonsphere(origCoords,nearest_nodes,warpTris,[warpPrefix,'.smoothwm.funcnorm.asc']);
                  
                % make a spec file with warp_sphere, warp_pial, warp_smoothm in it
                fid = fopen([warpPrefix,'.spec'],'w');
                fprintf(fid,'\tGroup = %s_%s_funcnormed\n',sub,hem);
                fprintf(fid,'\tStateDef = sphere.reg\n');
                fprintf(fid,'\tStateDef = pial\n');
                fprintf(fid,'\tStateDef = smoothwm\n');
                fprintf(fid,'\n\n');
                write_spec_surface('smoothwm');
                fprintf(fid,'\n');
                write_spec_surface('pial');
                fprintf(fid,'\n');
                write_spec_surface('sphere.reg');
                fclose(fid);
        end
end

end

function [coordinates, triangles] = read_fs_ascii_mesh(filename)

fid=fopen(filename,'r');
if -1 == fid
    error('cannot open freesurfer ascii mesh file');
end

coordsonly = (nargout == 1);

line=fgetl(fid);
% should be comment line
if line(1) ~= '#'
    error('line 1 unrecognized mesh format');
end

line=fgetl(fid);
% should be num_nodes, num_triangles
nodes_tris=sscanf(line,'%d');
if length(nodes_tris) ~= 2
    error('line 2 unrecognized mesh format');
end

N = nodes_tris(1);
coordinates = fscanf(fid,'%f',[4,N]);
if ~coordsonly
    T = nodes_tris(2);
    triangles = fscanf(fid,'%f',[4,T]);
else
    triangles = [];
end
fclose(fid);

end

function make_warped_nonsphere(origCoords, nearest_nodes, warpTris, outputfile)
    
fid = fopen(outputfile,'w');
if -1 == fid
    error('cannot open output warped mesh for writing');
end

N = size(nearest_nodes,2);
T = size(warpTris,2);
fprintf(fid,'#\n');
fprintf(fif,'%d %d\n',N,T);
for n = 1:N
    fprintf(fid,'%f %f %f 0\n',...
        origCoords(nearest_nodes(n),1),...
        origCoords(nearest_nodes(n),2),...
        origCoords(nearest_nodes(n),3));
end
for t = 1:T
    fprintf(fid,'%f %f %f 0\n',...
        warpTris(n,1),...
        warpTris(n,2),...
        warpTris(n,3));
end
fclose(fid);

end

function write_spec_surface(fid,sub,hem,state)

fprintf(fid,'NewSurface\n');
fprintf(fid,'\tSurfaceFormat = ASCII\n');
fprintf(fid,'\tSurfaceType = FreeSurfer\n');
smoothwm_file = sprintf('%s_%s.smoothwm.funcnorm.asc\n',sub,hem);
fprintf(fid,'\tFreeSurferSurface = %s_%s.%s.funcnorm.asc\n',sub,hem,state);
if strcmp(state,'smoothwm')
    localDomainParent = 'SAME';
else
    localDomainParent = smoothwm_file;
end
fprintf(fid,'\tLocalDomainParent = %s\n',localDomainParent);
fprintf(fid,'\tSurfaceState = %s\n',state);
fprintf(fid,'EmbedDimension = 3\n');

end
