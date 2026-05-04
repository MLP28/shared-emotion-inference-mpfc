function SPM  = get_spmdesign(SPM,Sess)


% adapted from  spm_fMRI_design and spm_fmri_concatenate !!!!!


SVNid = '$Rev: 7018 $';


%-Session and scan number
%--------------------------------------------------------------------------
nscan = SPM.nscan;
nsess = length(nscan);


%==========================================================================
% - D E S I G N   P A R A M E T E R S
%==========================================================================
SPM.SPMid = spm('FnBanner',mfilename,SVNid);

%-High-pass filtering
%==========================================================================

%-Low frequency confounds
%--------------------------------------------------------------------------
try
    HParam     = [SPM.xX.K(:).HParam];
    if length(HParam) == 1 
        HParam = repmat(HParam,1,nsess);
    end
catch
    error('High-pass filter not specified.');
end

% High-pass filter (see spm_spm.m)
%--------------------------------------------------------------------------
s = cumsum([0 SPM.nscan]);

for i=1:numel(SPM.nscan)
    K(i) = struct('HParam', SPM.xX.K(1).HParam,...
                  'row',    s(i) + (1:SPM.nscan(i)),...
                  'RT',     SPM.xY.RT);
end
SPM.xX.K = spm_filter(K);

% Temporal non-sphericity (see spm_spm.m)
%--------------------------------------------------------------------------
switch lower(SPM.xVi.form)
    case {'ar(1)','ar(0.2)'}
        SPM.xVi.Vi   = spm_Ce('ar',SPM.nscan,0.2);
        SPM.xVi.form = 'AR(0.2)';
    case 'fast'
        SPM.xVi.Vi   = spm_Ce('fast',SPM.nscan,SPM.xY.RT);
    case {'i.i.d', 'none'}
    otherwise
        warning('Unhandled temporal non-sphericity.');
end



%==========================================================================
% - C O N F I G U R E   D E S I G N
%==========================================================================

spm('Pointer','Watch');

%-Get image files
%==========================================================================

%-Map files
%--------------------------------------------------------------------------
fprintf('%-40s: ','Mapping files')                                      %-#
VY    = spm_data_hdr_read(SPM.xY.P);
fprintf('%30s\n','...done')                                             %-#

%-Check internal consistency of images
%--------------------------------------------------------------------------
spm_check_orientations(VY);

%-Place mapped files in xY
%--------------------------------------------------------------------------
SPM.xY.VY = VY;


%-Compute Global variate
%==========================================================================
GM    = 100;
q     = length(VY);
g     = zeros(q,1);
fprintf('%-40s: ','Calculating globals')                                %-#
spm_progress_bar('Init',q,'Calculating globals');
if spm_mesh_detect(VY)
    for i = 1:q
        dat = spm_data_read(VY(i));
        g(i) = mean(dat(~isnan(dat)));
        spm_progress_bar('Set',i)
    end
else
    for i = 1:q
        g(i) = spm_global(VY(i));
        spm_progress_bar('Set',i)
    end
end
spm_progress_bar('Clear');
fprintf('%30s\n','...done')                                             %-#

%-Scale if specified (otherwise session specific grand mean scaling)
%--------------------------------------------------------------------------
gSF   = GM./g;
if strcmpi(SPM.xGX.iGXcalc,'none')
    for i = 1:nsess
        gSF(Sess(i).row) = GM./mean(g(Sess(i).row));
    end
end

%-Apply gSF to memory-mapped scalefactors to implement scaling
%--------------------------------------------------------------------------
for i = 1:q
    SPM.xY.VY(i).pinfo(1:2,:) = SPM.xY.VY(i).pinfo(1:2,:) * gSF(i);
    if spm_mesh_detect(VY)
        SPM.xY.VY(i).private.private.data{1}.data.scl_slope = ...
            SPM.xY.VY(i).private.private.data{1}.data.scl_slope * gSF(i);
        SPM.xY.VY(i).private.private.data{1}.data.scl_inter = ...
            SPM.xY.VY(i).private.private.data{1}.data.scl_inter * gSF(i);
    else
        SPM.xY.VY(i).private.dat.scl_slope = ...
            SPM.xY.VY(i).private.dat.scl_slope * gSF(i);
        SPM.xY.VY(i).private.dat.scl_inter = ...
            SPM.xY.VY(i).private.dat.scl_inter * gSF(i);
    end
end

%-Place global variates in xGX
%--------------------------------------------------------------------------
SPM.xGX.sGXcalc = 'mean voxel value';
SPM.xGX.sGMsca  = 'session specific';
SPM.xGX.rg      = g;
SPM.xGX.GM      = GM;
SPM.xGX.gSF     = gSF;


%-Masking
%==========================================================================

%-Masking threshold, as proportion of globals
%--------------------------------------------------------------------------
try
    gMT = SPM.xM.gMT;
catch
    gMT = spm_get_defaults('mask.thresh');
end
TH = g.*gSF*gMT;

%-Place masking structure in xM
%--------------------------------------------------------------------------
SPM.xM = struct(...
    'T',   ones(q,1),...
    'TH',  TH,...
    'gMT', gMT,...
    'I',   0,...
    'VM',  {[]},...
    'xs',  struct('Masking','analysis threshold'));


%-Design description - for saving and display
%==========================================================================
xs = struct(...
    'Global_calculation',   SPM.xGX.sGXcalc,...
    'Grand_mean_scaling',   SPM.xGX.sGMsca,...
    'Global_normalisation', SPM.xGX.iGXcalc);
for fn=(fieldnames(xs))', SPM.xsDes.(fn{1}) = xs.(fn{1}); end


%==========================================================================
% - S A V E   A N D   D I S P L A Y
%==========================================================================

%-Save SPM.mat
%--------------------------------------------------------------------------
%if ~nargout
    fprintf('%-40s: ','Saving SPM configuration')                       %-#
    save('SPM.mat', 'SPM', spm_get_defaults('mat.format'));
    fprintf('%30s\n','...SPM.mat saved')                                %-#
%end


%-Display design report
%--------------------------------------------------------------------------
if ~spm('CmdLine') && ~isempty(spm_figure('FindWin','Graphics'))
    fprintf('%-40s: ','Design reporting')                               %-#
    try,   fname = reshape(cellstr(SPM.xY.P),size(SPM.xY.VY));
    catch, fname = {}; end
    spm_DesRep('DesMtx',SPM.xX,fname,SPM.xsDes)
    fprintf('%30s\n','...done')                                         %-#
end

spm('Pointer','Arrow');



