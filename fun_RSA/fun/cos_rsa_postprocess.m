function cos_rsa_postprocess(D)
% COS_RSA_POSTPROCESS  Normalise and smooth individual RSA searchlight maps (SPM12)
%
% After the searchlight RSA, individual correlation maps are in native space.
% This function applies two post-processing steps:
%   1. Spatial normalisation to MNI space (using subject's deformation field)
%   2. Spatial smoothing (Gaussian kernel)
%
% The deformation field (y_*.nii) is expected in the subject's anatomical
% directory, as produced by SPM12 segmentation during univariate preprocessing.
%
% Called by cos_rsa_run_firstlevel.m - do not call directly.
%
% INPUT
%   D : configuration structure (set in cos_rsa_run_firstlevel.m)
%       Required fields:
%           D.raw_dir               - subject univariate directory (contains anat/)
%           D.sub_dir               - subject RSA directory
%           D.glm_model             - GLM folder name (used to locate RSA maps)
%           D.searchlight.size      - searchlight radius (mm), used in filenames
%           D.searchlight.model_name     - label used in RSA output filenames
%           D.searchlight.rsm_of_interest  - cell array of RSM names
%           D.smooth_fwhm           - smoothing kernel [x y z] in mm
%           
% -------------------------------------------------------------------------

%% -------------------------------------------------------------------------
% LOCATE INPUT FILES
% -------------------------------------------------------------------------

spher    = D.searchlight.size;
mod      = D.searchlight.model_name;
cond     = D.searchlight.rsm_of_interest;

% Searchlight input directory
searchlight_dir  = fullfile(D.sub_dir, 'searchlight', sprintf('rsa_%s_%dmm',  D.searchlight.model_name, D.searchlight.size),D.glm_model);

% Deformation field from SPM12 segmentation (y_*.nii in anat folder)
norm_file = get_files(fullfile(D.raw_dir, 'anat'), 'y_*.nii');
if isempty(norm_file)
    error('Deformation field not found in %s', fullfile(D.sub_dir, 'anat'));
end

% RSA map files to normalise (one per RSM of interest)
resample_files = {};
for nS = 1:length(cond)
    fname = fullfile(searchlight_dir, sprintf('sphere%d-cov-%s_%s.nii', spher, mod, cond{nS}));
    if ~exist(fname, 'file')
        error('RSA map not found: %s', fname);
    end
    resample_files{nS} = fname;
end


%% -------------------------------------------------------------------------
% SPM BATCH: NORMALISE + SMOOTH
% -------------------------------------------------------------------------

spm_jobman('initcfg');
matlabbatch = {};

% --- Step 1: Normalise to MNI space ---------------------------------------
matlabbatch{1}.spm.spatial.normalise.write.subj.def = {norm_file};

for nS = 1:length(resample_files)
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample(nS,1) = {resample_files{nS}};
end

matlabbatch{1}.spm.spatial.normalise.write.woptions.bb     = [-78 -112 -70; 78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox    = [1 1 1];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

% --- Step 2: Smooth normalised maps ---------------------------------------
% Input: normalised images from step 1 (cfg_dep links the two steps)
matlabbatch{2}.spm.spatial.smooth.data = cfg_dep('Normalise: Write: Normalised Images (Subj 1)',substruct('.','val','{}',{1},'.','val','{}',{1},'.','val','{}',{1},'.','val','{}',{1}), substruct('()',{1},'.','files'));
matlabbatch{2}.spm.spatial.smooth.fwhm   = D.smooth_fwhm;
matlabbatch{2}.spm.spatial.smooth.dtype  = 0;
matlabbatch{2}.spm.spatial.smooth.im     = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = 's';


%% -------------------------------------------------------------------------
% RUN
% -------------------------------------------------------------------------

fprintf('  Normalising and smoothing %d RSA map(s) (FWHM = %d %d %d mm)...\n', ...
        length(cond), D.smooth_fwhm(1), D.smooth_fwhm(2), D.smooth_fwhm(3));

spm_jobman('run', matlabbatch);

end