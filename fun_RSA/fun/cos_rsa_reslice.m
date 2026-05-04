function cos_rsa_reslice(D)
% COS_RSA_RESLICE  Copy and reslice functional images for RSA analysis (SPM12)
%
% Functional images preprocessed during the univariate analysis (afcos*.nii,
% in native space) are copied from the univariate results directory to the
% RSA results directory, then resliced to ensure voxel-wise alignment.
% This step is required before first-level GLM estimation for RSA.
%
% Processing per session:
%   1. Copy afcos*.nii from results_univariate/<group>/<sub>/run_X
%      to    results_RSA/<group>/<sub>/preprocessing/run_X
%   2. Reslice all volumes using the first volume as reference
%      -> produces rafcos*.nii in the preprocessing folder
%
% Called by cos_rsa_run_firstlevel.m - do not call directly.
%
% INPUT
%   D : configuration structure (set in cos_rsa_run_firstlevel.m)
%       Required fields:
%           D.raw_dir  - univariate results directory for this subject
%                        (source of preprocessed afcos*.nii)
%           D.sub_dir  - RSA results directory for this subject
%                        (destination: sub_dir/preprocessing/run_X/)
%           D.filt     - regexp to select copied afcos*.nii in destination
%                        folder for spm_reslice (e.g. '^afcos.*\.nii$')
%           D.nscans   - number of scans per run (cell array)
%
% OUTPUT
%   Resliced images (rafcos*.nii) written to:
%       D.sub_dir/preprocessing/run_X/
%
% -------------------------------------------------------------------------

allsess_files   = [];

sessdir_out = fullfile(D.sub_dir, 'preprocessing');
mkdir(sessdir_out);
    
for ses = 1:length(D.nscans{1})

    % --- Define source and destination directories ------------------------
    sessdir_in  = fullfile(D.raw_dir, sprintf('run_%d', ses));

    % --- Copy .nii files to RSA preprocessing folder ----------------------
    nii_files = dir(fullfile(sessdir_in, 'afcos*.nii'));
    if isempty(nii_files)
        warning('No afcos*.nii found in %s  skipping session %d.', sessdir_in, ses);
    else
        copyfile(fullfile(sessdir_in, nii_files.name), sessdir_out);
    end  

    

end

% --- Select copied images in destination folder -----------------------
func_files = spm_select('ExtFPList', sessdir_out, D.filt{1});

    
% --- Reslice ----------------------------------------------------------
% Use default values
% Output: rafcos*.nii written in sessdir_out (prefix 'r' added by SPM).
spm_reslice(func_files);

end

