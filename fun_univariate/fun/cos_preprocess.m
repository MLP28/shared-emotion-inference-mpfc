function cos_preprocess(D)
% COS_PREPROCESS  fMRI preprocessing pipeline using SPM12 (single subject)
%
% This function performs the full preprocessing pipeline for one subject.
% % It is designed to be called internally by cos_run_preproc_1lev.m and
% should not be used as a standalone entry point unless the structure `D`
% is correctly specified.
%
% -------------------------------------------------------------------------
% USAGE
% -------------------------------------------------------------------------
%   cos_preprocess(D)
%
% INPUT
%   D : structure containing all subject-specific parameters and paths.
%       This structure is created and populated in cos_preproc_1lev_batch.
%
%       Required fields include:
%           - D.raw_dir        : path to raw data (func/ and anat/)
%           - D.sub_dir        : output directory (results)
%           - D.nSessions      : total number of fMRI runs
%           - D.n_sess_task    : number of runs per task
%           - D.tr             : repetition time (in seconds)
%           - D.nslice         : number of slices
%           - D.sliceorder     : slice acquisition order
%           - D.realign        : realignment option (1 = realign+unwarp, 2 = realign only)
%           - D.torun.copyraw  : whether to copy raw data to results directory
%
% -------------------------------------------------------------------------
% PIPELINE OVERVIEW
% -------------------------------------------------------------------------
% The following preprocessing steps are applied using SPM12:
%
%   0. Data preparation
%       - Create session directories
%       - Copy functional and anatomical data (optional)
%
%   1. Slice timing correction
%   2. Motion correction
%       - Realign & Unwarp (if D.realign == 1)
%       - Realign only (if D.realign == 2)
%   3. Anatomical segmentation (tissue probability maps)
%   4. Skull stripping (via imcalc)
%   5. Coregistration (EPI * anatomical space)
%   6. Normalization to MNI space
%       - Functional images (3 mm isotropic)
%       - Skull-stripped anatomical (1 mm isotropic)
%   7. Spatial smoothing (Gaussian kernel, 10 mm FWHM)
%
% All steps are executed via SPM batch (spm_jobman).
%
% -------------------------------------------------------------------------
% INPUT DATA REQUIREMENTS
% -------------------------------------------------------------------------
% Raw data must follow this structure:
%
%   D.raw_dir/
%       func/
%           *run_1*.nii
%           *run_2*.nii
%           ...
%       anat/
%           *anat_t1.nii
%
% Functional files must match the pattern:
%   '^fcos.*\.nii$' (after copying)
%
% -------------------------------------------------------------------------
% OUTPUTS
% -------------------------------------------------------------------------
% Preprocessed data are written in:
%
%   D.sub_dir/
%       run_1/, run_2/, ...     % preprocessed functional runs
%       anat/                  % anatomical + segmentation outputs
%
% Output prefixes follow SPM conventions:
%   a  : slice timing corrected
%   r  : realigned
%   u  : unwarped (if applicable)
%   w  : normalized
%   s  : smoothed
%
% Example:
%   swafcos*.nii = fully preprocessed functional images
%
% -------------------------------------------------------------------------
% DEPENDENCIES
% -------------------------------------------------------------------------
%   - SPM12 (spm_jobman, spm_select, etc.)
%   - get_files.m (custom file selection utility)
%
% -------------------------------------------------------------------------
% IMPORTANT NOTES
% -------------------------------------------------------------------------
% - This function assumes that file naming conventions are strictly respected.
% - No quality control (QC) metrics are computed automatically.
% - No parallelization is implemented.
% - Existing files may be overwritten depending on upstream settings.
%
% -------------------------------------------------------------------------
% create session directory and copy raw data
count = 0;
sessdir = {};

% Create session-wise directories and optionally copy raw functional data
for sessI = 1:D.nSessions
    count = count + 1;
        
    files   = {};
    % Select raw functional images for each run using naming convention
    files   = get_files(fullfile(D.raw_dir,'func'),[sprintf('*run_%d*.nii',sessI)]);
    sessdir{count} = fullfile(D.sub_dir,[sprintf('run_%d',sessI)]);
    
    if D.torun.copyraw
        if exist(sessdir{count})==0
            mkdir(sessdir{count})
        end
        copyfile(files,sessdir{count})
    end
end

% Select structural T1 image (assumes single file per subject)
files   = {};
files   = get_files(fullfile(D.raw_dir,'anat'),'*anat_t1.nii');
mridir  = fullfile(D.sub_dir,'anat');
if exist(mridir)==0
   if D.torun.copyraw
   mkdir(mridir)
   copyfile(files,mridir)
   end
end

anat_file = get_files(mridir,'acos*anat_t1.nii');

%------------------------------------------------------------------------
%% Run preprocessing 
% -------------------------------------------------------------------------
spm_jobman('initcfg');

matlabbatch = {};
% =========================================================================
% Step  1 - Slice Timing 
% =========================================================================

% Build SPM slice timing input structure (one cell per session)
% Each session contains all functional volumes
for nS=1:sum(D.n_sess_task)
        files   = [];
        [files,dirs] = spm_select('ExtFPListRec',sessdir{nS},'^fcos.*\.nii$',Inf);
        FI = {};
        for nF = 1:size(files,1)
            FI{nF} = strtrim(files(nF,:));
        end
        matlabbatch{1}.spm.temporal.st.scans{nS} = FI';
end

matlabbatch{1}.spm.temporal.st.nslices  = D.nslice;
matlabbatch{1}.spm.temporal.st.tr       = D.tr;
matlabbatch{1}.spm.temporal.st.ta       = D.tr-(D.tr/D.nslice);
matlabbatch{1}.spm.temporal.st.so       = D.sliceorder;
matlabbatch{1}.spm.temporal.st.refslice = D.nslice - 1; % which slice to align the timings to
matlabbatch{1}.spm.temporal.st.prefix   = 'a';


% =========================================================================
% Step  2 - Realign (and Unwarp) EPI data
% =========================================================================

% Choice of motion correction strategy:
% 1 = Realign + Unwarp (accounts for susceptibility-by-movement interaction)
% 2 = Classical realignment only
if D.realign == 1
    for nS = 1:sum(D.n_sess_task)
        matlabbatch{2}.spm.spatial.realignunwarp.data{nS}(1) = cfg_dep(sprintf('Slice Timing: Slice Timing Corr. Images (Sess %d)',nS), substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{nS}, '.','files'));
    end
    
    matlabbatch{2}.spm.spatial.realignunwarp.eoptions.quality   = 0.9;
    matlabbatch{2}.spm.spatial.realignunwarp.eoptions.sep       = 4;
    matlabbatch{2}.spm.spatial.realignunwarp.eoptions.fwhm      = 5;
    matlabbatch{2}.spm.spatial.realignunwarp.eoptions.rtm       = 0;
    matlabbatch{2}.spm.spatial.realignunwarp.eoptions.einterp   = 2;
    matlabbatch{2}.spm.spatial.realignunwarp.eoptions.ewrap     = [0 0 0];
    matlabbatch{2}.spm.spatial.realignunwarp.eoptions.weight    = '';
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.basfcn  = [12 12];
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.regorder= 1;
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.lambda  = 100000;
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.jm      = 0;
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.fot     = [4 5];
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.sot     = [];
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.uwfwhm  = 4;
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.rem     = 1;
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.noi     = 5;
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.expround= 'Average';
    matlabbatch{2}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
    matlabbatch{2}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
    matlabbatch{2}.spm.spatial.realignunwarp.uwroptions.wrap    = [0 0 0];
    matlabbatch{2}.spm.spatial.realignunwarp.uwroptions.mask    = 1;
    matlabbatch{2}.spm.spatial.realignunwarp.uwroptions.prefix  = 'u';

elseif D.realign == 2
    
    for nS = 1:sum(D.n_sess_task)
        matlabbatch{2}.spm.spatial.realign.estimate.data{nS}(1) = cfg_dep(sprintf('Slice Timing: Slice Timing Corr. Images (Sess %d)',nS), substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{nS}, '.','files'));
    end
    
    matlabbatch{2}.spm.spatial.realign.estimate.eoptions.quality = 0.9;
    matlabbatch{2}.spm.spatial.realign.estimate.eoptions.sep = 4;
    matlabbatch{2}.spm.spatial.realign.estimate.eoptions.fwhm = 5;
    matlabbatch{2}.spm.spatial.realign.estimate.eoptions.rtm = 1;
    matlabbatch{2}.spm.spatial.realign.estimate.eoptions.interp = 2;
    matlabbatch{2}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
    matlabbatch{2}.spm.spatial.realign.estimate.eoptions.weight = '';
    
end


% =========================================================================
% Step  3 - Segment anat
% =========================================================================
matlabbatch{3}.spm.spatial.preproc.channel.vols             = {anat_file};
matlabbatch{3}.spm.spatial.preproc.channel.biasreg          = 0.0001;
matlabbatch{3}.spm.spatial.preproc.channel.biasfwhm         = 60;
matlabbatch{3}.spm.spatial.preproc.channel.write            = [0 1];
ngaus  = [2 2 2 3 4 2];
native = [1 1 1 0 0 0];

% Tissue classes follow SPM TPM order:
% GM, WM, CSF, bone, soft tissue, air/background
for i = 1:6
    matlabbatch{3}.spm.spatial.preproc.tissue(i).tpm    = {sprintf('%s,%d',D.tpm_path,i)};
    matlabbatch{3}.spm.spatial.preproc.tissue(i).ngaus  = ngaus(i);
    matlabbatch{3}.spm.spatial.preproc.tissue(i).native = [native(i) native(i)];
    % execute a job
    if i == 1
        matlabbatch{3}.spm.spatial.preproc.tissue(i).warped = [1 0];
    else
        matlabbatch{3}.spm.spatial.preproc.tissue(i).warped = [0 0];
    end
end

matlabbatch{3}.spm.spatial.preproc.warp.mrf             = 1;
matlabbatch{3}.spm.spatial.preproc.warp.cleanup         = 1;
matlabbatch{3}.spm.spatial.preproc.warp.reg             = [0 0.001 0.5 0.05 0.2];
matlabbatch{3}.spm.spatial.preproc.warp.affreg          = 'mni';
matlabbatch{3}.spm.spatial.preproc.warp.fwhm            = 0;
matlabbatch{3}.spm.spatial.preproc.warp.samp            = 3;
matlabbatch{3}.spm.spatial.preproc.warp.write           = [1 1];

% get path
% ---------------------------------
matlabbatch{4}.cfg_basicio.file_dir.cfg_fileparts.files(1) = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));

% =========================================================================
% Step  5 - Skull strip anat 
% =========================================================================
matlabbatch{5}.spm.util.imcalc.input(1)         = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
matlabbatch{5}.spm.util.imcalc.input(2)         = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
matlabbatch{5}.spm.util.imcalc.input(3)         = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
matlabbatch{5}.spm.util.imcalc.input(4)         = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
matlabbatch{5}.spm.util.imcalc.output           = 'Brain';
matlabbatch{5}.spm.util.imcalc.outdir(1)        = cfg_dep('Get Pathnames: Directories (unique)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','up'));
matlabbatch{5}.spm.util.imcalc.expression       = '(i1 + i2 + i3) .* i4'; % Skull stripping: combine GM+WM+CSF probability maps multiplied by bias-corrected T1 image
matlabbatch{5}.spm.util.imcalc.var              = struct('name', {}, 'value', {});
matlabbatch{5}.spm.util.imcalc.options.dmtx     = 0;
matlabbatch{5}.spm.util.imcalc.options.mask     = 0;
matlabbatch{5}.spm.util.imcalc.options.interp   = 1;
matlabbatch{5}.spm.util.imcalc.options.dtype    = 4;

% =========================================================================
% Step  6 - Coregister the structural to EPI run1 and apply to all
% other EPI
% =========================================================================

matlabbatch{6}.spm.spatial.coreg.estimate.ref(1)            = cfg_dep('Image Calculator: ImCalc Computed Image: Brain', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));%{brmri};%cfg_dep('Image Calculator: Imcalc Computed Image', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{6}.spm.spatial.coreg.estimate.source(1)         = cfg_dep('Realign: Estimate: Realigned Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','cfiles'));

for nS = 1:sum(D.n_sess_task)
    matlabbatch{6}.spm.spatial.coreg.estimate.other(nS) = cfg_dep(sprintf('Realign: Estimate: Realigned Images (Sess %d)',nS), substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{nS}, '.','cfiles'));
end

matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2];
matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];


% =========================================================================
% Step  6 - Write images using normalisation parameter
% =========================================================================

% Normalise Functionnal image
matlabbatch{7}.spm.spatial.normalise.write.subj.def(1)          = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));

for nS = 1:sum(D.n_sess_task)
    matlabbatch{7}.spm.spatial.normalise.write.subj.resample(nS) =  cfg_dep(sprintf('Realign: Estimate: Realigned Images (Sess %d)',nS), substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{nS}, '.','cfiles'));
end

matlabbatch{7}.spm.spatial.normalise.write.woptions.bb      = [-78 -112 -70
    78 76 85];
matlabbatch{7}.spm.spatial.normalise.write.woptions.vox     = [3 3 3];
matlabbatch{7}.spm.spatial.normalise.write.woptions.interp  = 4;

% Normalise Brain image
matlabbatch{8}.spm.spatial.normalise.write.subj.def(1)        = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{8}.spm.spatial.normalise.write.subj.resample(1)   = cfg_dep('Image Calculator: Imcalc Computed Image', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{8}.spm.spatial.normalise.write.woptions.vox       = [1 1 1];


% =========================================================================
% Step  7 - Smooth to the EPI images
% =========================================================================
matlabbatch{9}.spm.spatial.smooth.data(1)   = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{9}.spm.spatial.smooth.fwhm      = repmat(10,1,3); % Gaussian smoothing kernel (10mm FWHM)
matlabbatch{9}.spm.spatial.smooth.dtype     = 0;
matlabbatch{9}.spm.spatial.smooth.im        = 0;
matlabbatch{9}.spm.spatial.smooth.prefix    = 's';


spm_jobman('run', matlabbatch)
% spm_jobman('interactive', matlabbatch)