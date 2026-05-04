function cos_run_ppi(S)
% COS_RUN_PPI  Entry point for PPI first-level analyses (SPM12)
%
% This function runs psychophysiological interaction (PPI) analyses
% for a set of subjects and a set of seed regions.
%
% For each subject and each seed ROI, the function:
%   1. Loads the subject's first-level SPM.mat
%   2. Extracts the seed ROI signal (coord_cont) from the seed mask
%   3. Extracts target voxels (coord_target) from the brain mask
%   4. Calls run_ppi to compute PPI maps
%
% run_ppi performs the following steps:
%   - Removes confounds (motion, session effects) from seed and target signals
%   - Deconvolves the seed signal to recover the underlying neural signal
%   - Computes the PPI interaction term (neural signal x psychological regressor)
%     for each condition (acc_dom, acc_comp, unc_dom, unc_comp, Emo, NN)
%   - Estimates a GLM at each target voxel (ppi_glm) yielding one beta map
%     per condition
%   - Computes contrast maps (ME_acc, ME_unc, CvsD_acc, DvsC_acc,
%     CvsD_unc, DvsC_unc) by linear combination of condition betas
%   - Smooths all output maps (FWHM 8mm) and saves them as .nii files
%   - Saves the full beta matrix and data structure as a .mat file
%
% -------------------------------------------------------------------------
% USAGE
% -------------------------------------------------------------------------
%   cos_run_ppi(S)
%
% INPUT
%   S : vector of subject indices (e.g., 1:10)
%       Indices refer to subjects detected automatically in the data folder.
%
% EXAMPLE
%   cos_run_ppi(1:10)
%   cos_run_ppi(3)
%
% -------------------------------------------------------------------------
% PREREQUISITES
% -------------------------------------------------------------------------
%   - First-level GLM must be estimated beforehand (cos_run_preproc_1lev.m)
%   - Seed ROI masks (.nii) must exist in D.path.seeds
%   - SPM.mat must be present in each subject's GLM directory
%   - SF.mat (stick functions) must exist in each subject's GLM directory
%     (saved by cos_estimate_SPM_video.m)

% -------------------------------------------------------------------------
% CONDITIONS AND CONTRASTS
% -------------------------------------------------------------------------
%   Defined in USER SETTINGS via structure C:
%       C.condname : cell array of condition names : must match the regressor
%                    names used in the first-level GLM (cos_estimate_SPM_video.m)
%       C.conlabel : cell array of contrast names
%       C.con      : contrast matrix [nContrasts x nConditions]
%                    rows = contrasts, columns = conditions (same order as C.condname)
%
%   Default contrasts:
%       ME_acc   = [1  1  0  0  0  0]  main effect accuracy
%       ME_unc   = [0  0  1  1  0  0]  main effect uncertainty
%       CvsD_acc = [-1 1  0  0  0  0]  competing vs dominant (accuracy)
%       DvsC_acc = [1 -1  0  0  0  0]  dominant vs competing (accuracy)
%       CvsD_unc = [0  0 -1  1  0  0]  competing vs dominant (uncertainty)
%       DvsC_unc = [0  0  1 -1  0  0]  dominant vs competing (uncertainty)
%
%
% =========================================================================
% USER SETTINGS  (edit this section only)
% =========================================================================

% --- SPM path -------------------------------------------------------------
D.spm_path = '/qcyceron/NIMH/vol1/servier1/MULTIBRAIN/fmri_tool/spm12_last/';

% --- Paths ----------------------------------------------------------------
D.path.rootPath = '/qcyceron/NIMH/vol7/COSIMAGE/cosimage_ml/activation/matlab_codes_for_publication';
D.path.seeds    = fullfile(D.path.rootPath, 'data', 'seed_masks'); % folder containing seed ROI .nii files

% --- Seed ROIs ------------------------------------------------------------
% Names must match .nii filenames in D.path.seeds (e.g. vmPFC_ROI_fig.nii)
D.seeds = {'vmPFC', 'dmPFC'};

% --- GLM model name -------------------------------------------------------
D.glm_model = 'glm_emo';   % first-level GLM folder name

% --- Subjects and groups --------------------------------------------------
D.groupname = {'COST', 'COSSA'};

% --- Conditions and contrasts ---------------------------------------------
% Must match 
C.condname = {'acc_dom_PR','acc_comp_PR','unc_dom_PR','unc_comp_PR','Emo_PR','NN_PR'};
C.conlabel    = {'ME_acc','ME_unc','CvsD_acc','DvsC_acc','CvsD_unc','DvsC_unc'};
C.con         = [1 1 0 0 0 0; 0 0 1 1 0 0; -1 1 0 0 0 0;1 -1 0 0 0 0; 0 0 -1 1 0 0; 0 0 1 -1 0 0];

%% -------------------------------------------------------------------------
% INITIALISATION
% -------------------------------------------------------------------------

addpath(D.spm_path);

D.path.PPIfunction = fullfile(D.path.rootPath, 'fun_PPI', 'fun');
addpath(genpath(D.path.PPIfunction));
D.path.functionRoot = fullfile(D.path.rootPath,'fun');
addpath(genpath(D.path.functionRoot)); 

% Error log directory
error_dir = fullfile(D.path.rootPath, 'fun_PPI', 'error');
mkdir(error_dir);

% Build subject list across groups
SubCODE  = {};
groupidx = [];
subidix  = [];

for g = 1:length(D.groupname)
    pth    = fullfile(D.path.rootPath, 'results_univariate', D.groupname{g});
    subpth = get_files(pth, 'COS*');
    groupidx = [groupidx; repmat(g, size(subpth,1), 1)];
    for s = 1:size(subpth, 1)
        subidix = [subidix, s];
        [~, na, ~] = fileparts(strtrim(subpth(s,:)));
        SubCODE{g}{s,1} = na;
    end
end


%% -------------------------------------------------------------------------
% SEED AND SUBJECT LOOPS
% -------------------------------------------------------------------------

for R = 1:length(D.seeds)

    % --- Load seed ROI mask and extract MNI coordinates ------------------
    seed_name = D.seeds{R};
    seed_file = fullfile(D.path.seeds, sprintf('%s.nii',seed_name));

    if ~exist(seed_file, 'file')
        warning('Seed mask not found, skipping: %s', seed_file);
    end

    maskvol    = spm_vol(seed_file);
    [M, xyz]   = spm_read_vols(maskvol);
    coord_cont = xyz(:, M > 0);   % MNI coordinates of seed voxels [3 x N]

    fprintf('\n=== Seed: %s (%d voxels) ===\n', seed_name, size(coord_cont,2));

    % --- Subject loop -----------------------------------------------------
    for subI = S

        D.subcode = SubCODE{groupidx(subI)}{subidix(subI)};
        gpcode    = D.groupname{groupidx(subI)};

        D.sub_dir  = fullfile(D.path.rootPath, 'results_univariate', gpcode, D.subcode);
        D.glmpath  = fullfile(D.sub_dir, D.glm_model);
        storepath  = fullfile(D.path.rootPath, 'results_PPI', D.subcode, seed_name);

        fprintf('\n[%s]  Subject: %s  |  Seed: %s\n', datestr(now,'HH:MM:SS'), D.subcode, seed_name);

        try

            % Load first-level SPM.mat
            SPM_file = fullfile(D.glmpath, 'SPM.mat');
            if ~exist(SPM_file, 'file')
                warning('SPM.mat not found for %s, skipping.', D.subcode);
            end
            load(SPM_file);
            cd(SPM.swd);

            % Extract target voxels from brain mask
            [M, xyz]      = spm_read_vols(SPM.VM);
            coord_target  = xyz(:, M > 0);   % MNI coordinates of all brain voxels

            % Output directory
            mkdir(storepath);

            % Run PPI
            analysisname = 'ppi';
            run_ppi(SPM, coord_target, coord_cont, storepath, analysisname,C);

        catch
            fprintf('  !! Error with %s | seed %s\n', D.subcode, seed_name);
            err = lasterror;
            fid = fopen(fullfile(error_dir, sprintf('%s_%s_error.txt', D.subcode, seed_name)), 'w');
            fprintf(fid, '%s', err.message);
            fclose(fid);
        end

    end % subject loop

end % seed loop

end