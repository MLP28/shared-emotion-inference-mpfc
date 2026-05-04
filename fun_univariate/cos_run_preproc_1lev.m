function cos_run_preproc_1lev(S)
% COS_RUN_PREPROC_1LEV  Preprocessing + First-level GLM estimation (SPM12)
%
% This function runs the full fMRI pipeline for a set of subjects:
%   (1) Preprocessing 
%   (2) First-level GLM specification and estimation (task-based, with 
%       parametric modulators)
%
% The pipeline is implemented using SPM12 and loops over selected subjects.
%
% -------------------------------------------------------------------------
% USAGE
% -------------------------------------------------------------------------
%   cos_run_preproc_1lev(S)
%
% INPUT
%   S : vector of subject indices (e.g., S = 1:10)
%       Indices refer to subjects detected automatically in the data folder
%       (see section "DATA ORGANIZATION").
%
% EXAMPLE
%   cos_run_preproc_1lev(1:5)
%
% -------------------------------------------------------------------------
% PIPELINE OVERVIEW
% -------------------------------------------------------------------------
% For each subject:
%   - Load subject-specific behavioral/onset file
%   - Run preprocessing (optional, controlled by flags)
%   - Run first-level GLM estimation
%
% Preprocessing is skipped if already completed unless forced.
%
% -------------------------------------------------------------------------
% DATA ORGANIZATION (REQUIRED)
% -------------------------------------------------------------------------
% Root directory (D.path.rootPath) must contain:
%
%   data/
%       COSSA/
%           COSXX/
%               (raw MRI data)
%       COST/
%           COSXX/
%               (raw MRI data)
%
%   data/behav/
%       trialinfo_MRIstudy/
%           <MRI_sub>_trialinfo.mat
%       emointensity_prestudy/
%           <prestudy1_sub>_intensite_emo_ctxt.mat
%       infRT_prestudy/
%           item_sujet_RT_ctxt_pretests2.xlsx
%       infotask/
%           infoemo.xlsx
%
%   results_univariate/
%       (created automatically)
%
% Each subject must have:
%   - Functional MRI runs
%   - A corresponding trialinfo .mat file with onsets
%
% -------------------------------------------------------------------------
% OUTPUTS
% -------------------------------------------------------------------------
% For each subject:
%   results_univariate/<group>/<subject>/
%       run_*/                  % preprocessed functional data
%       anat/                   % preprocessed anatomical data
%       glm_<task>/             % first-level SPM model (SPM.mat, beta maps, etc.)
%
% Error logs (if any):
%   fun_univariate/error/<SUBJECT>_error.txt
%
% -------------------------------------------------------------------------
% DEPENDENCIES
% -------------------------------------------------------------------------
% Required toolboxes / functions:
%   - SPM12
%   - cos_preprocess.m (in univariate directory)
%   - cos_estimate_SPM_video.m (in univariate directory)
%   - get_files.m (in fun directory)
%
% IMPORTANT:
%   - SPM orthogonalization has been manually disabled in:
%       spm_get_ons.m (line 116)
%       spm_fMRI_design.m (lines 242244)
%   See:
%   http://imaging.mrc-cbu.cam.ac.uk/imaging/ParametricModulations
%
% -------------------------------------------------------------------------
% NOTES FOR USERS
% -------------------------------------------------------------------------
% - Ensure all paths are correctly set before running.
% - Check that behavioral files match subject IDs.
% - This script assumes a fixed naming convention ("COS*").
% - Parallelization is not implemented; subjects are processed sequentially.
%
% -------------------------------------------------------------------------
% AUTHOR
% -------------------------------------------------------------------------
% Pierre Gagnepain & Marine Le Petit, U1077 NIMH Caen
% 2026
%
% -------------------------------------------------------------------------
%%
% =========================================================================
% USER SETTINGS 
% =========================================================================

% Define version of spm
D.spm_path = '/qcyceron/NIMH/vol1/servier1/MULTIBRAIN/fmri_tool/spm12_last/';

% Note: COMMENTED SPM orthogonalisation in spm_get_ons (line 116) and in
% spm_fMRI_design lines 242-244
% see http://imaging.mrc-cbu.cam.ac.uk/imaging/ParametricModulations

% % set path and directory
D.path.rootPath      = '/qcyceron/NIMH/vol7/COSIMAGE/cosimage_ml/activation/matlab_codes_for_publication';

% participant codes                
D.groupname             = {'COSSA'  'COST'};

% MRI parameters
D.tr          = 2.382;
D.nslice       = 42;
D.sliceorder    = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42]; % interleaved ascending

% control what to run
D.torun.prepro      = 1; % 1 = run preprocessing; 0 = do not run preprocessing
forcepreprocess     = 0; % 1 = do preprocessing even if it has been done
D.realign           = 2; % 1 = realign & unwarp (if fieldmap); 2 = realign & reslice (classical)
D.torun.copyraw     = 1; % 1 = copyraw data to results folder, 0 = do not copy rawdata

% task info
D.task        = {'emo'};
D.n_sess_task = [4];
D.nSessions   = sum(D.n_sess_task);

% scans number for each run
D.nscans{1}            	= [275 275 275 275]; % 1-4 = task runs

% =========================================================================
% OTHER SETTINGS 
% =========================================================================

% find TPMs and add SPM path
D.tpm_path = fullfile(D.spm_path, 'tpm', 'TPM.nii');
addpath(D.spm_path);

D.path.functionRoot = fullfile(D.path.rootPath,'fun');
D.path.univfunction = fullfile(D.path.rootPath,'fun_univariate');
% add function path
addpath(genpath(D.path.functionRoot)); % Catch sight of the functions
addpath(genpath(D.path.univfunction)); 

D.path.behav    = fullfile(D.path.rootPath,'data','behav'); % contains folders for onset info, probability rating from prestudy 1 and inference timing from prestudy 2

error_dir = fullfile(D.path.univfunction,'error');
mkdir(error_dir);

SubCODE     = {};
groupidx    = [];
subidix     = [];

% workout participant
% Build subject list across groups and assign global indexing (subI)
% groupidx maps each subject to its group (COSSA / COST)
% subidix stores subject position within group
for g = 1:length(D.groupname)
    
    pth         = fullfile(D.path.rootPath,'data',D.groupname{g});
    subpth      = get_files(pth,'COS*');

    groupidx    = [groupidx;repmat(g,size(subpth,1),1)];
    for s = 1:size(subpth,1)
        subidix = [subidix,s];
        [pa,na,~] = fileparts(strtrim(subpth(s,:)));
        SubCODE{g}{s,1} = na;
    end
    
end

Nsub = length(groupidx);

D.filt          = {'^swafcos.*\.nii$'};
  
%% Subject and task loops

% Loop over selected subjects (indices defined externally)
% Each iteration runs preprocessing + first-level GLM
for subI = S;

    % subject path
    D.subcode     = SubCODE{groupidx(subI)}{subidix(subI)};
    gpcode      = D.groupname{groupidx(subI)};
    subCount    = D.subcode(end-1:end);
    % Define subject-specific input (raw) and output (results) directories
    D.sub_dir   = fullfile(D.path.rootPath,'results_univariate',gpcode,D.subcode);
    D.raw_dir   = fullfile(D.path.rootPath,'data',gpcode,D.subcode);
    D.mridir    = fullfile(D.sub_dir,'anat');
    
    D.subdata   = fullfile(D.path.behav,'trialinfo_MRIstudy',sprintf('%s%s_trialinfo.mat',D.groupname{groupidx(subI)},subCount)); % individual file trial info with onsets
    
    try

        % =========================================================================
        % PREPROCESSING
        % =========================================================================
        
        % Run preprocessing only if enabled in configuration
        % Can be skipped if preprocessing already performed
        if D.torun.prepro

            % Check whether preprocessing outputs already exist for run 1
            % Used to avoid redundant preprocessing
            findimg = get_files(fullfile(D.sub_dir,'run_1'),'swafcos*.nii');

            % preprocessing pipeline
            if isempty(findimg) || forcepreprocess == 1
                cos_preprocess(D) %
            end

        end

        
        % =========================================================================
        % FIRST LEVEL
        % =========================================================================

        % glm path
        D.glmpath = fullfile(D.sub_dir,sprintf('glm_%s',D.task{1}));
        mkdir(D.glmpath);
        
        cos_estimate_SPM_video(D)
       
    catch
        
        % Capture subject-level error and write log file for debugging
        disp(['Error with ',D.subcode])
        err = lasterror;
        fail = fopen(fullfile(error_dir,[D.subcode,'_error.txt']), 'w');
        fprintf(fail,'%s',err.message);
        fclose(fail);
        
        return;
    end

end

