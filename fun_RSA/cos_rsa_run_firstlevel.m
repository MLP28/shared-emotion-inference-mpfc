function cos_rsa_run_firstlevel(S, mode)
% COS_RSA_RUN_FIRSTLEVEL  Entry point for RSA first-level analyses (SPM12 + Viviani toolbox)
%
% This function runs the RSA first-level pipeline for a set of subjects.
% It acts as a central controller for:
%   - First-level GLM estimation (item-wise beta maps for RSA)
%   - RDM computation and RSM construction
%   - Searchlight RSA (via Viviani toolbox)
%   - Post-processing: normalisation and smoothing of individual maps
%
% -------------------------------------------------------------------------
% USAGE
% -------------------------------------------------------------------------
%   cos_rsa_run_firstlevel(S, mode)
%
% INPUTS
%   S    : vector of subject indices (e.g., 1:10)
%          Indices refer to subjects detected automatically in the data folder.
%
%   mode : string controlling which steps to run:
%       'rdm'         - RDM/RSM computation only (requires SPM already estimated)
%       'reslice'     - reslice afunc
%       'spm'         - first-level GLM estimation only
%       'searchlight' - searchlight RSA only (requires RDMs already computed)
%       'postprocess' - normalisation + smoothing only (requires searchlight done)
%       'all'         - reslice -> spm -> searchlight -> postprocess (RDMs must be computed)
%
% EXAMPLE
%   cos_rsa_run_firstlevel(1:10, 'all')
%   cos_rsa_run_firstlevel(3, 'searchlight')
%
% -------------------------------------------------------------------------
% PREREQUISITES
% -------------------------------------------------------------------------
%   - SPM12 installed and path set below
%   - Viviani RSA toolbox installed in dependancies
%   - Behavioral/onset files available in D.path.behav
%   - For mode 'all': RDMs must be computed beforehand (mode 'rdm')
%   - For modes 'spm', 'searchlight', 'postprocess': previous steps must
%     be completed in order (reslice -> spm -> searchlight -> postprocess)
%
%% =========================================================================
% USER SETTINGS  (edit this section only)
% =========================================================================

% --- Toolbox paths --------------------------------------------------------
D.spm_path = '/qcyceron/NIMH/vol1/servier1/MULTIBRAIN/fmri_tool/spm12_last/';

% --- Root and output paths ------------------------------------------------
D.path.rootPath      = '/qcyceron/NIMH/vol7/COSIMAGE/cosimage_ml/activation/matlab_codes_for_publication';

% --- Task and GLM model ---------------------------------------------------
D.task      = {'emo'};
D.glm_model = 'glm_rsanative';   % name of the first-level GLM folder

% --- MRI acquisition parameters -------------------------------------------
D.tr          = 2.382;
D.nslice      = 42;
D.sliceorder  = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 ...
                 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42];

% --- Number of scans per run (one entry per run) --------------------------
D.nscans{1} = [275 275 275 275];

% --- Subjects and groups --------------------------------------------------
D.groupname = {'COST', 'COSSA'};

% --- Searchlight parameters -----------------------------------------------
D.searchlight.size   = 8;                  % searchlight radius in mm
D.searchlight.method = 'Pearson';          % 'Pearson', 'Spearman', or 'regression'
D.searchlight.brainmap_type = 'cov';       % 'sscp', 'cov', or 'cor'
D.searchlight.model_name    = 'interference_mod'; % label used in output filenames
D.searchlight.rsm_of_interest = {'compromise_sumcomp'}; % RSM(s) to test
D.searchlight.rsm_partial     = {'BCov', 'SCov', 'compromise_norm_RT', 'compromise_domint'}; % partial correlates

% --- Post-processing parameters -------------------------------------------
D.smooth_fwhm = [8 8 8];   % smoothing kernel in mm

%% -------------------------------------------------------------------------
% INITIALISATION
% -------------------------------------------------------------------------

% Add toolbox path
D.viviani_path = fullfile(D.path.rootPath, 'dependency','rsa-rsm-master');
addpath(D.spm_path);
addpath(genpath(D.viviani_path));

% Derive internal paths from root
D.path.functionRoot = fullfile(D.path.rootPath, 'fun');
D.path.behav        = fullfile(D.path.rootPath, 'data', 'behav');
D.path.rsm          = fullfile(D.path.rootPath, 'results_RSA', 'rsm');
D.filt              = {'^afcos.*\.nii$'}; % use non normalized and unsmoothed functional images for reslice
D.filt_rsa          = {'^rafcos.*\.nii$'}; % use resliced functional images for first level

addpath(genpath(D.path.functionRoot));
addpath(genpath(fullfile(D.path.rootPath, 'fun_RSA')));

% Error log directory
error_dir = fullfile(D.path.rootPath,'fun_RSA', 'error');
mkdir(error_dir);

% Build subject list across groups
SubCODE  = {};
groupidx = [];
subidix  = [];

for g = 1:length(D.groupname)
    pth    = fullfile(D.path.rootPath, 'results_univariate', D.groupname{g}); % select sub with preprocessed data (from univariate preprocessing)
    subpth = get_files(pth, 'COS*');
    groupidx = [groupidx; repmat(g, size(subpth,1), 1)];
    for s = 1:size(subpth, 1)
        subidix = [subidix, s];
        [~, na, ~] = fileparts(strtrim(subpth(s,:)));
        SubCODE{g}{s,1} = na;
    end
end

%% -------------------------------------------------------------------------
% CREATE RSMS FROM PRE STUDIES DATA
% -------------------------------------------------------------------------
if strcmpi(mode, 'rdm')
    fprintf('-> Computing RDMs and building RSMs (study-level, run once)...\n')
    cos_rsa_build_rdm(D);
    return;   % rien d'autre à faire
end

    
%% -------------------------------------------------------------------------
% SUBJECT LOOP
% -------------------------------------------------------------------------

for subI = S

    % --- Subject identifiers and paths ------------------------------------
    D.subcode  = SubCODE{groupidx(subI)}{subidix(subI)};
    gpcode     = D.groupname{groupidx(subI)};
    subCount   = D.subcode(end-1:end);

    D.sub_dir  = fullfile(D.path.rootPath, 'results_RSA', gpcode, D.subcode);
    D.raw_dir  = fullfile(D.path.rootPath, 'results_univariate', gpcode, D.subcode);
    D.mridir   = fullfile(D.sub_dir, 'anat');
    D.glmpath  = fullfile(D.sub_dir, D.glm_model);
    D.subdata  = fullfile(D.path.behav, 'trialinfo_MRIstudy', sprintf('%s%s_trialinfo.mat', gpcode, subCount));

    fprintf('\n[%s]  Subject: %s  |  Mode: %s\n', datestr(now,'HH:MM:SS'), D.subcode, mode);

    try

        switch lower(mode)
            
            case 'reslice'
                fprintf('  -> Reslicing functional images...\n')
                cos_rsa_reslice(D);
                
                % --------------------------------------------------------------
            case 'spm'
                % First-level GLM: item-wise beta maps for RSA
                fprintf('  -> Estimating first-level GLM...\n')
                mkdir(D.glmpath);
                cos_rsa_estimate_spm(D);

            % --------------------------------------------------------------
            case 'searchlight'
                % Run searchlight RSA using the Viviani toolbox
                fprintf('  -> Running searchlight RSA...\n')
                cos_rsa_searchlight(D);

            % --------------------------------------------------------------
            case 'postprocess'
                % Normalise and smooth individual RSA maps
                fprintf('  -> Post-processing (normalise + smooth)...\n')
                cos_rsa_postprocess(D);
                
            case 'all'
                fprintf('  -> Reslicing functional images...\n')
                cos_rsa_reslice(D);
                fprintf('  -> Estimating first-level GLM...\n')
                mkdir(D.glmpath);
                cos_rsa_estimate_spm(D);
                fprintf('  -> Running searchlight RSA...\n')
                cos_rsa_searchlight(D);
                fprintf('  -> Post-processing (normalise + smooth)...\n')
                cos_rsa_postprocess(D);

            % --------------------------------------------------------------
            otherwise
                error('Unknown mode. Use: ''rdm'', ''reslice'', ''spm'', ''searchlight'', ''postprocess'', or ''all''.');
        end

    catch
        % Log subject-level errors for debugging
        fprintf('  !! Error with %s\n', D.subcode)
        err = lasterror;
        fid = fopen(fullfile(error_dir, sprintf('%s_%s_error.txt', D.subcode, mode)), 'w');
        fprintf(fid, '%s', err.message);
        fclose(fid);
    end

end % subject loop

end