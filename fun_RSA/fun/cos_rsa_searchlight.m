function cos_rsa_searchlight(D)
% COS_RSA_SEARCHLIGHT  Build input structure and run searchlight RSA (Viviani toolbox)
%
% This function constructs the 'args' input structure required by the
% Viviani toolbox function rsa_rsm_searchlight, then launches the
% searchlight analysis.
%
% The args structure specifies:
%   - First-level GLM directories (one per subject)
%   - RSM file (built by cos_rsa_build_rdm.m)
%   - Beta image selection (emotional items only, neutral excluded)
%   - Searchlight parameters (radius, method, brain map type)
%   - Output directory
%
% Note: the Viviani toolbox requires a lowercase spm.mat file. If only
% SPM.mat exists in the GLM directory, a copy is created automatically.
%
% Called by cos_rsa_run_firstlevel.m : do not call directly.
%
% INPUT
%   D : configuration structure (set in cos_rsa_run_firstlevel.m)
%       Required fields:
%           D.path.behav                  - behavioural data directory
%           D.path.rsm                    - directory containing rsm_compromise.mat
%           D.sub_dir                     - subject output directory
%           D.glm_model                   - GLM folder name
%           D.searchlight.size            - searchlight radius (mm)
%           D.searchlight.method          - similarity method ('Pearson', etc.)
%           D.searchlight.brainmap_type   - pattern type ('cov', 'sscp', 'cor')
%           D.searchlight.model_name      - label for output filenames
%           D.searchlight.rsm_of_interest - cell array of RSM names to test
%           D.searchlight.rsm_partial     - cell array of RSM names for partial corr
%
% -------------------------------------------------------------------------


%% -------------------------------------------------------------------------
% IDENTIFY NON-NEUTRAL BETA INDICES
% -------------------------------------------------------------------------
% The GLM contains one beta per item (120 total).
% Neutral items (category 5) are excluded from the RSA pattern matrix.

[~, listemo] = xlsread(fullfile(D.path.behav, 'infotask','infoemo.xlsx'), 'E2:E121');

emolabel = {'Colère', 'Surprise', 'Embarras', 'Fierté', 'Neutre'};
idxemo   = zeros(length(emolabel), 24);
for e = 1:length(emolabel)
    idx        = strcmp(listemo(:,1), emolabel{e});
    idxemo(e,:) = find(idx);
end

% beta_idx: linear indices of non-neutral items (used to select beta images)
neut_idx = idxemo(5,:);                      % row indices of neutral items
all_idx  = 1:size(listemo,1);
beta_idx = setdiff(all_idx, neut_idx)';      % keep only emotional items


%% -------------------------------------------------------------------------
% BUILD ARGS STRUCTURE (Viviani toolbox input)
% -------------------------------------------------------------------------

args = struct();

% First-level GLM directory for this subject
args.Directories = fullfile(D.sub_dir, D.glm_model);

% RSM file (contains all RSMs built by cos_rsa_build_rdm)
args.MapFile = fullfile(D.path.rsm, 'rsm_compromise.mat');
load(args.MapFile);   % loads: rsm_set

% Pattern similarity method and brain map type
args.Method       = D.searchlight.method;        % 'Pearson' recommended
args.BrainMapType = D.searchlight.brainmap_type; % 'cov' recommended

% RSM selection: primary contrast and partial correlation regressors
args.MapSel       = D.searchlight.rsm_of_interest;
args.MapSelPcorr  = D.searchlight.rsm_partial;

% Beta image selection: indices of emotional items only
args.BetaIdx = beta_idx;
args.BetaFiles = [];   % auto-selected from args.Directories

% Mask: use SPM-generated mask.nii from first-level directory
args.MaskFile     = [];   % [] = use mask.nii in each subject's GLM directory
args.MaskConfound = [];   % specific mask for BCov/SCov (optional)

% Searchlight definition
args.SearchlightDef  = 'sphere';
args.SearchlightSize = D.searchlight.size;

% Off-diagonal offset (0 = standard; requires betas in temporal order otherwise)
args.OffDiagOffset = 0;

% Output naming and directory
args.ModelName  = D.searchlight.model_name;
args.OutputDir  = fullfile(D.sub_dir, 'searchlight', sprintf('rsa_%s_%dmm',  D.searchlight.model_name, D.searchlight.size));
mkdir(args.OutputDir);

%% -------------------------------------------------------------------------
% CREATE spm.mat
% -------------------------------------------------------------------------
% Viviani toolbox requires lowercase spm.mat : create a copy if needed
spm_upper = fullfile(D.sub_dir, D.glm_model, 'SPM.mat');
spm_lower = fullfile(D.sub_dir, D.glm_model, 'spm.mat');
if exist(spm_upper, 'file') && ~exist(spm_lower, 'file')
    copyfile(spm_upper, spm_lower);
end
%% -------------------------------------------------------------------------
% RUN SEARCHLIGHT RSA
% -------------------------------------------------------------------------

fprintf('  Searchlight radius : %d mm\n',  D.searchlight.size);
fprintf('  Method             : %s\n',     D.searchlight.method);
fprintf('  RSM of interest    : %s\n',     strjoin(D.searchlight.rsm_of_interest, ', '));
fprintf('  Output directory   : %s\n',     args.OutputDir);

rsa_rsm_searchlight(args);

end