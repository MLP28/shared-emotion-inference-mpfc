function run_ppi(SPM,coord_target,coord_cont,storepath,analysisname,C)
% RUN_PPI  Compute PPI maps for one subject and one seed region
%
% For a given seed ROI and a set of target voxels, this function:
%   1. Builds the confound/filtering structure from the SPM design
%   2. Extracts and averages the seed signal; extracts target signals
%   3. Removes confounds from both signals
%   4. Deconvolves the seed signal to recover the neural response
%   5. Estimates a PPI GLM at each target voxel (ppi_glm)
%   6. Writes one beta map per condition and one contrast map per contrast
%      All maps are smoothed (FWHM 8mm)
%   7. Saves the beta matrix and data structure to disk
%
% Called by cos_run_ppi.m - do not call directly.
%
% INPUT
%   SPM          : SPM structure loaded from subject's first-level SPM.mat
%   coord_target : [3 x N] MNI coordinates (mm) of target voxels
%   coord_cont   : [3 x M] MNI coordinates (mm) of seed ROI voxels
%   storepath    : output directory for PPI maps and .mat file
%   analysisname : label used in output filenames (e.g. 'ppi')
%   C            : contrast structure with fields:
%                      C.conlabel : cell array of contrast names
%                      C.con      : contrast matrix [nContrasts x nConds]
%                      C.condname : cell array of condition names
%
% -------------------------------------------------------------------------
 
%% -------------------------------------------------------------------------
% BUILD CONFOUND AND FILTERING STRUCTURE
% -------------------------------------------------------------------------
 
% prepare counfound and filtering structure
data            = [];
data.RT         = SPM.xY.RT;                % TR (seconds)
data.T          = data.RT*1000;                    % microtime resolution (number of time bins per scan)
data.T0         = data.T/2;                     % reference time bin (middle of TR)
data.dt         = .001;                    % time bin length (1ms)

% Whitening matrix and confound matrix (motion + session effects)
% KWX0: filtered and whitened confound regressors, used to clean the signal
data.W          = SPM.xX.W;
nc              = length(SPM.Sess.U);   % number of task regressors
data.KWX0       = spm_filter(SPM.xX.K,SPM.xX.W*SPM.xX.X(:,nc+1:end));
data.K          = vertcat(SPM.xX.K(:).X0); % high-pass filter basis

smk = [8 8 8]; % smoothing kernel (mm), applied to all output maps
 
% Load stick functions from univariate analysis (SF must be in SPM.swd)
% Each column is a binary onset vector for one condition, used to build
% the psychological regressor of the PPI interaction term
load('SF') 
data.StickFunction = [];
for c = 1:6
    input = SF(:,c);
    [sdf kernel] = msdf(input,'Gauss',3000); % smooth to approximate neural signal
    data.StickFunction(:,c) = [sdf];% TH/IN/NI

end
data.cond2inc   = 1:length(C.condname); % indices of conditions

%% -------------------------------------------------------------------------
% EXTRACT SIGNALS
% -------------------------------------------------------------------------

affmat          = SPM.xY.VY(1).mat;  % voxel-to-mm affine matrix
% Target voxels: convert MNI coords to voxel indices and extract time series
XYZ_target      = mm2voxel(coord_target,affmat,2);
ID_target       = sub2ind(SPM.xY.VY(1).dim,XYZ_target(:,1),XYZ_target(:,2),XYZ_target(:,3));
Ytarget         = spm_get_data(SPM.xY.VY,XYZ_target'); % [nScans x nTargetVoxels]

% Seed ROI: convert MNI coords, extract and average across ROI voxels
XYZ_cont        = mm2voxel(coord_cont,affmat,2);
Ycont           = spm_get_data(SPM.xY.VY,XYZ_cont');
Ycont           = nanmean(Ycont,2); % [nScans x 1] mean seed signal

%% -------------------------------------------------------------------------
% CONFOUND REMOVAL AND DECONVOLUTION
% -------------------------------------------------------------------------
 
% Project out confounds (motion parameters, session effects, HRF filter)
Ycont_c         = remove_confound(SPM,Ycont,nc);
Ytarget_c       = remove_confound(SPM,Ytarget,nc);

% Deconvolve seed signal: recover underlying neural activity
% (removes HRF to get the neural signal driving the BOLD response)
XNcont          = get_deconvolvedsignal(data,Ycont_c);

%% -------------------------------------------------------------------------
% PPI GLM
% -------------------------------------------------------------------------
% For each target voxel, estimate a GLM with:
%   - PPI interaction regressors (deconvolved neural signal x condition stick functions)
%   - Physiological regressor (seed signal)
%   - Psychological regressors (condition stick functions)
% Returns Bmat [nConds x nVoxels]: PPI beta per condition per target voxel

data.Ytarget_c  = Ytarget_c;
data.Ycont_c    = Ycont_c;
data.XNcont     = XNcont;

[Bmat,data] = ppi_glm(data);

%% -------------------------------------------------------------------------
% WRITE CONDITION BETA MAPS
% -------------------------------------------------------------------------
% One map per condition: voxel values = PPI beta for that condition
 
for i = 1:length(data.cond2inc)
    
    % Place beta values at target voxel indices
    I = zeros(SPM.xY.VY(1).dim);
    I(ID_target) = Bmat(i,:);
    
    % Build output volume header from first functional image
    nv          = SPM.xY.VY(1);
    if isfield(nv,'dat')
        nv          = rmfield(nv, {'private', 'dat'});
    else
        nv          = rmfield(nv, {'private'});
    end
    nv.fname    = fullfile(storepath,sprintf('%s_%s.nii',C.condname{i},analysisname));
    nv.descrip  = 'ppi map';
    
    spm_write_vol(nv,I);
    
    % Smooth
    in = nv.fname;
    out = fullfile(storepath, sprintf('%s_%s_sm8.nii',C.condname{i},analysisname));
    spm_smooth(in,out,smk)
end

% Save beta matrix and full data structure for later use
fname    = fullfile(storepath,sprintf('%s.mat',analysisname));
save(fname,'Bmat','data')

%% -------------------------------------------------------------------------
% WRITE CONTRAST MAPS
% -------------------------------------------------------------------------
% One map per contrast: voxel values = linear combination of condition betas
 
for i = 1:length(C.conlabel)
    
    % Apply contrast weights across conditions
    I = zeros(SPM.xY.VY(1).dim);
    I(ID_target) = C.con(i,:)*Bmat(data.cond2inc,:);
    
    % vol info
    nv          = SPM.xY.VY(1);
    if isfield(nv,'dat')
        nv          = rmfield(nv, {'private', 'dat'});
    else
        nv          = rmfield(nv, {'private'});
    end
    nv.fname    = fullfile(storepath,sprintf('%s_%s.nii',C.conlabel{i},analysisname));
    nv.descrip  = 'ppi map';
    
    spm_write_vol(nv,I);
    
    % Smooth
    in = nv.fname;
    out = fullfile(storepath,sprintf('%s_%s_sm8.nii',C.conlabel{i},analysisname));
    spm_smooth(in,out,smk)
end

end
