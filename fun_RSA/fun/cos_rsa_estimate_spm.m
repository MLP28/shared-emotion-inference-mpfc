function cos_rsa_estimate_spm(D)
% COS_RSA_ESTIMATE_SPM  First-level GLM estimation for item-wise RSA (SPM12)
%
% Builds and estimates a first-level GLM where each context item has its own
% regressor. The resulting beta maps (one per item) serve as input patterns
% for the searchlight RSA.
%
% The design includes:
%   - One regressor per context item (120 items: emotional + neutral)
%   - Three nuisance regressors: Neutral response, Context duration, Face duration
%   - Six motion parameters per session (rigid-body realignment)
%   - Session-specific intercepts
%
% Called by cos_rsa_run_firstlevel.m - do not call directly.
%
% INPUT
%   D : configuration structure (set in cos_rsa_run_firstlevel.m)
%       Required fields:
%           D.path.behav   - behavioural data directory
%           D.sub_dir      - subject output directory
%           D.raw_dir      - subject raw data directory
%           D.subdata      - path to subject trialinfo .mat file
%           D.glmpath      - path where SPM.mat will be saved
%           D.nscans       - number of scans per run (cell array)
%           D.filt_rsa     - functional file filter (regexp)
%           D.raw_dir      - univariate results directory
%                            (source of motion parameter files rp_*.txt)
%
% -------------------------------------------------------------------------

%% -------------------------------------------------------------------------
% LOAD BEHAVIOURAL AND STIMULUS FILES
% -------------------------------------------------------------------------
% extract data from pre-study 2 (inference timing) to later include it 
% during creation of fMRI regressors

datapath = D.path.behav;

% Stimulus lists and pre-experimental ratings
[~, ctxlist]  = xlsread(fullfile(datapath, 'infRT_prestudy', 'item_sujet_RT_ctxt_pretests2.xlsx'),    1, 'A2:A121');
[~, listemo]  = xlsread(fullfile(datapath, 'infotask','infoemo.xlsx'),                          'E2:E121');
[rtCxt, ~]    = xlsread(fullfile(datapath, 'infRT_prestudy', 'item_sujet_RT_ctxt_pretests2.xlsx'),    1, 'B2:BE121');

% Load subject trial info (contains onsets, accuracy, etc.)
load(D.subdata);   % loads variable: trial_info

% Emotion labels and category indices
emolabel = {'Colère', 'Surprise', 'Embarras', 'Fierté', 'Neutre'};
idxemo   =  [];

for e = 1:length(emolabel)
    idx             = strcmp(listemo(:,1), emolabel{e});
    idxemo(e,:)   = find(idx);

end


%% -------------------------------------------------------------------------
% HRF AND TIMING PARAMETERS
% -------------------------------------------------------------------------

TR        = D.tr;
nscan     = D.nscans{1};
Ns        = sum(nscan);
dt        = 0.001;              % HRF sampling resolution (1 ms)
fMRI_T    = TR * 1000;          % number of ms per TR
fMRI_T0   = fMRI_T / 2;        % reference time bin (middle of TR)

% Canonical HRF (SPM defaults)
p   = [6 16 1 1 6 0 32];
hrf = spm_hrf(dt, p, fMRI_T);

% Total experiment duration in ms
msduration = round((Ns * TR) * 1000);
starvec    = zeros(1, msduration);

% Session layout: 4 sessions x 30 items
sessorder = [1:30; 31:60; 61:90; 91:120];

% Video durations
[ctxdur,ab]    = xlsread(fullfile(D.path.behav,'infotask','infoemo.xlsx'),'J2:J121'); % extract list of trial duration (inference videos)
facedur = 4000;                          % face duration fixed at 4000 ms

% Item order and session index per trial
itemorder = vertcat(trial_info{:,3});
sesidx    = zeros(size(itemorder));
for ses = 1:4
    sesidx(ismember(itemorder, sessorder(ses,:))) = ses;
end


%% -------------------------------------------------------------------------
% BUILD DESIGN MATRIX REGRESSORS
% -------------------------------------------------------------------------
% Each context item gets its own binary stick function (response-time locked).
% Three additional regressors capture shared variance: neutral response,
% context duration, and face duration.

DC = starvec;   % context duration (shared regressor)
DF = starvec;   % face duration (shared regressor)
SF = [];        % all stick functions [time x regressors]

% --- Emotional context items --------------
P      = repmat(starvec, size(ctxdur,1), 1);   % item x time
onscxt = vertcat(trial_info{:,4});
onsvsg = vertcat(trial_info{:,5});

for i = 1:length(onscxt)
    
    nS        = sesidx(i);
    extratime = sum(nscan(1:nS-1), 2) * TR;
    
    % Context response stick function (response-time locked)
    ons          = onscxt(i) + round(extratime * 1000);
    sp           = rtCxt(i,:) + ons;
    sp(isnan(sp)) = [];
    P(i, sp)     = 1;
    
    % Context duration regressor (shared)
    sp_dur       = repmat(ons, 1, ctxdur(i));
    DC(sp_dur)   = 1;
    
    % Face duration regressor (shared)
    ons_face     = onsvsg(i) + round(extratime * 1000);
    sp_face      = repmat(ons_face, 1, facedur);
    DF(sp_face)  = 1;
    
end

SF = [SF, P'];

% --- Neutral items: response regressor ------------------------------------
NN_P     = starvec;
idxcontx = idxemo(5,:);
onscxt   = vertcat(trial_info{idxcontx,4});

for i = 1:size(idxcontx,2)
    nS        = sesidx(idxcontx(i));
    extratime = sum(nscan(1:nS-1), 2) * TR;
    ons           = onscxt(i) + round(extratime * 1000);
    sp            = rtCxt(idxcontx(i),:) + ons;
    sp(isnan(sp)) = [];
    NN_P(sp)      = 1;
end

SF = [SF, NN_P', DC', DF'];


%% -------------------------------------------------------------------------
% CONVOLVE STICK FUNCTIONS WITH HRF
% -------------------------------------------------------------------------

Xx = zeros(Ns, size(SF,2));

for e = 1:size(SF,2)
    % Smooth to approximate neural signal, then convolve with HRF
    [sdf, ~] = msdf(SF(:,e), 'Gauss', 3000);
    x1       = conv(sdf, hrf);
    x1       = x1((0:(Ns-1)) * fMRI_T + fMRI_T0, :);
    Xx(:,e)  = x1;
end


%% -------------------------------------------------------------------------
% BUILD SPM DESIGN STRUCTURE
% -------------------------------------------------------------------------

SPM = [];
SPM.nscan = nscan;
nses      = length(SPM.nscan);

% --- Item regressors (one per context item) --------------------------------
inc     = 0;
allname = {};

for c = 1:length(ctxlist)
    inc  = inc + 1;
    cname = sprintf('Sn(%d)_cxt%d_img%d_*bf(1)', sesidx(c), c, itemorder(c));
    allname{inc}                       = cname;
    SPM.Sess(1).U(inc).name            = {cname};
    SPM.Sess(1).U(inc).P(1).name       = 'none';
    SPM.Sess(1).Fc(inc).i              = inc;
    SPM.Sess(1).Fc(inc).name           = cname;
    SPM.Sess(1).Fc(inc).p              = 1;
end

% --- Nuisance regressors (neutral response, context duration, face duration)
cnam = {'Neutral', 'Context', 'Face'};
for c = 1:length(cnam)
    inc  = inc + 1;
    allname{inc}                       = cnam{c};
    SPM.Sess(1).U(inc).name            = {cnam{c}};
    SPM.Sess(1).U(inc).P(1).name       = 'none';
    SPM.Sess(1).Fc(inc).i              = inc;
    SPM.Sess(1).Fc(inc).name           = cnam{c};
    SPM.Sess(1).Fc(inc).p              = 1;
end


%% -------------------------------------------------------------------------
% MOTION PARAMETERS AND SESSION EFFECTS
% -------------------------------------------------------------------------

rnam    = {'X','Y','Z','x','y','z'};
M       = [];
Xb      = [];
bkname  = {};
mvtname = {};
Sess    = [];

for ses = 1:nses

    sessdir = fullfile(D.raw_dir, sprintf('run_%d', ses));
    fn      = get_files(sessdir, 'rp_*.txt');
    [r1,r2,r3,r4,r5,r6] = textread(fn, '%f%f%f%f%f%f');
    C       = [r1 r2 r3 r4 r5 r6];

    % Detrended motion parameters (block diagonal)
    M = blkdiag(M, spm_detrend(C));

    % Session intercepts
    Sess(ses).row = size(Xb,1) + (1:nscan(ses));
    B             = ones(nscan(ses), 1);
    Xb            = blkdiag(Xb, B);

    bkname{ses}   = sprintf('SesEffect%d', ses);
    for r = 1:length(rnam)
        mvtname = [mvtname, sprintf('%s%d', rnam{r}, ses)];
    end

end

Xx = [Xx, M];


%% -------------------------------------------------------------------------
% ASSEMBLE AND ESTIMATE SPM MODEL
% -------------------------------------------------------------------------

SPM.xX.X    = [Xx, Xb];
SPM.xX.iH   = [];
SPM.xX.iC   = 1:size(Xx,2);
SPM.xX.iB   = (1:size(Xb,2)) + size(Xx,2);
SPM.xX.iG   = [];
SPM.xX.name = [allname, mvtname, bkname];

% Global normalisation and autocorrelation
SPM.xGX.iGXcalc = 'None';
SPM.xVi.form    = 'AR(1)';

% High-pass filter (one per session)
for ses = 1:nses
    SPM.xX.K(ses).HParam = 128;
end

% Functional image files (one for 4 sessions)
sessdir    = fullfile(D.sub_dir, 'preprocessing');
func_files = spm_select('ExtFPList', sessdir, D.filt_rsa{1});
    
SPM.xY.P  = cat(1, func_files);
SPM.xY.RT = TR;

% Configure design and estimate
SPM.swd = D.glmpath;
cd(SPM.swd)
SPM = get_spmdesign(SPM, Sess);

% Remove pre-existing SPM.mat and mask to avoid SPM dialog prompts
spm_mat_dir = fullfile(D.glmpath, 'SPM.mat');
maskfile    = fullfile(D.glmpath, 'mask.nii');
if exist(spm_mat_dir, 'file'); delete(spm_mat_dir); end
if exist(maskfile,    'file'); delete(maskfile);    end

SPM = spm_spm(SPM);

end