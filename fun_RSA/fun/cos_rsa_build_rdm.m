function cos_rsa_build_rdm(D)
% COS_RSA_BUILD_RDM  Compute confound RDMs and build RSM structure (Viviani toolbox format)
%
% This function computes four representational dissimilarity matrices (RDMs)
% from pre-experimental behavioural data, and converts them into the RSM
% structure required by the Viviani searchlight toolbox.
%
% RDMs computed:
%   1. Inference timing (RT)        - control: response time similarity
%                                     (from prestudy 2)
%   2. Emotional competition        - RDM of interest: cosine distance between
%                                     competing emotion profiles across items
%   3. Dominant emotion             - control: intensity of the dominant emotion
%   4. Interference of competing    - RDM of interest: sum of competing emotion
%                                     intensities
%
% The RSM structure is saved to D.path.rsm for use in cos_rsa_searchlight.m.
%
% Called by cos_rsa_run_firstlevel.m - do not call directly.
%
% INPUT
%   D : configuration structure (set in cos_rsa_run_firstlevel.m)
%       Required fields:
%           D.path.behav - behavioural data directory
%           D.path.rsm   - output directory for RSM .mat file
%
% OUTPUT
%   RSM .mat file saved to D.path.rsm/rsm_compromise.mat
%
% -------------------------------------------------------------------------

datapath  = D.path.behav;
mkdir(D.path.rsm);


%% -------------------------------------------------------------------------
% LOAD STIMULUS INFORMATION
% -------------------------------------------------------------------------

[~, listemo] = xlsread(fullfile(datapath,'infotask', 'infoemo.xlsx'), 'E2:E121');
[rtCxt, ~]   = xlsread(fullfile(datapath, 'infRT_prestudy','item_sujet_RT_ctxt_pretests2.xlsx'), 1, 'B2:BE121');

% Emotion category indices (one row per emotion, columns = item indices)
emolabel = {'Colère', 'Surprise', 'Embarras', 'Fierté', 'Neutre'};
idxemo   = [];  
for e = 1:length(emolabel)
    idx          = strcmp(listemo(:,1), emolabel{e});
    idxemo(e,:)  = find(idx);
end

% Individual emotion intensity ratings from pre-experiment
ctxfile  = get_files(fullfile(datapath,'emointensity_prestudy'), '*ctxt.mat');
emoindiv = [];
for suj = 1:size(ctxfile,1)
    load(strtrim(ctxfile(suj,:)));   % loads: emo_ctxt_120
    emoindiv(:,:,suj) = emo_ctxt_120;
end

% Mapping between pilot emotion order and canonical order
emocat         = repmat((1:5)', 1, 24);
orderemo_pilot = [4 1 5 2 3];   % pilot order: Fear, Anger, Neutral, Surprise, Embarrassment


%% -------------------------------------------------------------------------
% RDM 1 : INFERENCE TIMING (RT)
% -------------------------------------------------------------------------
% Dissimilarity based on response time profiles across items.
% Neutral items are excluded. RTs are normalised by context video duration.

% Load context durations and remove neutral items
[cxtDuration, ~] = xlsread(fullfile(datapath,'infotask', 'infoemo.xlsx'), 1, 'J2:J121');
rtCxt_clean      = rtCxt;
rtCxt_clean(idxemo(5,:), :) = [];
cxtDuration(idxemo(5,:))    = [];

% Normalise by duration
rtCxt_clean = rtCxt_clean ./ cxtDuration;

distMats = [];
count    = 0;
for s = 1:size(rtCxt_clean, 2)
    rtVec = rtCxt_clean(:,s);
    rtVec = zscore(rtVec, 0, 'omitnan');
    % Exclude participants with too many missing values
    if sum(isnan(rtVec)) < 48
        count            = count + 1;
        distMats(:,:,count) = squareform(pdist(rtVec, 'euclidean'));
    end
end

% Compromise RDM across participants (MDS-based, see getcompromiseRDM)
[compromiseF, ~, ~, ~] = getcompromiseRDM(distMats, 1);
rdm_RT = squareform(pdist(compromiseF(:,1:2)));


%% -------------------------------------------------------------------------
% RDM 2 : EMOTIONAL COMPETITION
% -------------------------------------------------------------------------
% Cosine distance between competing emotion profiles of pairs of items.
% Only emotional (non-neutral) items are included.

C   = [];
inc = 0;

for i = 1:120
    idimg   = find(idxemo(:,:) == i);
    img_emo = find(orderemo_pilot == emocat(idimg));
    if img_emo ~= 3   % exclude neutral (category 3)
        inc        = inc + 1;
        % Extract ratings for all 4 emotional categories (fixed order: 1 2 4 5)
        comp       = squeeze(emoindiv(i, [1 2 4 5], :));
        dom_idx    = find([1 2 4 5] == img_emo);
        comp(dom_idx,:) = NaN;   % set dominant emotion to NaN
        C(inc,:,:) = comp;
    end
end

nimg      = size(C,1);
np        = nchoosek(1:nimg, 2);   % all pairwise item combinations
rdm_indiv = zeros(nimg, nimg, size(C,3));
compprop  = zeros(1, size(C,3));

% construct RDM
warning('off');
for s = 1:size(C,3)
    resp        = C(:,:,s);
    compprop(s) = mean(nansum(resp,2) > 0);
    rdm         = zeros(nimg, nimg);
    for p = 1:size(np,1)
        jug1  = C(np(p,1),:,s);
        jug2  = C(np(p,2),:,s);
        comid = ~isnan(sum([jug1; jug2], 1));
        jug1  = jug1(comid);
        jug2  = jug2(comid);
        rdm(np(p,1), np(p,2)) = pdist([jug1; jug2], 'cosine');
        rdm(np(p,2), np(p,1)) = rdm(np(p,1), np(p,2));
    end
    rdm_indiv(:,:,s) = rdm;
end
warning('on');

% Exclude participants with too few competitive responses
torm = compprop < 0.1;
rdm_indiv(:,:,torm) = [];

[compromiseF, ~, ~, ~] = getcompromiseRDM(rdm_indiv, 1);
rdm_comp = squareform(pdist(compromiseF(:,1:2)));


%% -------------------------------------------------------------------------
% RDM 3 : DOMINANT EMOTION INTENSITY
% -------------------------------------------------------------------------
% Euclidean distance between dominant emotion intensity values across items.

C   = [];
inc = 0;

for i = 1:120
    idimg   = find(idxemo(:,:) == i);
    img_emo = find(orderemo_pilot == emocat(idimg));
    if img_emo ~= 3 % neutral
        inc        = inc + 1;
        C(inc,:,:) = squeeze(emoindiv(i, img_emo, :));
    end
end

rdm_indiv = zeros(size(C,1), size(C,1), size(C,3));
for s = 1:size(C,3)
    jug                = C(:,s);
    rdm_indiv(:,:,s)   = squareform(pdist(jug, 'euclidean'));
end

[compromiseF, ~, ~, ~] = getcompromiseRDM(rdm_indiv, 1);
rdm_domint = squareform(pdist(compromiseF(:,1:2)));


%% -------------------------------------------------------------------------
% RDM 4 : SUM OF COMPETING EMOTION INTENSITIES
% -------------------------------------------------------------------------
% Euclidean distance between total competing emotion load across items.

C   = [];
inc = 0;

for i = 1:120
    idimg    = find(idxemo(:,:) == i);
    img_emo  = find(orderemo_pilot == emocat(idimg));
    otheremo = setdiff([1 2 4 5], img_emo);
    if img_emo ~= 3
        inc        = inc + 1;
        C(inc,:,:) = squeeze(emoindiv(i, otheremo, :));
    end
end

Csum      = squeeze(nansum(C, 2));   % sum competing intensities across categories
rdm_indiv = zeros(size(Csum,1), size(Csum,1), size(Csum,2));
prop      = zeros(1, size(C,3));

for s = 1:size(C,3)
    jug              = Csum(:,s);
    prop(s)          = sum(nansum(jug,2) > 0);
    rdm_indiv(:,:,s) = squareform(pdist(jug, 'euclidean'));
end

% Exclude participants with too few observations
rdm_indiv(:,:, prop < 10) = [];

[compromiseF, ~, ~, ~] = getcompromiseRDM(rdm_indiv, 1);
rdm_sumcomp = squareform(pdist(compromiseF(:,1:2)));


%% -------------------------------------------------------------------------
% BUILD RSM STRUCTURE (Viviani toolbox format)
% -------------------------------------------------------------------------
% RSM = 1 - normalised RDM (similarity, bounded [0,1])

rsm_RT      = 1 - rdm_RT      ./ max(rdm_RT(:));
rsm_comp    = 1 - rdm_comp    ./ max(rdm_comp(:));
rsm_domint  = 1 - rdm_domint  ./ max(rdm_domint(:));
rsm_sumcomp = 1 - rdm_sumcomp ./ max(rdm_sumcomp(:));


rsm_set = struct();

rsm_set(1).name   = 'compromise_sumcomp';        % RDM of interest: sum of competing emotions
rsm_set(1).type   = 'similarity';
rsm_set(1).RSM    = rsm_sumcomp;
rsm_set(1).idx    = [];
rsm_set(1).data   = [];
rsm_set(1).distf  = 'euclidean';

rsm_set(2).name   = 'compromise_compet_cos';     % RDM of interest: emotional competition
rsm_set(2).type   = 'similarity';
rsm_set(2).RSM    = rsm_comp;
rsm_set(2).idx    = [];
rsm_set(2).data   = [];
rsm_set(2).distf  = 'cosine';

rsm_set(3).name   = 'compromise_domint';               % control: dominant emotion intensity
rsm_set(3).type   = 'similarity';
rsm_set(3).RSM    = rsm_domint;
rsm_set(3).idx    = [];
rsm_set(3).data   = [];
rsm_set(3).distf  = 'euclidean';

rsm_set(4).name   = 'compromise_norm_RT';              % control: inference timing
rsm_set(4).type   = 'similarity';
rsm_set(4).RSM    = rsm_RT;
rsm_set(4).idx    = [];
rsm_set(4).data   = [];
rsm_set(4).distf  = 'euclidean';


%% -------------------------------------------------------------------------
% SAVE
% -------------------------------------------------------------------------

outfile = fullfile(D.path.rsm, 'rsm_compromise.mat');
save(outfile, 'rsm_set');
fprintf('  RSM saved to: %s\n', outfile);

end