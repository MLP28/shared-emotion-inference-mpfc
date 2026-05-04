function [Bmat,data] = ppi_glm(data)

% This function calculate an interaction term between each condition stick
% functions and neural predictors (bold deconvolved neural activity and optimal control energy).
% This interaction term is then convolved with the HRF to create PPI
% regressors for .

% Then a new GLM is created and estimated to estimate the relationship
% between control and targets nodes (assumed here to be simply nodes not
% defined as control) using PPI regressors derived from bold deconvolution
% and control energy

% In addition to PPI regressors, this new GLM include a standard main effect of source region bold
% signal and a main effect of "psychological" conditions.

% V3 build 2 distincts GLM for XN & EN while V2 combines those 2 predictors
% into a single GLM

% Setup variables
%--------------------------------------------------------------------------

RT              = data.RT;
dt              = data.dt;
NT              = RT/dt;
fMRI_T0         = data.T0;
N               = size(data.Ytarget_c,1);
k               = 1:NT:N*NT;
U               = data.StickFunction; % onset events
Cond            = data.cond2inc;


% HRF
%--------------------------------------------------------------------------
hrf = spm_hrf(dt);

% Source activity (Physiological regressor)
PHYS = data.Ycont_c;

Bmat = [];

% Psychological and PPI regressor
PSY     = zeros(N*NT,length(Cond));
PSYHRF  = [];
PPI_XN  = [];
for c = 1:length(Cond)
    
    % Psychological regressor
    PSY(:,c)    = spm_detrend(full(U(1:end,Cond(c))));
    psyconv     = conv(PSY(:,c) ,hrf);
    PSYHRF(:,c) = psyconv((k-1) + fMRI_T0);
    
    % PPI-Bold deconvolded (XN)
    PSYxn       = PSY(:,c).*data.XNcont;
    ppi         = conv(PSYxn,hrf);
    ppi         = ppi((k-1) + fMRI_T0);
    PPI_XN(:,c) = spm_detrend(ppi);
    
end

% GLM XN
X       = [PPI_XN PSYHRF PHYS ones(size(PHYS,1),1)];
xX      = spm_sp('Set',X);
xX.pX   = spm_sp('x-',xX);
Bmat    = xX.pX*data.Ytarget_c;
data.X  = X;
    
           




