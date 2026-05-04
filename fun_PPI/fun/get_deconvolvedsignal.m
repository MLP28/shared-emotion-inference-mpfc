function XN = get_deconvolvedsignal(data,Yc)


% This function has been adapted from the original spm_peb_ppi.m function
% to deconvolve the bold response and infer underlying neural response


% Setup variables
%--------------------------------------------------------------------------
RT      = data.RT;
dt      = data.dt;
NT      = RT/dt;
Ns      = size(Yc,1);
k       = 1:NT:Ns*NT; % microtime to scan time indices
W       = data.W;
X0      = data.KWX0; % filtered and whitened confounds
M       = size(X0,2);

% Create basis functions and hrf in scan time and microtime
%--------------------------------------------------------------------------
spm('Pointer','watch')
hrf = spm_hrf(dt);


% Create convolved explanatory {Hxb} variables in scan time
%--------------------------------------------------------------------------
xb  = spm_dctmtx(Ns*NT + 128,Ns);
Hxb = zeros(Ns,Ns);
for i = 1:Ns
    Hx       = conv(xb(:,i),hrf);
    Hxb(:,i) = Hx(k + 128);
end
xb = xb(129:end,:);


% BOLD DECONV (from
% https://github.com/Masharipov/TMFC_simulations/blob/main/deconvolution/matlab/bold_deconvolution.m)
%--------------------------------------------------------------------------
alpha   = .005; % alpha parameter for ridge regression
XN      = bold_deconvolution(Yc,RT,alpha,NT,0,xb,Hxb);
XN      = spm_detrend(XN);   



