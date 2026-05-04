function Yc = remove_confound(SPM,Y,nc)

%
% filter and whithen data, and remove confound
% --------------------------------------------

% 1) filter and Withen data (KWY)
KWY            = spm_filter(SPM.xX.K,SPM.xX.W*Y);

% 2) filter and whithening confound (already filtered in create_counfound)
X0      = SPM.xX.X(:,nc:end);
KWX0 = spm_filter(SPM.xX.K,SPM.xX.W*X0);

% 3) add filter to confound
K   = vertcat(SPM.xX.K(:).X0);
X0w = [K,KWX0];

% 4) remove confound from data
Yc    = KWY - X0w*inv(X0w'*X0w)*X0w'*KWY;
