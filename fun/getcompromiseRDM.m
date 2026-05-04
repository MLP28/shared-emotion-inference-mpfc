function [compromiseF,compromise,RF,G] = getcompromiseRDM(rdm,type)

% type 1 = distance matrix; type 2 = corr
nimg    = size(rdm,1);
njuges  = size(rdm,3);
m       = ones(nimg,1)./nimg;
cent    = eye(nimg,nimg)-ones(nimg,1)*m';
Dsum    = zeros(nimg,nimg);

for k=1:njuges;
    D      = rdm(:,:,k);
    Dsum   = Dsum+D;

    % fix when nan values 
    D(find(isnan(D))) = 1;
    
    % transform D -> SCP
    if type == 1
        r_S=-(1/2)*cent*D*cent;
    elseif type == 2
        r_S=(1/2)*cent*D*cent;
    end
    
    [PS,lS]=eigen(r_S);
    % normalise by first eig
    S=r_S.*(lS(1).^(-1));% normalize
    eval(['S',int2str(k),'=S;']);
    eval(['r_S',int2str(k),'=r_S;']);
end

% Compute the RV coefficients
% --------------------------------------------------------------------------

ncells     = nimg*nimg;
Grosse_mat = zeros(ncells,njuges);
for k=1:njuges;
    eval(['la_c=S',int2str(k),'(:);']);
    Grosse_mat(:,k)=la_c;
end
num_RV  = Grosse_mat'*Grosse_mat;
la_diag = diag(num_RV);
nd      = length(la_diag);
den_RV  = (repmat(la_diag',nd,1).*repmat(la_diag,1,nd)).^(1/2);
mat_RV  = num_RV./den_RV;

% PCA of the RV matrix
% -----------------------------------------------------------------------------

[P,phi] = eigen(mat_RV);
G       = P*diag(phi.^(1/2));
tau_RV  = round(100 *(phi./sum(phi)));

%  Compute the alpha weights
% -----------------------------------------------------------------------------
weights=P(:,1)/sum(abs(P(:,1)))  ;

% Commpute the compromise
% -----------------------------------------------------------------------------
compromise = zeros(nimg,nimg);
for k=1:njuges;
    eval(['compromise=compromise+weights(k)*S',...
        int2str(k),';'])
end

% PCA of the compromise
% -----------------------------------------------------------------------------
[pc,lc]=eigen(compromise);
compromiseF=pc*diag(sqrt(lc));

% Projection of indiv obs on the compromise
% -----------------------------------------------------------------------------
RProj=pc*diag(lc.^(-1/2)) ;
RF = [];
for k=1:njuges;
    eval(['tmpproj','=S',int2str(k),'* RProj;']);
    
    RF(:,:,k) = tmpproj;
end

