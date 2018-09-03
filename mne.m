function [M,snr,ev] = mne(L,num_dip,C,lamda,gamma)
%MNE Minimum Norm Estimates (MNE)
%   Detailed explanation goes here

num_ch = size(C,1);
num_dim = size(L,2)/num_dip;

% Whitening
[U,S] = svd(C,'econ');
ev = diag(S);
W = bsxfun(@times,U,1./sqrt(ev)');
L = W'*L;

% Depth weight
weight = sum(reshape(L.*L,[num_ch*num_dim,num_dip]),1).^-gamma;
weight = sqrt(reshape(repmat(weight,[num_dim,1]),[1,num_dim*num_dip]));

% Solving MNE problem
[U,S,V] = svd(bsxfun(@times,L,weight),'econ');
gamma = diag(S);
snr = lamda^-2*sum(gamma.^2)/num_ch; % Power SNR
gamma = gamma./(gamma.^2+lamda^2);
M = bsxfun(@times,weight',V*bsxfun(@times,gamma,U')*W');

end

