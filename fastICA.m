function [IC,W,U,mu,conv,loop] = fastICA(X,varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

num_ic = size(X,1);
real_f = isreal(X);
alg = 'unit';
max_ite = 50;
tol = 1e-2;
rdct_mode = 0;
contfcn = 1;
batch_size = 0;

if nargin < 1
    error('Invarid input arguments\n');
end

if nargin > 1
    num_ic = varargin{1};
end

if nargin > 2
    alg = varargin{2};
end

if nargin > 3
    max_ite = varargin{3};
end

if nargin > 4
    tol = varargin{4};
end

if nargin>5&&mod(nargin-5,2)==0
    for i=5:2:nargin-1
        switch varargin{i}
            case 'reduction1'
                rdct_mode = 1;
                rdct_param = varargin{i+1};
            case 'reduction2'
                rdct_mode = 2;
                rdct_param = varargin{i+1};
            case 'reduction3'
                rdct_mode = 3;
                rdct_param = varargin{i+1};
            case 'ContrastFcn'
                contfcn = varargin{i+1};
            case 'batch_size'
                batch_size = varargin{i+1};
        end
    end
else
    error('Invarid input arguments\n');
end

if isa(contfcn,'numeric')&&contfcn==1
    if real_f
        g = @(x) (exp(x)-exp(-x))./(exp(x)+exp(-x));
        g_dash = @(x) 4./((exp(x)+exp(-x)).^2);
    else
        g = @(x) 1./(2*sqrt(0.1+x));
        g_dash = @(x) -1./(4*sqrt((0.1+x).^3));
    end
elseif isa(contfcn,'numeric')&&contfcn==2
    if real_f
        g = @(x) x.*exp(-x.*x/2);
        g_dash = @(x) (1-x.*x).*exp(-x.*x/2);
    else
        g = @(x) log(0.1+x);
        g_dash = @(x) -1./((0.1+x).^2);
    end
elseif iscell(contfcn)&&numel(contfcn)==2 ...
    &&isa(contfcn{1},'function_handle')...
    &&isa(contfcn{2},'function_handle')
    g = contfcn{1};
    g_dash = contfcn{2};
else
    if real_f
        g = @(x) (exp(x)-exp(-x))./(exp(x)+exp(-x));
        g_dash = @(x) 4./((exp(x)+exp(-x)).^2);
    else
        g = @(x) 1./(2*sqrt(0.1+x));
        g_dash = @(x) -1./(4*sqrt((0.1+x).^3));
    end
end

% Dimension reduction and whitening
mu = mean(X,2);
X_w = bsxfun(@minus,X,mu);
[U,D] = svd(X_w*X_w'/size(X_w,2),'econ');
ev = diag(D);
switch rdct_mode
    case 1
        lv = cumsum(ev/sum(ev))<=rdct_param;
        U = U(:,lv);
        ev = ev(lv);
    case 2
        lv = (ev/max(ev))>rdct_param;
        U = U(:,lv);
        ev = ev(lv);
    case 3
        U = U(:,1:rdct_param);
        ev = ev(1:rdct_param);
end
U = bsxfun(@times,U,1./sqrt(ev.'));
X_w = U'*X_w;

if num_ic > size(X_w,1)
    num_ic = size(X_w,1);
    
    if rdct_mode < 1
        warning('Number of independent components was forced to be number of observations (%d -> %d).',num_ic,size(X_w,1));
    end
end

if strcmp(alg,'unit')
    if real_f
        [W,conv,loop] = fastICA_unit(X_w,num_ic,g,g_dash,max_ite,tol);
    else
        [W,conv,loop] = cfastICA_unit(X_w,num_ic,g,g_dash,max_ite,tol);
    end
elseif strcmp(alg,'sym')
    if real_f && batch_size>0
        [W,conv,loop] = fastICA_sym_online(X_w,num_ic,g,g_dash,batch_size,max_ite,tol);
    elseif real_f
        [W,conv,loop] = fastICA_sym(X_w,num_ic,g,g_dash,max_ite,tol);
    else 
        [W,conv,loop] = cfastICA_sym(X_w,num_ic,g,g_dash,max_ite,tol);
    end
end

IC = W' * X_w;
W = U*W;
end

function [W,conv,loop] = fastICA_unit(X,num_ic,g,g_dash,max_ite,tol)
num_obs = size(X,1);

W = randn(num_obs,num_ic);
W = bsxfun(@rdivide,W,sqrt(sum(W.*W,1)));

loop = zeros(1,num_ic);
conv = false(1,num_ic);
for i=1:num_ic
    while (1)
        old = W(:,i);
        y = X'*W(:,i);
        W(:,i) = mean(bsxfun(@times,X',g(y)),1)' - mean(g_dash(y))*W(:,i);
        if i>1
            W(:,i) = W(:,i) - sum(bsxfun(@times,W(:,1:i-1),W(:,i)'*W(:,1:i-1)),2);
        end
        W(:,i) = W(:,i)/sqrt(W(:,i)'*W(:,i));
        
        loop(i) = loop(i) + 1;
        
        % stop criterion
        tmp = abs(W(:,i)'*old);
        if tmp<=1+tol && tmp>=1-tol
            conv(i) = true;
            break;
        end
        
        if max_ite>0 && loop(i)>=max_ite
            break;
        end
    end
end
end


function [W,conv,loop] = fastICA_sym(X,num_ic,g,g_dash,max_ite,tol)
num_obs = size(X,1);
num_samp = size(X,2);

W = randn(num_obs,num_ic);
W = decorrelation(W);

loop = 0;

while true
    old = W;
    Y = X'*W;
    beta = mean(Y.*g(Y),1);
    alpha = 1./(beta - mean(g_dash(Y),1));
    
    W = W + W*bsxfun(@times,alpha,squeeze(...
        mean(bsxfun(@times,Y,reshape(g(Y),[num_samp,1,num_ic])),1))...
        -diag(beta));
    
    % Decorrelation
    W = decorrelation(W);
    
%     W2 = W*W';
%     norm_w2 = sqrt(max(sum(abs(W2),1)));
%     W = W/norm_w2;
%     W2 = W*W';
%     
%     while max(abs(W2(:) - I(:)))>tol
%         W = 1.5*W - 0.5*W2*W;
%         W2 = W*W';
%     end
    
    loop = loop + 1;
    
    % stop criterion
    tmp = max(abs(abs(sum(old.*W,1))-1));
    
    if tmp<=tol
        conv = true;
        break;
    end
    
    if max_ite>0&&loop>=max_ite
        conv = false;
        break;
    end
end
end

function [W,conv,loop] = fastICA_sym_online(X,num_ic,g,g_dash,batch_size,max_ite,tol)

num_obs = size(X,1);
num_samp = size(X,2);
num_div = fix(num_samp/batch_size);
ep = 1e-10;

W = randn(num_obs,num_ic);
W = decorrelation(W);

loop = 0;

while true
    Ind = batch(num_samp,batch_size);
    old = W;
    
    for i=1:num_div
        Y = X(:,Ind(:,i))'*W;
        beta = mean(Y.*g(Y),1);
        alpha = 1./(beta - mean(g_dash(Y),1)+ep);
        W = W + W*bsxfun(@times,alpha,squeeze(...
            mean(bsxfun(@times,Y,reshape(g(Y),[batch_size,1,num_ic])),1))...
            -diag(beta));
        
        % Decorrelation
        W = decorrelation(W);
    end
    
    loop = loop + 1;
    
    % stop criterion
    tmp = max(abs(abs(sum(old.*W,1))-1));
    
    if tmp<=tol
        conv = true;
        break;
    end
    
    if max_ite>0&&loop>=max_ite
        conv = false;
        break;
    end
    
end
end


function Ind = batch(num_sample,batch_size)
num_div = fix(num_sample/batch_size);
ind = randperm(num_sample);
Ind = sort(reshape(ind(1:num_div*batch_size),[batch_size,num_div]));
end


function W = decorrelation(W)
[U,D] = svd(W'*W,'econ');
ev = diag(D).^-0.5;
W = W*U*bsxfun(@times,U',ev);
end


function [W,conv,loop] = cfastICA_unit(X,num_ic,g,g_dash,max_ite,tol)
num_obs = size(X,1);

W = randn(num_obs,num_ic)+1i*randn(num_obs,num_ic);
W = bsxfun(@rdivide,W,sqrt(sum(W.*conj(W),1)));

loop = zeros(1,num_ic);
conv = false(1,num_ic);
for i=1:num_ic
    while (1)
        old = W(:,i);
        y = W(:,i)'*X;
        y2 = y.*conj(y);
        W(:,i) = mean(bsxfun(@times,X,conj(y).*g(y2)),2) - mean(g(y2)+y2.*g_dash(y2))*W(:,i);
        if i>1
            W(:,i) = W(:,i) - sum(bsxfun(@times,W(:,1:i-1),(W(:,1:i-1)'*W(:,i)).'),2);
        end
        W(:,i) = W(:,i)/sqrt(W(:,i)'*W(:,i));
        
        loop(i) = loop(i) + 1;
        
        % stop criterion
        tmp = abs(W(:,i)'*old);
        if (tmp <= 1+tol)&&(tmp >= 1-tol)
            conv = true;
            break;
        end
        if max_ite>0&&loop(i)>=max_ite
            break;
        end
    end
end
end

function [W,conv,loop] = cfastICA_sym(X,num_ic,g,g_dash,max_ite,tol)
num_obs = size(X,1);
len = size(X,2);

W = randn(num_obs,num_ic) + 1i*randn(num_obs,num_ic);
[U,D] = svd(W'*W,'econ');
W = W*U*bsxfun(@rdivide,U',sqrt(diag(D)));
% W = randn(num_obs,num_ic)+1i*randn(num_obs,num_ic);
% W = bsxfun(@rdivide,W,sqrt(sum(W.*conj(W),1)));

loop = 0;
while true
    old = W;
    Y = W'*X;
    Y2 = Y.*conj(Y);
    beta = mean(Y2.*g(Y2),2);
    W = W + W*bsxfun(@rdivide,squeeze(mean(bsxfun(@times,...
        reshape(Y,[num_ic 1 len]),reshape(conj(Y).*g(Y2),[1 num_ic len])),3))...
        -diag(beta),(beta-mean(g(Y2)+Y2.*g_dash(Y2),2)).');
    [U,D] = svd(W'*W,'econ');
    W = W*U*bsxfun(@rdivide,U',sqrt(diag(D)));
    
    loop = loop + 1;
    
    % stop criterion
    tmp = max(abs(abs(sum(old.*conj(W),1))-1));
    if tmp<=tol
        conv = true;
        break;
    end
    
    if max_ite>0&&loop>=max_ite
        conv = false;
        break;
    end
end
end

