function nnschalfprior3( X, A, p, fname )
% nnschalfprior - Do hyperspectral unmixing using L1/2 method
%
% SYNTAX:
% nnschalfprior( X, A p, fname )
%
% INPUT:    
% X       the non-negative data (in columns)
% A       initial matrix of A
% fname   name of file to write
% p       algorithm options:
%
% p.sources      number of components
% p.derta        value of derta
% p.lambda       value of lambda
dims = size(X,1);
samples = size(X,2);
sources = p.sources;
derta = p.derta;
lambda = p.lambda;
times = 0;
A_est = A;
% Initializing S
S = p.S;

% These will store the history of the objective function
objhistory = [];
iterhistory = [];

% Loop indefinitely
iter = 0;
tt = 0;
obj = 0.5 * sum(sum((X-A_est*S).^2)) + lambda * sum(sum(sqrt(S)));

while 1,
    objhistory = [objhistory obj];
    iterhistory = [iterhistory iter];
    if rem(iter,100)==0,
    % Update activity measures
    activations = sqrt(mean(S'.^2));
    fprintf(['\nSaving file: ' fname '...']);
%     save(fname,'p','A_est','S','objhistory','iterhistory','activations');
    fprintf('DONE!\n');
    end  
    Xf = [X;derta * ones([1 samples])];
    Af = [A_est;derta * ones([1 sources])]; 
    
    S = S.*(Af'*Xf)./((Af'*Af*S) + 0.5 * lambda ./ sqrt(S+1e-19)+1e-19 );
    A_est = A_est.*(X*S')./(A_est*S*S' );
    objnew = 0.5 * sum(sum((X-A_est*S).^2)) + lambda * sum(sum(sqrt(S))); 
    toltemp = abs(obj - objnew)/obj;
%     if toltemp < 0.0001
    if toltemp < 0.001
        times = times + 1;
    end
    
    obj = objnew;
    error(iter + 1) = obj;
    if rem(iter,100)==0,
        display(['iter: ', num2str(iter), ' error: ', ...
        num2str(obj), ' tol: ', num2str(toltemp)]);
    end
    iter = iter+1;
%     if icclose allter > 1000 || times>=1%%
    if iter > 300 || times>=1%%
        figure;plot(error);
        display(['finnay iter: ', num2str(iter), ' error: ', ...
        num2str(obj), ' tol: ', num2str(toltemp)]);
        fprintf(['\nSaving file: ' fname '...']);
        save(fname,'p','A_est','S','objhistory','iterhistory');
        fprintf('DONE!\n');
        break;
    end
end

