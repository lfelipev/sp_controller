function [u_k,yD_k,k_new,Useq_new] = fcn(Yseq,y_k,Useq,k,gamma)
%#codegen

% number of elements in sequence
Nseq = numel(Yseq);

% counter
k_new = k+1;
% reset if we've been round the sequence
if k_new>Nseq,
    k_new=1;
end

% extract current demand signal
yD_k = Yseq(k);

% and the current control
u_k = Useq(k);

% form the error
e_k = yD_k - y_k;

% copy the control sequence over
Useq_new = Useq;

% which entry to update
k_learn = k-1;
if k_learn<1,
    k_learn=Nseq;
end

% do the learning update
Useq_new(k_learn) = Useq_new(k_learn) + gamma*e_k;

