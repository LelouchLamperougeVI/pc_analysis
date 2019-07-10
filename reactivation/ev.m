function [EV,REV, e1, e2] = ev(varargin)
% Run explained variance analysis
% Takes as input variables of types 'basic_ensembles', 'deconv' and NxN
% corr matrix
% If all 3 epochs are provided as 'deconv', need to specify Gaussian 
% smoothing parameter ['sig', value] and sampling rate ['fs', value]
% [EV, REV] = ev(pre, exp, post, 'sig', val, 'fs', val)

[Rpre, Rexp, Rpost, clust1, clust2] = parse_inputs(varargin);
[EV, REV] = calc_ev(Rpre, Rexp, Rpost);

e1 = zeros(length(clust1), 2);
for i = 1:length(clust1)
    [e1(i,1), e1(i,2)] = calc_ev(Rpre(clust1{i},clust1{i}), Rexp(clust1{i},clust1{i}), Rpost(clust1{i},clust1{i}));
end
e2 = zeros(length(clust2), 2);
for i = 1:length(clust2)
    [e2(i,1), e2(i,2)] = calc_ev(Rpre(clust2{i},clust2{i}), Rexp(clust2{i},clust2{i}), Rpost(clust2{i},clust2{i}));
end



function [EV, REV] = calc_ev(Rpre, Rexp, Rpost)
Rpre = triu(Rpre,1); Rpre = Rpre(Rpre~=0);
Rexp = triu(Rexp,1); Rexp = Rexp(Rexp~=0);
Rpost = triu(Rpost,1); Rpost = Rpost(Rpost~=0);

PrePost = corr(Rpre, Rpost);
ExpPre = corr(Rexp, Rpre);
ExpPost = corr(Rexp, Rpost);

EV = ( (ExpPost - ExpPre*PrePost) / (sqrt((1 - ExpPre^2) * (1 - PrePost^2))) )^2;
REV = ( (ExpPre - ExpPost*PrePost) / (sqrt((1 - ExpPost^2) * (1 - PrePost^2))) )^2;


function [Rpre, Rexp, Rpost, clust1, clust2] = parse_inputs(inputs)
if length(inputs) < 3
    error('You need to provide at least 3 variables: pre, exp, post');
end

sig=[];
fs=[];

epochs = cell(3,1);
for ii = 1:3
    epochs{ii} = idtype(inputs{ii});
    if length(epochs{ii})>1
        if ~isempty(sig)
            if abs(epochs{ii}{2} - sig) > .1 || abs(epochs{ii}{3} - fs) > .1
                error('The ''sig'' or ''fs'' values are inconsistent. Don''t compare apples with oranges!');
            end
        end
        sig = epochs{ii}{2};
        fs = epochs{ii}{3};
        if ii == 1
            clust1 = epochs{ii}{4};
        elseif ii == 3
            clust2 = epochs{ii}{4};
        end
    end
    epochs{ii} = epochs{ii}{1};
end

idx=4;
while(idx<length(inputs))
    switch lower(inputs{idx})
        case 'sig'
            idx=idx+1;
            sig=inputs{idx};
        case 'fs'
            idx=idx+1;
            fs=inputs{idx};
        otherwise
            error(['''' inputs{idx} ''' is not a valid parameter']);
    end
    idx=idx+1;
end

for ii = 1:3
    if isempty(epochs{ii})
        epochs{ii} = dec2corr(inputs{ii}, sig, fs);
    end
end
Rpre = epochs{1};
Rexp = epochs{2};
Rpost = epochs{3};


function rval = idtype(input)
issym = @(A) ( sum(sum( A == A' )) / numel(A) ) > .3; %computers are only as precise as themselves... values returned from corr() don't always return TRUE in issymmetric()

switch class(input)
    case 'ensemble'
        try
            rval{1} = input.R;
            rval{2} = input.ops.sig;
            rval{3} = input.twop.fs;
            rval{4} = input.clust;
        catch
            if isempty(input.R)
                error('The ''basic_ensemble'' does not have a correlation matrix constructed');
            else
                error('Shit, an undefined error has occured :(');
            end
        end
    case 'double'
        if min(size(input)) < 2
            error('''Deconv'' or ''corr'' must be a matrix');
        end
        if ~ismatrix(input)
            error('Input matrix must not contain more than 2 dimensions');
        end
        if size(input,1) == size(input,2) && issym(input)
            rval{1} = input;
        else
            rval{1} = [];
        end
    otherwise
        error('The input variable is neither of class ''double'', nor a ''basic_ensembles''')
end


function R = dec2corr(deconv,sig,fs)
if isempty(sig) || isempty(fs)
    error('You need to manually provide ''sig'' and ''fs''');
end
deconv=(deconv-mean(deconv,'omitnan'))./std(deconv,'omitnan'); %zscore
deconv=fast_smooth(deconv,sig*fs);
deconv(isnan(sum(deconv,2)),:)=[];

R=corr(deconv);