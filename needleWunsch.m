function [eScore,M,trace,traceback]=needleWunsch(seq1,seq2,traceFlag)
% a modified Needleman-Wunsch algorithm for global alignment of replay
% sequences
%
% by HaoRan Chang

M=zeros(length(seq1)+1,length(seq2)+1);
traceback=M;
M(:,1)=0:-1:-length(seq1);
M(1,:)=0:-1:-length(seq2);
traceback(:,1)=2;
traceback(1,:)=1;

simM=getSimM(seq1,seq2);

for i=2:size(M,1)
    for j=2:size(M,2)
        M(i,j)=max([M(i-1,j-1)+simM(i-1,j-1) M(i-1,j)-1 M(i,j-1)-1]);
        arrow=M(i,j)==[M(i-1,j-1)+simM(i-1,j-1) M(i-1,j)-1 M(i,j-1)-1];
        traceback(i,j)=base2dec(replace(num2str(arrow),' ',''),2);
    end
end

eScore=M(end,end);

if nargin==3 && traceFlag
    trace=trackIt(M,traceback);
    trace=unpackCell(trace);
end

function simM=getSimM(seq1,seq2)
% compute similarity matrix
simM=seq1'*seq2;
simM=sqrt(simM);
simM=double(simM==floor(simM));
simM(~simM)=-1;

function trace=trackIt(M,traceback)
% scoring function
scoring=true;
i=size(M,1);
j=size(M,2);
trace=[];
while(scoring)
    switch traceback(i,j)
        case 1
            j=j-1;
            trace=['1' trace];
        case 2
            i=i-1;
            trace=['2' trace];
        case 4
            i=i-1;
            j=j-1;
            trace=['4' trace];
            
        case 3
            recursiveM=M(1:i,1:j); % first time in my life using a recursive function :)
            recursiveTrace={traceback(1:i,1:j),traceback(1:i,1:j)};
            recursiveTrace{1}(i,j)=2;
            recursiveTrace{2}(i,j)=1;
            trace={strcat(trackIt(recursiveM,recursiveTrace{1}),trace),strcat(trackIt(recursiveM,recursiveTrace{2}),trace)};
%             trace=trace+[trackIt(recursiveM,recursiveTrace{1}) trackIt(recursiveM,recursiveTrace{2})]-M(i,j);
            break
        case 5
            recursiveM=M(1:i,1:j);
            recursiveTrace={traceback(1:i,1:j),traceback(1:i,1:j)};
            recursiveTrace{1}(i,j)=4;
            recursiveTrace{2}(i,j)=1;
            trace={strcat(trackIt(recursiveM,recursiveTrace{1}),trace),strcat(trackIt(recursiveM,recursiveTrace{2}),trace)};
            break
        case 6
            recursiveM=M(1:i,1:j);
            recursiveTrace={traceback(1:i,1:j),traceback(1:i,1:j)};
            recursiveTrace{1}(i,j)=4;
            recursiveTrace{2}(i,j)=2;
            trace={strcat(trackIt(recursiveM,recursiveTrace{1}),trace),strcat(trackIt(recursiveM,recursiveTrace{2}),trace)};
            break
            
        case 7
            recursiveM=M(1:i,1:j);
            recursiveTrace={traceback(1:i,1:j),traceback(1:i,1:j),traceback(1:i,1:j)};
            recursiveTrace{1}(i,j)=4;
            recursiveTrace{2}(i,j)=2;
            recursiveTrace{3}(i,j)=1;
            trace={strcat(trackIt(recursiveM,recursiveTrace{1}),trace),strcat(trackIt(recursiveM,recursiveTrace{2}),trace),strcat(trackIt(recursiveM,recursiveTrace{3}),trace)};
            break
    end
%     trace=trace+M(i,j);
    if i==1 && j==1
        scoring=false;
    end
end

function C=unpackCell(C)
% unpack trace
if ~iscellstr(C)
    C=horzcat(C{:});
    for i=1:length(C)
        if iscell(C{i})
            tmp=unpackCell(C{i});
            C=horzcat(C,tmp);
        end
    end
end
idx=zeros(length(C));
for i=1:length(C)
    idx(i)=iscell(C{i});
end
C(logical(idx))=[];        