function deconv=adjacency_matrix(deconv,win)
% Activation at every time points t is defined by all adjacent active time points t-win/t+win

kernel=ones(2*win+1,size(deconv,2));
deconv=mat_conv(deconv,kernel);
deconv([1:win end-win+1:end],:)=[];
deconv=double(logical(deconv));