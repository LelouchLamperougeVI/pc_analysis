function h = violin(X, varargin)
% Make a violin plot
% by HaoRan Chang, PhD candidate
% Polaris Research Group, Canadian Centre for Behavioural Neuroscience
% University of Lethbridge, Alberta, Canada

h = figure;
patch([d -d(end:-1:1)], [x x(end:-1:1)], 'red')