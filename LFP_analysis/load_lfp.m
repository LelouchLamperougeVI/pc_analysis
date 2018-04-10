function [signal,fs]=load_lfp(ch,fn,path)
if nargin==0
    ch=7;
end
if nargin<2
    [fn,path]=uigetfile('*.abf','Select ABF file');
end
[signal,fs]=abfload([path fn]);
signal=signal(:,ch);
fs=fs/1000000; %sec
fs=1/fs; %hz