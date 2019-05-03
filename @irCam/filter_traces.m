function filter_traces(obj,freq)
%filter the ROI traces to remove low frequency fluctuations

if nargin < 2
    freq = .5;
end

[a,b,c,d]=butter(2,freq/obj.cam.FrameRate*2,'high');
[sos,g]=ss2sos(a,b,c,d);
obj.traces=zscore( filtfilt(sos,g,obj.original_traces) );