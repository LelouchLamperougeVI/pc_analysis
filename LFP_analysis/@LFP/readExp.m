function readExp(obj, fn)
% Parse the Experiment.xml file

fields = {'date', 'magnification', 'zstage', 'lsm', 'pmt', 'streaming', 'comments'};
s = parseXML(fn, fields);

obj.session.exp.date = datetime(s.ThorImageExperiment.Date.Attributes.date, 'inputformat', 'MM/dd/yyyy HH:mm:ss');
obj.session.objective.magnification = str2double(s.ThorImageExperiment.Magnification.Attributes.mag);
obj.session.objective.refracI = str2double(s.ThorImageExperiment.Magnification.Attributes.indexOfRefraction);
obj.session.comments = s.ThorImageExperiment.Comments.Attributes.text;

obj.session.pmt.enabled = str2double({s.ThorImageExperiment.PMT.Attributes.enableA, s.ThorImageExperiment.PMT.Attributes.enableB, s.ThorImageExperiment.PMT.Attributes.enableC, s.ThorImageExperiment.PMT.Attributes.enableD});
obj.session.pmt.gain = str2double({s.ThorImageExperiment.PMT.Attributes.gainA, s.ThorImageExperiment.PMT.Attributes.gainB, s.ThorImageExperiment.PMT.Attributes.gainC, s.ThorImageExperiment.PMT.Attributes.gainD});

obj.twop.planes.stepsize = str2double(s.ThorImageExperiment.ZStage.Attributes.stepSizeUM);

numPlanes = str2double(s.ThorImageExperiment.ZStage.Attributes.steps) + str2double(s.ThorImageExperiment.Streaming.Attributes.flybackFrames);
if obj.twop.planes.numplanes ~= numPlanes
    error(['There were ' num2str(obj.twop.planes.numplanes) ' planes inside the data folder. However, ' num2str(numPlanes) ' planes were identified in the Experiment.xml file. Please check consistency.']);
end

fs = str2double(s.ThorImageExperiment.LSM.Attributes.frameRate) / numPlanes;
if abs(fs - obj.twop.fs) > .1
    error('The framerate detected from the abf file is significantly different from the framerate specified in the Experiment.xml file. Please check consistency.');
end

if length(obj.twop.planes.planes) > numPlanes - str2double(s.ThorImageExperiment.Streaming.Attributes.flybackFrames)
    warning('You''re analysing more planes than the total number of imaged planes excluding flybacks.');
end

obj.topo.FOV = [str2double(s.ThorImageExperiment.LSM.Attributes.widthUM); str2double(s.ThorImageExperiment.LSM.Attributes.heightUM)];

if ~all([str2double(s.ThorImageExperiment.LSM.Attributes.pixelX) str2double(s.ThorImageExperiment.LSM.Attributes.pixelY)]...
        == [size(obj.topo.mimg, 1) size(obj.topo.mimg,2)])
    error('The pixel count of the loaded masks does not match that indicated in the Experiment.xml file.');
end