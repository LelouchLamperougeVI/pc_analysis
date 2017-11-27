function [behavior,deconv,tcs,maskNeurons,mimg]=load_data(animal,date,session,plane)

if nargin<4
    plane='plane1';
else
    plane=['plane' num2str(plane)];
end
% path=['F:\Dun\analysis\' animal '\' date '\' session '\' plane '\'];
path=['X:\scratch\Adam\analysis\' animal '\' date '\' session '\' plane '\'];
deconv=load([path 'deconv.mat']);
maskNeurons=load([path 'masks_neurons.mat']);
tcs=load([path 'timecourses.mat']);
mimg=load([path 'mean_img.mat']);
behavior=load([path '..\behavior.mat']);

deconv=deconv.deconv;
maskNeurons=maskNeurons.maskNeurons;
tcs=tcs.tcs;
behavior=behavior.behavior;
mimg=mimg.mimg;