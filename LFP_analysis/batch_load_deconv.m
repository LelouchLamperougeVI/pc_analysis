function deconv=batch_load_deconv(sessions)
deconv=[];
for i=sessions
    cd(num2str(i))
    cd('Plane1')
    temp=load('deconv.mat');
    
    deconv=[deconv;temp.deconv];
    cd ..
    cd ..
end