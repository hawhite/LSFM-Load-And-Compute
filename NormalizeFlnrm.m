function output = NormalizeFlnrm(act)

% Normalizes activity tracing of N worms 
% input act is array of neurons,timepoints,worms(epocs)
% calculation DeltF/Fo for each neuron, Fo is mean value for that neuron
% 


act=double(act);
lthnrm=length(act(:,1,1));
lthfrm=length(act(1,:,1));
lthep=length(act(1,1,:));


for ep=1:lthep 
    for nrm=1:lthnrm
    nrmsort=sort(act(nrm,:,ep));
    low=nrmsort(1:floor(lthfrm*0.01));  %cal Fo as mean of lowest 1% of values 
    Fo=mean(low);
    mean(act(nrm,:,ep));
    act(nrm,:,ep)=(act(nrm,:,ep)/Fo)-1; % normalize by Fo
    end

end


output= act;    
end
