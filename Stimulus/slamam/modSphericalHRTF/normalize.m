angle = 0;%-90:5:90;
distance = [1.5 3];%[0.5 1 1.5 3];

i=0;
for ds=1:length(distance) 
    for ag=1:length(angle)  
    i=i+1;    
    src_name=strcat(num2str(angle(ag)),'d_',num2str(distance(ds)*100),'cm');
    eval(sprintf('load human_hrtf_full_%s',src_name));
    H_OUT=H_OUT/max(max(abs(H_OUT)));
    
    filt=struct('Type','arbitrary','Fs',50000,'b',H_OUT(:,1),'a',ones(size(H_OUT(:,1))));
    filt_name=strcat(src_name,'_L');
    eval(sprintf('save human_hrtf_filt_%s filt',filt_name));
    
    filt=struct('Type','arbitrary','Fs',50000,'b',H_OUT(:,2),'a',ones(size(H_OUT(:,2))));
    filt_name=strcat(src_name,'_R');
    eval(sprintf('save human_hrtf_filt_%s filt',filt_name));
    
    end
end