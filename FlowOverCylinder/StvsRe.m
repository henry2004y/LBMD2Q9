Re_=[30 50 75 100 125 150 175 200 250 300];
for k=1:length(Re_)
   Re =Re_(k);
   clearvars -except Re k Re_ longstore_St longstore_v
   LB_Project_Auto
   N = 512;
   dfreqs = periodogram(storage_vprobe(end-N-1:end)-mean(storage_vprobe(end-N-1:end)),[],[],1);
   [trash, index] = max(dfreqs);
   freq = index/N;
   longstore_St(k) = freq*10/0.1;
end

plot(Re_,longstore_St)