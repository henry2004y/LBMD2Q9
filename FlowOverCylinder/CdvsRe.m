Re_=[1 5 10 25 50 100 150 200 300 400];
for k=1:length(Re_)
   Re =Re_(k);
   clearvars -except Re k Re_ longstore_Cd longstore_Fx longstore_St
   LB_Project_Auto
   Cd = storage_Fx/0.5/U^2/ceil(2*R);
   longstore_Cd(k) = mean(Cd(200:end));
   longstore_Fx(k) = mean(Fx(200:end));
   %longstore_St(k) = St;
end

plot(Re_,longstore_Cd)