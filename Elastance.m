function E = Elastance(T,step,HR,end_t)

En = [];
Tc = 60/HR; % Cycle time for one beat;
tmax = 0.2 + 0.1555*Tc; % Find time of systole from cycle time
points = Tc/step; % Number of samples of one cycle

r = ceil(end_t/Tc); % Find number of cycles needed, rounding to upper integer
% r = ceil(length(T)/points); % Find number of cycles needed, rounding to upper integer

for i = 1:points+1
    tn(i) = i*(step/tmax);
    t1 = (tn(i)/0.7)^1.9;
    t2 = 1 + t1;
    t3 = 1 + (tn(i)/1.173474)^21.9;
    Ecyc(i) = 1.553174*(t1/t2)*(1/t3);
end
for j = 1:r+1
    En = [En Ecyc];
end
En = En(1:length(T));
E = En;

% function E = Elastance(step,T)
% 
% En = [];
% 
% HR = 90; 
% Tc = 60/HR;
% tmax = 0.2 + 0.1555*Tc;
% points = Tc/step;
% r = ceil(length(T)/points);
% 
% for j = 1:r+1
%     for i = 1:points
%         tn(i) = i*(step/tmax);
%         t1 = (tn(i)/0.7)^1.9;
%         t2 = 1 + t1;
%         t3 = 1 + (tn(i)/1.173474)^21.9;
%         Ecyc(i) = 1.553174*(t1/t2)*(1/t3);
%     end
%     En(length(En)+1:length(En)+points) = Ecyc;
% end
% 
% En = En(1:length(T));
% E = En;