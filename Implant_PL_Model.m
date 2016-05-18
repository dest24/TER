function [pow_received, Path_Loss, snr]=Implant_PL_Model(d,N_point,5)
%==================================================
% Inputs
% d: distance (mm) 50mm <= d <= 500mm
% N_point: Number of path loss realizations
% model_id:
% 1 - Deep implant to on-body
% 2 - Near-surface implant to on-body
% 3 - Deep implant to Implant
% 4 - Near-surface implant to Implant
%==================================================
% Output
% Path_Loss: N_point realizations
%==================================================
% if(d<0.01 | d>2500)
% error('ERROR: Please specify proper distance');
% end
d0=0.25; % 250mm = 0.25 m
P_Transmtd = -10; %in dBm  = 0.1mW 
Noise_Floor = -100; %in dBm
if(model_id==1)
alpha=4.26;
PL_d0=47.14;
sigma_s=7.85;
elseif(model_id==2)
alpha=4.22;
PL_d0=49.81;
sigma_s=6.81;
elseif(model_id==3)
alpha=6.26;
PL_d0=35.04;
sigma_s=8.18;
elseif(model_id==4)
alpha=4.99;
PL_d0=40.94;
sigma_s=9.05;
end
% Upper and lower limit settings for Shadowing as multiples 
% of the standard deviation (i.e. uplimit & lowlimit)
% 1*Standard Deviation (68%); 2*Standard Deviation (95%);
% 3*Standard Deviation (99.7%)
% (Default Truncation is within 2 standard deviation)
uplimit=2; lowlimit=-uplimit;
myarray = randn(1,N_point);
S_normal=sum(myarray)/N_point;
%disp(S_normal);
index=find(S_normal>uplimit);
S_normal(index)=uplimit;
index=find(S_normal<lowlimit);
S_normal(index)=lowlimit;
Path_Loss = PL_d0+10*alpha*log10(d/d0)+sigma_s*S_normal;
pow_received = P_Transmtd - Path_Loss; 
snr = pow_received - Noise_Floor ; % in dB
disp(pow_received);
end 


