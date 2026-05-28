
time1 = load('timesmoothPlusDam_order2_lambda0p1_nu1p0.csv');
time2 = load('timesmoothPlusDam_order5_lambda0p1_nu1p0.csv');
time3 = load('timesmoothPlusDam_adaptive_couplingPBC_nu1p0_lambda0p1.csv');
time4 = load('timesmoothPlusDam_adaptive_couplingCIF_nu1p0_lambda0p1.csv');

% time1 = load('timesmoothPlusDam_order1_lambda0p1_nu0p1.csv');
% time2 = load('timesmoothPlusDam_order5_lambda0p1_nu0p1.csv');
% time3 = load('timesmoothPlusDam_adaptive_couplingPBC_nu0p1_lambda0p1.csv');
% time4 = load('timesmoothPlusDam_adaptive_couplingCIF_nu0p1_lambda0p1.csv');
% 
% time1 = load('timesmoothPlusDam_order1_lambda1p0_nu0p1.csv');
% time2 = load('timesmoothPlusDam_order5_lambda1p0_nu0p1.csv');
% time3 = load('timesmoothPlusDam_adaptive_couplingPBC_nu0p1_lambda1p0.csv');
% time4 = load('timesmoothPlusDam_adaptive_couplingCIF_nu0p1_lambda1p0.csv');

times_rel = [time1,time2,time3,time4]/time2;

times_rel