clear all
close all
clc
format long

relaxation_time = '0p5';
order_class_low = 2;
order_class_high= 12;

smooth_par = 200;
toldown = '0p0003';
tolup = '0p00045';
time = '0p8';

name2 = strcat('shockPlusSmooth','_order',string(2),...
    '_t',time,'_Kn',string(relaxation_time));type2='hme';
name4 = strcat('shockPlusSmooth','_order',string(4),...
    '_t',time,'_Kn',string(relaxation_time));type4='hme';
name6 = strcat('shockPlusSmooth','_order',string(6),...
    '_t',time,'_Kn',string(relaxation_time));type6='hme';
name8 = strcat('shockPlusSmooth','_order',string(8),...
    '_t',time,'_Kn',string(relaxation_time));type8='hme';
name10 = strcat('shockPlusSmooth','_order',string(10),...
    '_t',time,'_Kn',string(relaxation_time));type10='hme';
name12 = strcat('shockPlusSmooth','_order',string(12),...
    '_t',time,'_Kn',string(relaxation_time));type12='hme';
nameA = strcat('shockPlusSmooth_','smoothpar',string(smooth_par),'_toldown',toldown,...
    '_tolup',tolup,'_t',string(time),'_Kn',relaxation_time);typeA='hmeAdaptive';

[x2,rho2,u2,T2,f32,f42,f52,f62,f72,f82,f92,f102,f112,f122,moments2] = readDataAdaptiveHME1D(name2,type2);
[x4,rho4,u4,T4,f34,f44,f54,f64,f74,f84,f94,f104,f114,f124,moments4] = readDataAdaptiveHME1D(name4,type4);
[x6,rho6,u6,T6,f36,f46,f56,f66,f76,f86,f96,f106,f116,f126,moments6] = readDataAdaptiveHME1D(name6,type6);
[x8,rho8,u8,T8,f38,f48,f58,f68,f78,f88,f98,f108,f118,f128,moments8] = readDataAdaptiveHME1D(name8,type8);
[x10,rho10,u10,T10,f310,f410,f510,f610,f710,f810,f910,f1010,f1110,f1210,moments10] = readDataAdaptiveHME1D(name10,type10);
[x12,rho12,u12,T12,f312,f412,f512,f612,f712,f812,f912,f1012,f1112,f1212,moments12] = readDataAdaptiveHME1D(name12,type12);
[xA,rhoA,uA,TA,f3A,f4A,f5A,f6A,f7A,f8A,f9A,f10A,f11A,f12A,momentsA] = readDataAdaptiveHME1D(nameA,typeA);

time2_arr = load(strcat('time',name2,'.csv'));
time4_arr = load(strcat('time',name4,'.csv'));
time6_arr = load(strcat('time',name6,'.csv'));
time8_arr = load(strcat('time',name8,'.csv'));
time10_arr = load(strcat('time',name10,'.csv'));
time12_arr = load(strcat('time',name12,'.csv'));
timeA_arr = load(strcat('time',nameA,'.csv'));

time2 = time2_arr;
time4 = time4_arr;
time6 = time6_arr;
time8 = time8_arr;
time10 = time10_arr;
time12 = time12_arr;
timeA = timeA_arr;

relTime2 = time2/time12;
relTime4 = time4/time12;
relTime6 = time6/time12;
relTime8 = time8/time12;
relTime10 = time10/time12;
relTime12 = time12/time12;
relTimeA = timeA/time12;

err2 = norm(rho2-rho12)/norm(rho12)+norm(u2-u12)/norm(u12)+norm(T2-T12)/norm(T12);
err4 = norm(rho4-rho12)/norm(rho12)+norm(u4-u12)/norm(u12)+norm(T4-T12)/norm(T12);
err6 = norm(rho6-rho12)/norm(rho12)+norm(u6-u12)/norm(u12)+norm(T6-T12)/norm(T12);
err8 = norm(rho8-rho12)/norm(rho12)+norm(u8-u12)/norm(u12)+norm(T8-T12)/norm(T12);
err10 = norm(rho10-rho12)/norm(rho12)+norm(u10-u12)/norm(u12)+norm(T10-T12)/norm(T12);
err12 = norm(rho12-rho12)/norm(rho12)+norm(u12-u12)/norm(u12)+norm(T12-T12)/norm(T12);
errA = norm(rhoA-rho12)/norm(rho12)+norm(uA-u12)/norm(u12)+norm(TA-T12)/norm(T12);

% err2 = norm((rho2-rho12)./rho12)+norm((u2-u12)./u12)+norm((T2-T12)./T12);
% err4 = norm((rho4-rho12)./rho12)+norm((u4-u12)./u12)+norm((T4-T12)./T12);
% err6 = norm((rho6-rho12)./rho12)+norm((u6-u12)./u12)+norm((T6-T12)./T12);
% err8 = norm((rho8-rho12)./rho12)+norm((u8-u12)./u12)+norm((T8-T12)./T12);
% err10 = norm((rho10-rho12)./rho12)+norm((u10-u12)./u12)+norm((T10-T12)./T12);
% err12 = norm((rho12-rho12)./rho12)+norm((u12-u12)./u12)+norm((T12-T12)./T12);
% errA = norm((rhoA-rho12)./rho12)+norm((uA-u12)./u12)+norm((TA-T12)./T12);

relTimes = [relTime2 relTime4 relTime6 relTime8 relTime10 relTime12 relTimeA]
errors = [err2 err4 err6 err8 err10 err12 errA]