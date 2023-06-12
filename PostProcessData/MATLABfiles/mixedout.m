UpMeasure = csvread('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu4percent/DataCSV/UpstreamMeasurementPoint_Tu4.csv',1,0);
DownMeasure = csvread('/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu4percent/DataCSV/DownStreamMeasurementPoint_Tu4.csv',1,0);

mixedOutUpMeasure = mixedOutQuantity(UpMeasure(:,1),UpMeasure(:,4),UpMeasure(:,8),UpMeasure(:,9),UpMeasure(:,14));
mixedOutDownMeasure = mixedOutQuantity(DownMeasure(:,1),DownMeasure(:,4),DownMeasure(:,8),DownMeasure(:,9),DownMeasure(:,14));
PressureLoss = (mixedOutUpMeasure(2) - mixedOutDownMeasure(2))/(mixedOutUpMeasure(2)-mixedOutDownMeasure(1))
filename5 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/46_10/postProcessing/';
filename4 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu4percent/postProcessing/';

[mixedOutQuant] = getMixedOut(filename5,90,'up');


% plot((DownMeasure(:,8).^2+DownMeasure(:,9).^2).^0.5,DownMeasure(:,14));
% figure()
% plot(DownMeasure(:,8),DownMeasure(:,14));
filenameUp46_10_velocity = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/46_10/postProcessing/sets/90/downstream_U.xy';
filenameUp46_10_p_omega_k_rho = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/46_10/postProcessing/sets/90/downstream_p_omega_k_rho_T.xy';
UpMeasure46_10_velocity = textread(filenameUp46_10_velocity);
UpMeasure46_10_p_omega_k_rho = textread(filenameUp46_10_p_omega_k_rho);
