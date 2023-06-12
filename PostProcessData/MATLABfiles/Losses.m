%Calculate the pressure and wake losses, plotting the wake profile also. 
% n=90;
% filename1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/45/postProcessing/';
% filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/45_50/postProcessing/';
% filename3 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/45_75/postProcessing/';
% filename4 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu4percent/postProcessing/';
% filename5 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/46_10/postProcessing/';
% filename6 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/46_25/postProcessing/';
% filename7 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/46_50/postProcessing/';
% filename8 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/47/postProcessing/';
% 
% pressureLoss(filename1,90,n,45);
% pressureLoss(filename2,90,n,45.5);
% pressureLoss(filename3,90,n,45.75);
% pressureLoss(filename4,90,n,46);
% pressureLoss(filename5,90,n,46.1);
% pressureLoss(filename6,90,n,46.25);
% pressureLoss(filename7,90,n,46.5);
% pressureLoss(filename8,90,n,47);
n=98;
filename1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu1percent/postProcessing/';
 filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu4percent/postProcessing/';
 filename3 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu6percent/postProcessing/';
 filename4 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu8percent/postProcessing/';
pressureLoss(filename1,90,n,1);
pressureLoss(filename2,90,n,4);
pressureLoss(filename3,90,n,6);
pressureLoss(filename4,90,n,8);
 

fileDNSUp = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/DNSData/upstream.csv';
fileDNSDown = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/DNSData/downstream.csv';

DNSUp = csvread(fileDNSUp,1,0);
DNSDown = csvread(fileDNSDown,1,0);

mixedOutUpMeasure = mixedOutQuantity(DNSUp(:,1),DNSUp(:,5),DNSUp(:,2),DNSUp(:,3),DNSUp(:,9));
mixedOutDownMeasure = mixedOutQuantity(DNSDown(:,1),DNSDown(:,5),DNSDown(:,2),DNSDown(:,3),DNSDown(:,9));
PressureLoss = (mixedOutUpMeasure(2) - mixedOutDownMeasure(2))/(mixedOutUpMeasure(2)-mixedOutDownMeasure(1))
plot(4,PressureLoss,'r*','LineWidth',2,'DisplayName','DNS Tu 4%');
ylabel('\omega_M');
xlabel('Tu Intensity');
set(gca,'fontsize', 18)
% legend('RANS')
legend('-DynamicLegend');
% xlim([-0.02 1]);

