%%Script to generate various plots

%% Initialize variables.
k = 32;

%%VARIANCE OF INLET ANGLE
% filename1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/45/postProcessing/';
% filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/45_50/postProcessing/';
% % filename3 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/45_75/postProcessing/';
% filename4 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu4percent/postProcessing/';
% filename5 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/46_10/postProcessing/';
% % filename6 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/46_25/postProcessing/';
% filename7 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/46_50/postProcessing/';
% filename8 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/VaryInletAngle/47/postProcessing/';

% plotCpProfile(filename1,90,'g--',k,'IA 45\circ');
% plotCpProfile(filename2,90,'c--',k,'IA 45.50\circ');
% % plotCpProfile(filename3,90,'y--',k,'IA 45.75');
% plotCpProfile(filename4,90,'b--',k,'IA 46\circ');
% plotCpProfile(filename5,90,'r--',k,'IA 46.10\circ');
% % plotCpProfile(filename6,90,'c--',k,'IA 46.25');
% plotCpProfile(filename7,90,'y--',k,'IA 46.50\circ');
% plotCpProfile(filename8,90,'m--',k,'IA 47\circ');

%% VARIANCE IN TU INTENSITY
%  filename1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu1percent/postProcessing/';
 filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu4percent/postProcessing/';
%  filename3 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu6percent/postProcessing/';
%  filename4 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu8percent/postProcessing/';

 filename1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/SAModelLimitedLinear/postProcessing/';
%  filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu4percent/postProcessing/';
%  filename3 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu6percent/postProcessing/';
%  filename4 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/TuIntensity/Tu8percent/postProcessing/';
 
%% VARIANCE IN SCHEMES
 filename9 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/Other/LimitedLinearSchemes/Tu4/postProcessing/';
% filename10 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/Other/LinearSchemes/Tu4/postProcessing/';

%%VARIANCE IN THE MODELS
filename1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/SAModelLimitedLinear/postProcessing/';
filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/kEpsilon/postProcessing/';
filename3 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/FromHPC/OtherModels/realizablekEpsilon/postProcessing/';

% plotCpProfile(filename2,90,'r--',k,'RANS SST','SST');
plotCpProfile(filename1,60,'b--',k,'RANS SA','SA');
plotCpProfile(filename9,30,'r--',k,'RANS  SST','SST');
plotCpProfile(filename67,7.994,'m--',k,'RANS Tu k-\epsilon','Epsilon');
plotCpProfile(filename3,88,'g--',k,'RANS realizable k-\epsilon','Epsilon');
% plotCpProfile(filename8,90,'c--',k,'RANS Tu 8%');

filename45 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/PostProcessData/Tu_4percent.csv';
DNS=csvread(filename45,0,0);
figure(k);
plot(DNS(:,1),DNS(:,2),'k--','LineWidth',2,'DisplayName','DNS Tu 4%');
legend('-DynamicLegend');
xlabel('x/C_{ax}');
ylabel('C_p');
set(gca,'fontsize', 18)
xlim([-0.02 1]);
k=k+1;
