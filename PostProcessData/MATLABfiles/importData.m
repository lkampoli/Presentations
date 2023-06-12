%% Import data from text file.
% Script for importing data from the following text file:

%% Initialize variables.
filename1 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/TURB_OFF/postProcessing/surfaces/3.75/p_blade.raw';
filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/TURB_OFF/postProcessing/surfaces/5.9/p_blade.raw';
filename3 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/TURB_OFF/postProcessing/surfaces/13.9/p_blade.raw';

% filename2 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/TurbModelOnTu_4per/postProcessing/surfaces/4.15/p_blade.raw';
% filename3 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/Tu4percent/postProcessing/surfaces/60/p_blade.raw';
% filename4 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/Tu8percent/postProcessing/surfaces/30/p_blade.raw';
% filename5 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/ToHPC/Tu6percent/postProcessing/surfaces/30/p_blade.raw';

% filename3 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/DealWithUProblem/PStag7770Pout7100/postProcessing/surfaces/8/p_blade.raw';
% filename4 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/DealWithUProblem/PStag7870/postProcessing/surfaces/8/p_blade.raw';
% filename5 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/DealWithUProblem/PStag7920/postProcessing/surfaces/8/p_blade.raw';
% filename6 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/DealWithUProblem/PStag7670/postProcessing/surfaces/8/p_blade.raw';
m=5;
k=57;
delimiter = ' ';
startRow = 3;
endRow = 861;
StagP = [7770, 7770, 7770, 7770, 7770];
Pstat = [6950, 6950, 6950, 6950, 6950];
formatSpec = '%s%s%s%s%[^\n\r]';
 fileID1 = fopen(filename1,'r');
 fileID2 = fopen(filename2,'r');
 fileID3 = fopen(filename3,'r');
%  fileID4 = fopen(filename4,'r');
%  fileID5 = fopen(filename5,'r');

% for j =1:2  
%     strFile = sprintf('filename%d', j);
%     fileID = fopen(filename1,'r');
% 

    dataArray1 = textscan(fileID1, formatSpec, endRow, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
    dataArray2 = textscan(fileID2, formatSpec, endRow, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
    dataArray3 = textscan(fileID3, formatSpec, endRow, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
%     dataArray4 = textscan(fileID4, formatSpec, endRow, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
%     dataArray5 = textscan(fileID5, formatSpec, endRow, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
%     dataArray6 = textscan(fileID6, formatSpec, endRow, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false);


    %% Close the text file.
    fclose(fileID1);
    fclose(fileID2);
    fclose(fileID3);
%     fclose(fileID4);
%     fclose(fileID5);
    [x1,y1,z1,p1]= getDataFromRaw(dataArray1);
    [x2,y2,z2,p2]= getDataFromRaw(dataArray2);
    [x3,y3,z3,p3]= getDataFromRaw(dataArray3);
%     [x4,y4,z4,p4]= getDataFromRaw(dataArray4);
%     [x5,y5,z5,p5]= getDataFromRaw(dataArray5);

    %% Clear temporary variables
    clearvars filename delimiter endRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;
%     label2 = 'PStag = ' num2str(PStag(j))PStat = 
    plotCpProfile(x1,p1,StagP(1),Pstat(1),'g--',k,'Laminar 3750');
    plotCpProfile(x2,p2,StagP(2),Pstat(2),'r--',k,'Laminar 5900');
    plotCpProfile(x3,p3,StagP(3),Pstat(3),'y--',k,'Laminar 13900');
%     plotCpProfile(x4,p4,StagP(4),Pstat(4),'b--',k,'Tu=6%');
%     plotCpProfile(x5,p5,StagP(5),Pstat(5),'g--',k,'Tu=8%');


    filename45 = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/PostProcessData/Tu_4percent.csv';
    [xDNS, pDNS]=ImportCSV(filename45);
    figure(k);
    plot(xDNS,pDNS,'k--','LineWidth',2,'DisplayName','DNS Tu=4%');
    legend('-DynamicLegend');
    xlabel('x/C');
    ylabel('C_p');
    set(gca,'fontsize', 18)
    
% end

%recalculate these values and make sure that I have the right conditions. 
