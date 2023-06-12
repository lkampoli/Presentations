function [] = plotCpProfile(directory,time1,mark,n,label2,model)
% Plots the Cp for given pressure profiles over the blade. 
figure(n)
delimiter = ' ';
startRow = 3;
endRow = 861;
formatSpec = '%s%s%s%s%[^\n\r]';
filename = strcat(directory,'surfaces/',num2str(time1),'/p_blade.raw');
fileID1 = fopen(filename,'r');
dataArray1 = textscan(fileID1, formatSpec, endRow, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
fclose(fileID1);
[x,y,z,p]= getDataFromRaw(dataArray1);

mixedOutUp = getMixedOut(directory,time1,'up',model); 
mixedOutDown = getMixedOut(directory,time1,'down',model); 

p1 = [p(3:193)' p(433:861)' p(194:432)'];
x1 = [x(3:193)' x(433:861)' x(194:432)'];
% cp=((p1 - Pout)/(PStagInlet-Pout)); 
cp=((p1 - mixedOutDown(1))/(mixedOutUp(2)-mixedOutDown(1))); 

plot(x1/0.8597,cp,mark,'LineWidth',2,'DisplayName',label2); 
hold on;
legend('-DynamicLegend');
end