%Wall shear stress 
delimiter = ' ';
startRow = 3;
endRow = 861;
formatSpec = '%s%s%s%s%[^\n\r]';
% filenameT = strcat(directory,'surfaces/',num2str(time1),'/p_blade.raw');
filenameT = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/Tu4percent/postProcessing/surfaces/90/T_blade.raw';
filenameWallGradU = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/Tu4percent/postProcessing/surfaces/90/wallGradU_blade.raw';
filenameWallShearStress = '/home/harshal/OpenFOAM/harshal-2.4.0/run/RAS/T106A/LPT_DNSGrid/Steady/FunctioningCases/Tu4percent/postProcessing/surfaces/90/wallShearStress_blade.raw';

fileID1 = fopen(filenameT,'r');
fileID2 = fopen(filenameWallGradU,'r');
fileID3 = fopen(filenameWallShearStress,'r');

dataArray1 = textscan(fileID1, formatSpec, endRow, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
fclose(fileID1);
dataArray2 = textscan(fileID2, formatSpec, endRow, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
fclose(fileID2);
dataArray3 = textscan(fileID3, formatSpec, endRow, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
fclose(fileID3);


[x,y,z,T]= getDataFromRaw(dataArray1);
[x1,y1,z1,wallGradU, wallGradV,wallGradW]= getDataFromRaw6(filenameWallGradU,3,861);
[x3,y3,z3,wallShearX, wallShearY,wallShearZ]= getDataFromRaw6(filenameWallShearStress,3,861);

mu = (T.^1.5*1.458*10^-6)./(T+110.4);
TauWall = ((wallGradU.^2+wallGradV.^2).^0.5).*mu(3:861);
figure();
TauWall1 = [TauWall(1:191)' TauWall(431:859)' TauWall(192:430)'];
TauWall2 = 2.375*10^-4 *0.08 *((wallShearX.^2+wallShearZ.^2).^0.5).*mu(3:861);
x2 = [x(3:193)' x(433:861)' x(194:432)'];

figure(5)
plot(x2/0.86,TauWall1,'r--');
hold on;
% plot(x2/0.86,TauWall2,'b--');
xlabel('x/C_{ax}');
ylabel('\tau_W');
