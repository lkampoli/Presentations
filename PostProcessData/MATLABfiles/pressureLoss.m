function [] = pressureLoss(directory,time1,n,x)
figure(n)
mixedOutUp = getMixedOut(directory,time1,'up'); 
mixedOutDown = getMixedOut(directory,time1,'down'); 
PressureLoss = (mixedOutUp(2) - mixedOutDown(2))/(mixedOutUp(2)-mixedOutDown(1));
plot(x,PressureLoss,'k*','LineWidth',2); 
hold on;
% legend('-DynamicLegend');
end

