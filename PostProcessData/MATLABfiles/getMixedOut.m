function [mixedOutQuant] = getMixedOut(directory,time1,location,model)
fileVelo ='a';
fileOther = 'b';
if strcmp(model,'SA')
    if strcmp(location,'in')
        fileVelo = strcat(directory,'sets/',num2str(time1),'/inlet_U.xy');
        fileOther = strcat(directory,'sets/',num2str(time1),'/inlet_p_nuTilda_rho_T.xy');
    elseif strcmp(location,'out')
        fileVelo = strcat(directory,'sets/',num2str(time1),'/outlet_U.xy');
        fileOther = strcat(directory,'sets/',num2str(time1),'/outlet_p_nuTilda_rho_T.xy');
    elseif strcmp(location,'up')
        fileVelo = strcat(directory,'sets/',num2str(time1),'/upstream_U.xy');
        fileOther = strcat(directory,'sets/',num2str(time1),'/upstream_p_nuTilda_rho_T.xy');
    %     velocity = textscan(fileVelo,'%n');
    %     other = textscan(fileOther,'%n');
    elseif strcmp(location,'down')
        fileVelo = strcat(directory,'sets/',num2str(time1),'/downstream_U.xy');
        fileOther = strcat(directory,'sets/',num2str(time1),'/downstream_p_nuTilda_rho_T.xy');
    else 
        'Please enter correct location'
    end
elseif strcmp(model,'SST')
    if strcmp(location,'in')
        fileVelo = strcat(directory,'sets/',num2str(time1),'/inlet_U.xy');
        fileOther = strcat(directory,'sets/',num2str(time1),'/inlet_p_omega_k_rho_T.xy');
    elseif strcmp(location,'out')
        fileVelo = strcat(directory,'sets/',num2str(time1),'/outlet_U.xy');
        fileOther = strcat(directory,'sets/',num2str(time1),'/outlet_p_omega_k_rho_T.xy');
    elseif strcmp(location,'up')
        fileVelo = strcat(directory,'sets/',num2str(time1),'/upstream_U.xy');
        fileOther = strcat(directory,'sets/',num2str(time1),'/upstream_p_omega_k_rho_T.xy');
    %     velocity = textscan(fileVelo,'%n');
    %     other = textscan(fileOther,'%n');
    elseif strcmp(location,'down')
        fileVelo = strcat(directory,'sets/',num2str(time1),'/downstream_U.xy');
        fileOther = strcat(directory,'sets/',num2str(time1),'/downstream_p_omega_k_rho_T.xy');
    else 
        'Please enter correct location'
    end
else  
    if strcmp(location,'in')
        fileVelo = strcat(directory,'sets/',num2str(time1),'/inlet_U.xy');
        fileOther = strcat(directory,'sets/',num2str(time1),'/inlet_p_epsilon_k_rho_T.xy');
    elseif strcmp(location,'out')
        fileVelo = strcat(directory,'sets/',num2str(time1),'/outlet_U.xy');
        fileOther = strcat(directory,'sets/',num2str(time1),'/outlet_p_epsilon_k_rho_T.xy');
    elseif strcmp(location,'up')
        fileVelo = strcat(directory,'sets/',num2str(time1),'/upstream_U.xy');
        fileOther = strcat(directory,'sets/',num2str(time1),'/upstream_p_epsilon_k_rho_T.xy');
    %     velocity = textscan(fileVelo,'%n');
    %     other = textscan(fileOther,'%n');
    elseif strcmp(location,'down')
        fileVelo = strcat(directory,'sets/',num2str(time1),'/downstream_U.xy');
        fileOther = strcat(directory,'sets/',num2str(time1),'/downstream_p_epsilon_k_rho_T.xy');
    else 
        'Please enter correct location'
    end
end
    
velocity = textread(fileVelo);
other    = textread(fileOther);

if strcmp(model,'SA')     
    mixedOutQuant = mixedOutQuantity(other(:,4),other(:,2),velocity(:,2),velocity(:,3),velocity(:,1));
else
    mixedOutQuant = mixedOutQuantity(other(:,5),other(:,2),velocity(:,2),velocity(:,3),velocity(:,1));
end
end

