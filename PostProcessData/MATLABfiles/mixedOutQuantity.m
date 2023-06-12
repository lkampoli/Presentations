function [mixedOutQuant] = mixedOutQuantity(rho,p,u,v,y)
%function returns mixed out quantities
gamma=1.4;
gam1=gamma-1.0;
mdot=0.0;
m1_mix=0.0;
m2_mix=0.0;
energ=0.0;
l = abs(y(1)-y(end));

for j =2:length(p)
    Delta = y(j-1)-y(j);
    mdot=mdot+((rho(j)*u(j)+rho(j-1)*u(j-1))*Delta/2);
    m1_mix=m1_mix+((rho(j)*u(j)*u(j)+p(j)+ rho(j-1)*u(j-1)*u(j-1)+p(j-1))*Delta/2);
    m2_mix=m2_mix+((rho(j)*u(j)*v(j)+ rho(j-1)*u(j-1)*v(j-1))*Delta/2);
    energ=energ+((u(j)*(gamma/gam1*p(j)+(u(j)*u(j)+ v(j)*v(j))*rho(j)*0.5) + u(j-1)*(gamma/gam1*p(j-1)+(u(j-1)*u(j-1)+ v(j-1)*v(j-1))*rho(j-1)*0.5))*Delta/2);    
end  

m0_store=mdot/l; 
m1_mix_out=m1_mix/l;  
m2_mix_out=m2_mix/l;  
ener_mix_out=energ/l; 

neg_p_half = m1_mix_out/m0_store*gamma/(gamma+1);
neg_q = 2./m0_store*gam1/(gamma+1.)*(m2_mix_out^2/m0_store*0.5-ener_mix_out);
u_mix_out = neg_p_half-sqrt(neg_p_half^2+neg_q);
rho_mix_out = m0_store/u_mix_out;
v_mix_out=m2_mix_out/rho_mix_out/u_mix_out;
p_mix_out=m1_mix_out-rho_mix_out*u_mix_out*u_mix_out;

% self.t_mix_out=self.p_mix_out*1.4*self.in_file.xm**2/self.rho_mix_out
t_mix_out = p_mix_out/(rho_mix_out*287);
c_mix_out=sqrt(u_mix_out^2+v_mix_out^2);
ma_mix_out=c_mix_out/sqrt(gamma*p_mix_out/rho_mix_out);
pt_mix_out = p_mix_out*(1+0.2*ma_mix_out^2)^(1.4/0.4);
tt_mix_out = t_mix_out*(1+0.2*ma_mix_out^2);
mu_mix_out = t_mix_out^1.5*(1.458*10^-6)/(t_mix_out+110.4);
angle_mix_out = atan(v_mix_out/u_mix_out)*45/atan(1);

mixedOutQuant = [p_mix_out pt_mix_out rho_mix_out u_mix_out v_mix_out c_mix_out ma_mix_out t_mix_out tt_mix_out mu_mix_out angle_mix_out];
end

