data = importdata('matlab.mat');
plot(data(:,1),data(:,2));

%Initialising parameters
 Nb = 2;
 Cla = 5.73;
 R = 0.355;
 rpm = 1500;
 cd = 0.01;
 rho_air = 1.225;
 C = 0.032;
 
 %Calculating terms
 omega = rpm*pi/30;                              %Angular velocity
 sigma = (Nb*C)/(pi*R);                          %Advanced ratio
    
 dr = 1/100;
 r = 0.3:dr:1;                                    %with root cutout of 0.2
 theta = 2:1:15;   
 
 theta = theta.*(pi/180);                         %converitng to radians 
 lembda = zeros(length(theta),length(r));
 
%Loops to valculate AOA

 %initia lembda formula 
    for y = 1:length(theta)
        for x = 1:length(r)
            lembda(y,x) = sigma*(Cla/16)*(sqrt(1+(32*theta(y)*r(x))/(sigma*Cla))-1);
        end 
    end 
    
    
    %alpnha = theta - lembda/r
    
    for y = 1:length(theta)
        for x = 1:length(r)
            alpha(y,x) = theta(y) - lembda(y,x)/r(x);
        end 
    end 
%lookup table 
%interpolation
cl = interp1(data(:,1),data(:,2),(alpha.*(180/pi))); 

%again calculating lembda
for y = 1:length(theta)
    for x = 1:length(r)
        lembda1(y,x) = sqrt(((sigma*r(x))/8)*cl(y,x));
    end 
end 


%cost function 
cost1 = sum(abs((lembda1 - lembda)/length(lembda)),"all");
%it will not converge for all at one time 
lembda = lembda1;
cost0 = 0;
%Now we will try to converge Lembda by repeating the same steps again 
while cost1 > cost0
    cost0 = cost1;
    for y = 1:length(theta)
        for x = 1:length(r)
            alpha(y,x) = theta(y) - lembda1(y,x)/r(x);
            
        end 
    end
    
    alpha(isnan(alpha)) = 0; 
    
    cl = interp1(data(:,1),data(:,2),alpha.*(180/pi));
    
    cl(isnan(cl)) = 0;
    for y = 1:length(theta)
        for x = 1:length(r)
            lembda1(y,x) = sqrt(((sigma*r(x))/8)*cl(y,x));
        end 
    end
    lembda1(isnan(lembda1)) = 0;
     
    cost1 = sum(abs((lembda1 - lembda)/length(lembda)),"all");
    lembda = lembda1;
end
%Now we got lembda at every r, Now we will integrate over R 
%rectangular method integaration 
for y = 1:length(theta)
    for x = 1:length(r)
        dct(y,x) = 0.5*sigma*cl(y,x)*(r(x)^2)*dr;
    end 
end
ct = sum(dct,2);
T = rho_air*pi*(R^2)*(omega*R)^2.*ct;
plot(theta.*180/pi,T');
ylabel('Thrust(N)')
xlabel('Pitch Angle(deg)')