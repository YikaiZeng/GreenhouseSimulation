%% Initialization
% clear all;
clc;
close all;
%% Simulation Mode
SMode = 0;
% 1 - 24h
% 0 - continous
if SMode == 0
    DateStartM  = 7;
    DateStartD  = 10;
    DateEndM    = 7;
    DateEndD    = 25;
end
%% Construction Parameter
length      = 400;             % length of the greenhouse
width       = 250;              % width of the greenhouse
h_gutter    = 6;               % gutter height
A_floor     = length * width;  % floor square
theta       = 26/(180 * pi);   % angle of the roof
h_ridge     = h_gutter + 0.5 * width * tan(theta);
V_air       = A_floor * h_gutter;
V_airAS     = A_floor * 0.5 * (h_ridge - h_gutter);
%V_air      = 6840;
A_wall      = 2 * (length * h_gutter + width * h_gutter );
%% Cover
A_cov       = A_floor / cos(theta);% Surface of cover
%A_cov = 1842;                       % Other shape
lambda_cov  = 1.16;                 % Thermal heat conductivity of the cover W/m
d_cov       = 5e-3;                 % Thickness of the cover layer
%% Wall
k_glass = 0.88;
k_air   = 0.0274;
Lo = 3e-3;
Li = 6e-3;
ki = k_glass;
ko = k_air;
U  = 1/(2*Lo/ki + Li/ko);
%% Day of year of simulation for 24h
if SMode == 1
    month = 12;  % The month of simulation
    date  = 13; % The date of simulation
    tnow = DOY(month,date);
end
%% Load data 1
% Year Info for ground temperature
% this is supposed to be used for the ground temperature
% however the equation doesn't give a reliable result
%% Load data 2
%Radiation & temperature & Relative Humid
%load('date.mat');
if SMode == 1
date = csvread("date(20180124).csv",1,0);
%date = csvread("date.csv",1,0);
    Temp            = date(:,5);
    I               = date(:,3);
    I_t             = date(:,4);
    I_p             = date(:,2);
    RH              = date(:,6);
    Vwind           = date(:,7)/3.6;
else
    load("2018Yixing.mat");
    ts = DOY(DateStartM,DateStartD);
    te = DOY(DateEndM,  DateEndD);
    td = Overalldata(:,1);
    index1 = find(td <= te);
    index2 = find(td >= ts);
    index  = intersect(index1,index2);
    index  = [index;index(end)+1];
    Temp   = Overalldata(index,5);
    I      = Overalldata(index,7);
    RH     = Overalldata(index,6);
    Vwind  = Overalldata(index,8)/3.6; % into ,m/s
    I_t    = Overalldata(index(end):index(end)+23,7);
    I_p    = Overalldata((index(1) - 24):(index(1) - 1),7);
end
%% Weather data Humid air 
% Saturation pressure
MM_water = 18e-3;
R = 8.314;
clear length;
VP_SatO = zeros(length(Temp),1);
for i = 1:length(Temp)
    VP_SatO(i) = satVP(Temp(i));
end
VP = (RH/100).* VP_SatO;
C_H2O_out = VP * MM_water ./ (R .* (Temp + 273.15));
%% Energy Demand - Winter
T_set    = 15; % Maintain this duiring the night
Mode     = 0;  % Cooling : 0; Heating : 1
AC_ON    = 1;  % Turning off - 0 ; Turing on - 1;
RH_setAC = 85;
%% Roof
T_set_roof  = 30;
RH_set_roof = 85;
%% Weather data - Sky
if SMode == 1
    T_sky = SkyTemp(I_p,I_t,I,Temp,tnow);
else
    T_sky = zeros(length(Temp),1);
    I_sky = [I_p;I;I_t];
    for i = 1:floor(length(Temp)/24)
        Ip = I_sky((24*(i-1) + 1):(24*(i-1) + 24));
        Ic = I_sky(24*i+1:24*i+24);
        It = I_sky((24*(i+1) + 1):(24*(i+1) + 24));
        T  = Temp((24*(i-1) + 1):(24*(i-1) + 25));
        tnow = ts + (i-1);
        T_sky((24*(i-1) + 1):(24*(i-1) + 25)) = SkyTemp(Ip,It,Ic,T,tnow);
    end
end
%% CO2 concentration 
C_CO2 = 400; % ppm
%% Simulation Boundary
if AC_ON == 1 
T_start_cov   = T_set;
T_start_air   = T_set;
T_start_floor = T_set;
T_start_can   = T_set;
VP_sat_start  = satVP(T_set);
VP_start      = (RH_setAC/100) *  VP_sat_start;
C_H2O_start   = (VP_start * MM_water) / (R * (T_set + 273.15));
else % Natural Start
T_start_cov   = Temp(1);
T_start_air   = Temp(1);
T_start_floor = Temp(1);
T_start_can   = Temp(1);
C_H2O_start   = C_H2O_out(1); % For the start in the greenhouse 
end
LAI = 2; %Leaf Area Index
%% GND temperature
dfloor = 0.05;
dsoil1 = 0.08;
dsoil2 = 0.12;
dsoil3 = 0.2;
dsoil4 = 0.36;
dsoil5 = 0.72;
if AC_ON == 1
    Tsoil1 = GNDT(dfloor + 0.5 * dsoil1,T_set);
    Tsoil2 = GNDT(dfloor + dsoil1 + 0.5 * dsoil2,T_set);
    Tsoil3 = GNDT(dfloor + dsoil1 + dsoil2 + 0.5 * dsoil3,T_set);
    Tsoil4 = GNDT(dfloor + dsoil1 + dsoil2 + dsoil3 + 0.5 * dsoil4,T_set);
    Tsoil5 = GNDT(dfloor + dsoil1 + dsoil2 + dsoil3 + dsoil4 + 0.5 * dsoil5,T_set);
    %TGND   = GNDT(dfloor + dsoil1 + dsoil2 + dsoil3 + dsoil4 + dsoil5 ,T_set);
else
    Tsoil1 = GNDT(dfloor + 0.5 * dsoil1,Temp(1));
    Tsoil2 = GNDT(dfloor + dsoil1 + 0.5 * dsoil2,Temp(1));
    Tsoil3 = GNDT(dfloor + dsoil1 + dsoil2 + 0.5 * dsoil3,Temp(1));
    Tsoil4 = GNDT(dfloor + dsoil1 + dsoil2 + dsoil3 + 0.5 * dsoil4,Temp(1));
    Tsoil5 = GNDT(dfloor + dsoil1 + dsoil2 + dsoil3 + dsoil4 + 0.5 * dsoil5,Temp(1));
    %TGND   = GNDT(dfloor + dsoil1 + dsoil2 + dsoil3 + dsoil4 + dsoil5 ,Temp(1));
end
TGND = GNDT(dfloor + dsoil1 + dsoil2 + dsoil3 + dsoil4 + dsoil5 ,mean(Temp));
%% Stack the time vector
if SMode == 1
    Hour = (0:24)';
else
    Hour = (0:24*(te-ts+1))';
end
Hour = Hour * 3600;
T_out       = [Hour,Temp ];
T_sky       = [Hour,T_sky];
RH_out      = [Hour,RH];
C_H2O_out   = [Hour,C_H2O_out];
Vwind       = [Hour,Vwind];
I           = [Hour,I];
%% Apply to model change
simIn = Simulink.SimulationInput('ModelMarco_Full');
simIn = simIn.setModelParameter('StartTime','0','StopTime',num2str(Hour(end)));
simIn = setBlockParameter(simIn,'ModelMarco_Full/RHset_AC','Value',num2str(RH_setAC));
simIn = setBlockParameter(simIn,'ModelMarco_Full/Tset_vent','Value',num2str(T_set_roof));
simIn = setBlockParameter(simIn,'ModelMarco_Full/VairAS','Value',num2str(V_airAS));
simIn = setBlockParameter(simIn,'ModelMarco_Full/RHset_vent','Value',num2str(RH_set_roof));
% simIn = setBlockParameter(simIn,'ModelMarco_Full/Mode','Value',num2str(Mode));
simIn = setBlockParameter(simIn,'ModelMarco_Full/AC_ON','Value',num2str(AC_ON));
simIn = setBlockParameter(simIn,'ModelMarco_Full/Tset','Value',num2str(T_set));
simIn = setBlockParameter(simIn,'ModelMarco_Full/A_cov','Value',num2str(A_cov));
simIn = setBlockParameter(simIn,'ModelMarco_Full/A_wall','Value',num2str(A_wall));
simIn = setBlockParameter(simIn,'ModelMarco_Full/U','Value',num2str(U));
simIn = setBlockParameter(simIn,'ModelMarco_Full/A_floor','Value',num2str(A_floor));
simIn = setBlockParameter(simIn,'ModelMarco_Full/C_CO2','Value',num2str(C_CO2));
% simIn = setBlockParameter(simIn,'ModelMarco/Theta','Value',num2str(theta));
simIn = setBlockParameter(simIn,'ModelMarco_Full/T_start_cov','Value',num2str(T_start_cov));
simIn = setBlockParameter(simIn,'ModelMarco_Full/T_start_air','Value',num2str(T_start_air));
simIn = setBlockParameter(simIn,'ModelMarco_Full/T_start_can','Value',num2str(T_start_can));
simIn = setBlockParameter(simIn,'ModelMarco_Full/T_start_floor','Value',num2str(T_start_floor));
simIn = setBlockParameter(simIn,'ModelMarco_Full/LAI','Value',num2str(LAI));
simIn = setBlockParameter(simIn,'ModelMarco_Full/TGND','Value',num2str(TGND));
simIn = setBlockParameter(simIn,'ModelMarco_Full/T_start_Soil1','Value',num2str(Tsoil1));
simIn = setBlockParameter(simIn,'ModelMarco_Full/T_start_Soil2','Value',num2str(Tsoil2));
simIn = setBlockParameter(simIn,'ModelMarco_Full/T_start_Soil3','Value',num2str(Tsoil3));
simIn = setBlockParameter(simIn,'ModelMarco_Full/T_start_Soil4','Value',num2str(Tsoil4));
simIn = setBlockParameter(simIn,'ModelMarco_Full/T_start_Soil5','Value',num2str(Tsoil5));
%simIn = setBlockParameter(simIn,'ModelMarco/Tmean','Value',num2str(Tmean));
%simIn = setBlockParameter(simIn,'ModelMarco/Tamp','Value',num2str(Tamp));
%simIn = setBlockParameter(simIn,'ModelMarco/tnow','Value',num2str(tnow));
%simIn = setBlockParameter(simIn,'ModelMarco/tshift','Value',num2str(tshift));
simIn = setBlockParameter(simIn,'ModelMarco_Full/Vair','Value',num2str(V_air));
simIn = setBlockParameter(simIn,'ModelMarco_Full/CH2O_start','Value',num2str(C_H2O_start));
simIn.applyToModel;
%% Hint
disp("Simulation Start");
simOut = sim(simIn);
clc;
disp("Simulation Done");
%% Visulization
if SMode == 1
    figure(1);
    plot(T_out(:,1)/3600,T_out(:,2),'y');
    hold on;
    yyaxis left;
    plot(simOut.T_airin.time/3600,simOut.T_airin.data,'c');
    plot(simOut.T_can.time/3600,simOut.T_can.data,'g');
    xlim([0 24]);
    xticks(0:24);
    grid on;
    ylabel("Temperature /¡æ");
    xlabel("Hour of the day");
    yyaxis right
    plot(simOut.RH_in.time/3600,simOut.RH_in.data,'b');
    ylabel("Relative Humidity (%)");
    legend("Outsde Temperature","Inside Air Temperature","Canopy Temperature","Relative Humidity - Inside");
    title("24h Simulation in Greenhouse");
    hold off;
    %% Figure 2
    figure(2);
    yyaxis left;
    plot(simOut.P_in.time/3600,simOut.P_in.data/1000);
    hold on;
    ylabel("Power Input (kW)");
    xlabel("Hour of the day");
    xlim([0 24]);
    xticks(0:24);
    grid on;
    yyaxis right;
    plot(simOut.Q_in.time/3600,simOut.Q_in.data/(3600000000));
    ylabel("Accumulated Heat (MWh)");
    title("24h Simulation in Greenhouse - Energy");
    hold off;
    %% figure(3)
    figure(3);
    subplot(2,1,1);
    title("Local weather data 24h");
    yyaxis left;
    plot(T_out(:,1)/3600,RH);
    ylabel("Relative Humidity(%)");
    hold on;
    yyaxis right;
    plot(T_out(:,1)/3600,T_out(:,2));
    ylabel("Temperature(¡æ)");
    legend("Relative Humidity","Outside Temperature");
    grid on;
    xlim([0 24]);
    xticks(0:24);
    subplot(2,1,2);
    yyaxis left;
    plot(I(:,1)/3600,I(:,2));
    title("Local radiation & wind speed data 24h");
    ylabel("Global Radiation(W/m^2)");
    grid on;
    xlim([0 24]);
    xticks(0:24);
    yyaxis right;
    plot(Vwind(:,1)/3600,Vwind(:,2));
    ylabel("Wind speed(m/s)");
    xlabel("Hour of the day");
    hold off;
    %% Demand - winter
    if AC_ON == 1
        figure(4);
        yyaxis left;
        plot(simOut.Q_heating.time/3600,simOut.Q_heating.data/1000);
        ylabel("Heating Power / kW");
        grid on;
        yyaxis right;
        plot(simOut.Heating_Energy.time/3600,simOut.Heating_Energy.data/3600000);
        ylabel("Energy Used / kWh");
        xlabel("Hour of the day");
        xlim([0 24]);
        xticks(0:24);
        title("Energy Demand Simulation");
        disp("Energy Demand per day for " + num2str(T_set) + "¡æ lower bound = " + num2str(simOut.Heating_Energy.data(end)/3600000) + "kWh");
        disp("Water to remove per day = " + num2str(simOut.WaterAmount.data(end)) + "kg");
    end
    disp("24 Mean Temperature of canopy is " + num2str(mean(simOut.T_can)) + "¡æ");
    %% Thermal Screen
    figure(5);
    plot(simOut.T_airin.time/3600,simOut.T_airin.data);
    hold on;
    plot(simOut.T_airAS.time/3600,simOut.T_airAS.data);
    xlim([0 24]);
    xticks(0:24);
    ylabel("Temperature(¡æ)");
    xlabel("Hour of the day");
    legend("Inside - Below Screen","Inside - Above Screen");
    title("Thermal Screen Damping Effect");
    grid on;
    hold off;
else
    figure(1);
    subplot(2,1,1);
    plot(simOut.T_airin.time/3600,simOut.T_airin.data,'c');
    hold on;
    plot(simOut.T_can.time/3600,simOut.T_can.data,'g');
    plot(T_out(:,1)/3600,T_out(:,2));
    legend("Inside - Air Temperature","Canopy Temperature","Outside Temperature");
    ylabel("Temperature(¡æ)");
    xlabel("Hours Since "+num2str(DateStartM)+"-"+num2str(DateStartD));
    xlim([0 simOut.T_can.time(end)/3600]);
    grid on;
    subplot(2,1,2);
    plot(simOut.RH_in.time/3600,simOut.RH_in.data,'b');
    title("Inside Greenhouse");
    ylabel("Relative Humidity(%)");
    xlabel("Hours Since "+num2str(DateStartM)+"-"+num2str(DateStartD));
    xlim([0 simOut.T_can.time(end)/3600]);
    grid on;
    hold off;
    figure(2);
    plot(simOut.Q_heating.time/3600,simOut.Q_heating.data/1000);
    ylabel("Air Condition Power Output / kWh");
    xlabel("Hours Since "+num2str(DateStartM)+"-"+num2str(DateStartD));
    xlim([0 simOut.T_can.time(end)/3600]);
    disp("The Simulation Period is 2018 "+num2str(DateStartM)+"-"+num2str(DateStartD)+" to "+num2str(DateEndM)+"-"+num2str(DateEndD));
    disp("Energy Demand for the whole phase is " + num2str(T_set) + "¡æ lower bound = " + num2str(simOut.Heating_Energy.data(end)/3600000) + "kWh");
    disp("Energy Per day =  "+ num2str(simOut.Heating_Energy.data(end)/3600000/(te-ts+1))+ "kWh");
    disp("Water to remove per day = " + num2str(simOut.WaterAmount.data(end)/(te-ts+1)) + "kg");
end
%% for Extraterrestrial radiation
function [Ea] = ExtraRadiation(day)
    day = rem(day,365);
    E_sc = 1367;
    b = 2 * pi * day / 365;
    R = 1.00011 + 0.034221 * cos(b) + 0.00128 * sin(b) + 0.000719 * cos(2*b) + 0.000077 * sin(2*b);
    Ea = E_sc * R;
end
%% GND temperature
function [Tgnd] = GNDT(z,Tsurf)
% cp_soil = 1050;
% lambda_soil = 0.85;
% rho_soil = 1400;
% D = lambda_soil / (rho_soil * cp_soil);
% Tgnd = Tmean - Tamp * exp(-z * (pi * D / 365)^0.5) * cos((2 * pi / 365) * (tnow - tshift - (z/2) * (365*D/pi)^0.5));
Tgnd = Tsurf + z * 0.03;
end
%% VP_sat
function VP_sat = satVP(Temp)
    VP_h = 610.780 * exp(17.08085 * Temp/(234.175 + Temp));
    VP_l = 610.714 * exp(22.44294 * Temp/(272.440 + Temp));
    SF_H = 1/(1 + exp(-50 * Temp));
    SF_L = 1/(1 + exp(50 * Temp));
    VP_sat = VP_h * SF_H + VP_l * SF_L;
end
%%
function t = DOY(M,D)
Monthday    = [31,28,31,30,31,30,31,31,30,31,30,31]; 
if M == 1
    t = D;
else
    t = sum(Monthday(1:(M-1))) + D;
end
end
%% Sky Temperature
function T_sky = SkyTemp(Ip,It,I,Temp,tnow)
    CF_p = mean(Ip)/ ExtraRadiation(tnow - 1);
    CF_t = mean(It)/ ExtraRadiation(tnow + 1);
    CF   = mean(I) / ExtraRadiation(tnow);
    day_index = find(I > 0) ;
    day_start = day_index(1);
    day_end   = day_index(end);
    CF_time   = zeros(25,1);
    CF_time(day_start:day_end) = CF;
    grad_s    = (CF - CF_p)/(day_start -1);
    grad_e    = (CF_t - CF)/(24 - day_end) ;
    for i = 1 : (day_start - 1)
        CF_time(i) = CF_p + (i - 1) * grad_s;
    end
    for i = (day_end + 1) : 25
        CF_time(i) = CF + ( i - day_end - 1) * grad_e;
    end
    % Sky temp
    eplison_skyclear = 1;
    sigma = 5.67e-8;
    T_sky = zeros(24,1);
    for i = 1 : 25
        T_sky(i) = ((1 - CF_time(i)) * eplison_skyclear * ( Temp(i) + 273.15 )^4 + CF_time(i) * ((Temp(i) + 273.15)^4 - 9/sigma))^0.25 - 273.15;
    end
end