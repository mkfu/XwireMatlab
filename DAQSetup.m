%% DAQ Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs
%Analog
MotorOut.Channel = 'ao0'; MotorOut.Name = 'MotorOut';MotorOut.Range = [-10,10];
MotorOut.dev = 'Dev4';
%Digital
Pitot02.DChannel = 'port0/line0';   Pitot02Valve.DName = 'Pitot02Valve'; 
Pitot1.DChannel = 'port0/line1';    Pitot1Valve.DName = 'Pitot1Valve';  
Pitot5.DChannel = 'port0/line2';    Pitot5Valve.DName = 'Pitot5';   

ScaniHome.DChannel = 'port0/line5';  ScaniHome.DName = 'ScaniHome';
ScaniSkip.DChannel = 'port0/line6';  ScaniSkip.DName = 'ScaniSkip';   
ScaniPower.DChannel = 'port0/line7'; ScaniPower.DName = 'ScaniPower';       

%% Inputs
Temperature.Channel = 'ai0';        Temperature.cal = @(V) V*100+273.15; %kelvin
Temperature.Name = 'Temperature';   Temperature.Range = [-5,5];

Scanivalve.Channel = 'ai1';         Scanivalve.cal = @(V)   V.* 133.322; %1V/torr
Scanivalve.Name = 'Scanivalve';     Scanivalve.Range = [-5,5];

TunnelStatic.Channel = 'ai2';       TunnelStatic.cal = @(V) V.*4000/10.*6894.75729 ;
TunnelStatic.Name = 'TunnelStatic'; TunnelStatic.Range = [-10,10];

Dantec.Channel = 'ai3';             Dantec.cal =  @(P,V) polyval(P,V);
Dantec.Name = 'Dantec';             Dantec.Range = [-10,10];

Dantec.Channel = 'ai3';             Dantec.cal =  @(P,V) polyval(P,V);
Dantec.Name = 'Dantec';             Dantec.Range = [-10,10];

Pitot02.Channel = 'ai4';            Pitot02.cal = @(V) V*0.2/5*6894.75729; %0.2psi/5V
Pitot02.Name = 'Pitot02';           Pitot02.Range = [-5,5];

Pitot1.Channel = 'ai5';             Pitot1.cal = @(V) V*1/5*6894.75729; %1psi/5V
Pitot1.Name = 'Pitot1';             Pitot1.Range = [-5,5];

Pitot5.Channel = 'ai6';             Pitot5.cal = @(V) V*5/5*6894.75729; %1psi/5V
Pitot5.Name = 'Pitot5';             Pitot5.Range = [-5,5];

LimitSwitch.Channel = 'ai7';        LimitSwitch.cal = @(V) min(V)>1;
LimitSwitch.Name = 'LimitSwitch';   LimitSwitch.Range = [-10,10];

d = daq.getDevices; daqCal = daq.createSession('ni');