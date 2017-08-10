%% DAQ Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs
%Analog
MotorOut.Channel = 'ao0'; MotorOut.Name = 'MotorOut';MotorOut.Range = [-10,10];
MotorOut.dev = 'Dev4';
%Digital
Pitot02.DChannel = 'port0/line0';   Pitot02Valve.DName = 'Pitot02Valve';    
Pitot02.Ddev = 'Dev4';

Pitot1.DChannel = 'port0/line1';    Pitot1Valve.DName = 'Pitot1Valve';      
Pitot1.Ddev = 'Dev4';

Pitot5.DChannel = 'port0/line2';    Pitot5Valve.DName = 'Pitot5';           
Pitot5.Ddev = 'Dev4';

ScaniHome.DChannel = 'port0/line5';  ScaniHome.DName = 'ScaniHome';         
ScaniHome.Ddev = 'Dev4';

ScaniSkip.DChannel = 'port0/line6';  ScaniSkip.DName = 'ScaniSkip';         
ScaniSkip.Ddev = 'Dev4';

ScaniPower.DChannel = 'port0/line7'; ScaniPower.DName = 'ScaniPower';       
ScaniPower.Ddev = 'Dev4';     

%% Inputs
Temperature.Channel = 'ai0';        Temperature.cal = @(V) V*100+273.15; %kelvin
Temperature.Name = 'Temperature';   Temperature.Range = [-5,5];
Temperature.dev = 'Dev4';

Scanivalve.Channel = 'ai1';         Scanivalve.cal = @(V)   V.* 133.322; %1V/torr
Scanivalve.Name = 'Scanivalve';     Scanivalve.Range = [-5,5];
Scanivalve.dev = 'Dev4';

TunnelStatic.Channel = 'ai2';       TunnelStatic.cal = @(V) V.*4000/10.*6894.75729 ;
TunnelStatic.Name = 'TunnelStatic'; TunnelStatic.Range = [-10,10];
TunnelStatic.dev = 'Dev4';

hw1.Channel = 'ai3';                hw1.cal =  @(P,V) polyval(P,V);
hw1.Name = 'Dantec';                hw1.Range = [-10,10];
hw1.dev = 'Dev4';

hw2.Channel = 'ai7';                hw2.cal =  @(P,V) polyval(P,V);
hw2.Name = 'Dantec';                hw2.Range = [-10,10];
hw2.dev = 'Dev4';

Pitot02.Channel = 'ai4';            Pitot02.cal = @(V) V*0.2/5*6894.75729; %0.2psi/5V
Pitot02.Name = 'Pitot02';           Pitot02.Range = [-5,5];
Pitot02.dev = 'Dev4';

Pitot1.Channel = 'ai5';             Pitot1.cal = @(V) V*1/5*6894.75729; %1psi/5V
Pitot1.Name = 'Pitot1';             Pitot1.Range = [-5,5];
Pitot1.dev = 'Dev4';

Pitot5.Channel = 'ai6';             Pitot5.cal = @(V) V*5/5*6894.75729; %1psi/5V
Pitot5.Name = 'Pitot5';             Pitot5.Range = [-5,5];
Pitot5.dev = 'Dev4';

% LimitSwitch.Channel = 'ai7';        LimitSwitch.cal = @(V) min(V)>1;
% LimitSwitch.Name = 'LimitSwitch';   LimitSwitch.Range = [-10,10];
% LimitSwitch.dev = 'Dev4';

d = daq.getDevices; daqCal = daq.createSession('ni');