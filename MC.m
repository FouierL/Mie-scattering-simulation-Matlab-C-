% Monte Carlo Simulation
% output = MC(input)
% 
% for details about the data structure of input, see InitializeInput

% Bo Qiang
% last revised: 10/05/05

function [output] = MC(input)

setting.np = input.np;
setting.dz = input.dz;
setting.dr = input.dr;
setting.da = input.da;
setting.nz = input.nz;
setting.nr = input.nr;
setting.na = input.na;
setting.nl = input.nl;
setting.n_above = input.n_above;
setting.n_below = input.n_below;

n = [input.layers.n];
mua = [input.layers.mua];
mus = [input.layers.mus];
g = [input.layers.g];
d = [input.layers.z1];

rsp = 0;
rspp = libpointer('doublePtr',rsp);
rd_ra = zeros(input.nr,input.na);
rd_rap = libpointer('doublePtr',rd_ra);
ab_rz = zeros(input.nr,input.nz);
ab_rzp = libpointer('doublePtr',ab_rz);
tt_ra = zeros(input.nr,input.na);
tt_rap = libpointer('doublePtr',tt_ra);

loadlibrary('MC_dll.dll','MC_dll.h');

setting_Str=libstruct('SettingStruct',setting);

calllib('C:\Users\纳豆够巷\Desktop\生物医学光子学实验\MC_dll','MC',setting_Str,n,mua,mus,g,d,rspp,rd_rap,ab_rzp,tt_rap);

output.rsp = get(rspp,'value');
output.rd_ra = get(rd_rap,'value');
output.ab_rz = get(ab_rzp,'value');
output.tt_ra = get(tt_rap,'value');

output = ProcessOutput(input,output);

clear setting_Str
unloadlibrary MC_dll;