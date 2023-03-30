% Matlab Monte Carlo Simulation Tool (MMCST) v 1.2
%
% Sum and rescale the result
% [output] = ProcessOutput(input,output)
%
% will be called internally in MC
%
% Author: Qiang Bo
% qiang.2@osu.edu
% Last revised: 11/01/2004
% Last revised: 09/26/2005

function [output] = ProcessOutput(input,output)

output.rd_ra = output.rd_ra./input.np;
output.ab_rz = output.ab_rz./input.np;
ouptut.tt_ra = output.tt_ra./input.np;
output.rd_r = sum(output.rd_ra')';
output.rd_a = sum(output.rd_ra);
output.ab_r = sum(output.ab_rz')';
output.ab_z = sum(output.ab_rz);
output.tt_r = sum(output.tt_ra')';
output.tt_a = sum(output.tt_ra);
output.rd = sum(sum(output.rd_ra));
output.ab = sum(sum(output.ab_rz));
output.tt = sum(sum(output.tt_ra));

dr = input.dr;
da = input.da;
dz = input.dz;
nr = input.nr;
na = input.na;
nz = input.nz;

scale1 = 4.0*pi*pi*dr*sin(da/2)*dr;

for ir = 0:nr-1
    for ia = 0:na-1
        scale2 = 1.0/((ir+0.5)*sin(2.0*(ia+0.5)*da)*scale1);
        output.rd_ra(ir+1,ia+1) = output.rd_ra(ir+1,ia+1)*scale2;
        output.tt_ra(ir+1,ia+1) = output.tt_ra(ir+1,ia+1)*scale2;
    end
end


scale1 = 2.0*pi*dr*dr;

for ir = 0:nr-1
    scale2 = 1.0/((ir+0.5)*scale1);
    output.rd_r(ir+1) = output.rd_r(ir+1)*scale2;
    output.tt_r(ir+1) = output.tt_r(ir+1)*scale2;
end

scale1  = 4.0*pi*sin(da/2);

for ia = 0:na-1
    scale2 = 1.0/(sin((ia+0.5)*da)*scale1);
    output.rd_a(ia+1) = output.rd_a(ia+1)*scale2;
    output.tt_a(ia+1) = output.tt_a(ia+1)*scale2;
end

scale1 = 2.0*pi*dr*dr*dz;

for iz = 0:nz-1
    for ir = 0:nr-1
        output.ab_rz(ir+1,iz+1) = output.ab_rz(ir+1,iz+1)/((ir+0.5)*scale1);
    end
end

scale1 = 1.0/dz;
output.ab_z = output.ab_z*scale1;