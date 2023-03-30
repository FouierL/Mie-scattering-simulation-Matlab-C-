% Matlab Monte Carlo Simulation Tool (MMCST)
% Initialize the input structure for program MC (v1.1)
%
% input = InitializeInput(inputfilename)
% inputfilename: the file contains the input parameters for Monte Carlo
% Simulation
% input: the data structure which can be used as the input in program MC
%        details about the fields of input:
%               input.np: number of photon
%               input.dz,input.dr,input.da: the grid size in z, r and alpha
%               dimension
%               input.nz,input.nr,input.na: the number of grids in three
%               dimensions
%               input.nl: number of layers
%               input.n_above,input.n_below: the refractive index of the 
%               two ambients
%               input.layers[m].n: the refractive index of layer m
%               input.layers[m].mua: mua of the layer m
%               input.layers[m].mus: mus of the layer m
%               input.layers[m].g: anitroscopy of layer m
%               input.layers[m].z0 and z1: the beginning and ending z
%               coordinates of layer m
%               input.layers[m].cos_crit0 and cos_crit1: the cosines of the
%               critical angles of the two boundaries of layer m
% Also see: MC

% 
% Author: Qiang Bo
% Last revised: 10/05/2005

function [input] = InitializeInput(inputfilename)

run(inputfilename); % load the input parameters
% data structure for one layer
layer = struct('n',[],... % refractive index
               'mua',[],... % absorption coefficient
               'mus',[],... % scattering coefficient
               'g',[],... % anisotropy
               'z0',[],... % z coordinate of upper boundary
               'z1',[],... % z coordinate of lower boundary
               'cos_crit0',[],... % cosine of critical angle for upper boundary
               'cos_crit1',[]... % cosine of critical angle for lower boundary
               );
% Initialize data structure for multiple layers and set their values
% according to input
layers = repmat(layer,1,nl);
for ilayer = 1:nl
    z1 = sum(d(1:ilayer));
    z0 = z1-d(ilayer);
    n0 = n(ilayer);
    n1 = n(ilayer+1); % This is actually the current layer's n. 
                            % Remember there are (nl+2) elements in n.
    n2 = n(ilayer+2);
    if (n1 > n0)
        cos_crit0 = sqrt(1-(n0^2)/(n1^2));
    else
        cos_crit0 = 0;
    end
    if (n1 > n2)
        cos_crit1 = sqrt(1-(n2^2)/(n1^2));
    else
        cos_crit1 = 0;
    end
    layers(ilayer).n = n1;
    layers(ilayer).mua = mua(ilayer);
    layers(ilayer).mus = mus(ilayer);
    layers(ilayer).g = g(ilayer);
    layers(ilayer).z0 = z0;
    layers(ilayer).z1 = z1;
    layers(ilayer).cos_crit0 = cos_crit0;
    layers(ilayer).cos_crit1 = cos_crit1;
end
% data structure for input
input = struct('np',np,... % number of photons to be simulated
               'dz',dz,... % z grid separation
               'dr',dz,... % r grid separation
               'da',0.5*pi/na,... % alpha grid separation
               'nz',nz,... % array range 1..nz
               'nr',nr,... % array range 1..nr
               'na',na,... % array range 1..na
               'nl',nl,... % number of layers
               'n_above',n(1),... % refractive index for layer above
               'n_below',n(end),... % refractive index for layer below
               'layers',layers... % specifications for layers
                );           