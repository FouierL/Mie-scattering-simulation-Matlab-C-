% MieExample.m

clc,clear
fprintf('=============\n')

lambda = 0.643; % um, wavelength of light in vacuo
dia = 0.2; % um, diameter of sphere
density = 180; % #/um^3, density of spheres
npar = 1.59; % particle refractive index (polystyrene)
nmed = 1.33; % medium refractive index (water)

m = npar/nmed + 0i; 
x = pi*dia/(lambda/nmed);

u = Mie(m, x)'; % <----- Matlzer's subroutine

qsca = u(5);
g    = u(8);

A       = pi*dia^2/4; % um^2
sigma_s = qsca*A; % um^2
mus     = sigma_s*density*1e4; % cm^-1
musp    = mus*(1-g); 

fprintf('--1.-- choice:\n')
fprintf(sprintf('lambda  \t= %0.3f um\n', lambda))
fprintf(sprintf('dia     \t= %0.3f um\n', dia))
fprintf(sprintf('density \t= %0.3f um\n', density))
fprintf('--2.-- result:\n')
fprintf(sprintf('real(m) \t= %0.3f\n', u(1)))
fprintf(sprintf('imag(m) \t= %0.3e\n', u(2)))
fprintf(sprintf('x       \t= %0.3e\n', u(3)))
fprintf(sprintf('qext    \t= %0.3e\n', u(4)))
fprintf(sprintf('qsca    \t= %0.3e\n', u(5)))
fprintf(sprintf('qabs    \t= %0.3e\n', u(6)))
fprintf(sprintf('qb      \t= %0.3e\n', u(7)))
fprintf(sprintf('asy       \t= %0.4f\n', u(8)))
fprintf(sprintf('qratio  \t= %0.3e\n', u(9)))
fprintf('--3.-- optical properties:\n')
fprintf(sprintf('mus  \t= %0.3f cm^-1\n', mus))
fprintf(sprintf('g       \t= %0.4f\n', g))
fprintf(sprintf('musp  \t= %0.3f cm^-1\n', musp))