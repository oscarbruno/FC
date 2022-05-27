
% fc_der - Routine to compute the derivative of a function using Fourier
% continuation.
%
%
%
% Inputs:
%   fx : (real) values of the functions f at nodes x_i
%   der_coeffs : complex vector 2*1i*pi/prd * k
%   filter : real vector containing the spectral filtering coefficients
%   d : number of Gram interpolation points
%   C : number of continuation points
%   AQ_r: matrix product (A * Q) used to obtain the continuation from the
%   last d nodes
%   FAQ_lF: matrix product (A * Q) used to obtain the continuation from the
%   first d nodes
%   BC: 2x2 matrix specifying the type of boundary condition at each
%   interval end (Dirichlet or Neumann), and the associated boundary value
%   for a Neumann boundary condition
%   h : grid size
%
% Output:
%   fx_der : Real vector of the same size as fx, containing the values of
%   the spatial derivative of f at nodes x_i computed with the FC procedure
%
%
% Author:
%
%   Daniel Leibovici
%   email: dleibovi@caltech.edu
%

function [fx_der, fc_coeffs] = fc_der(fx, der_coeffs, filter, d, C, ... 
    AQ_r, FAQ_lF, BC, h)

fc_coeffs = fcont_gram_blend(fx, d, C, AQ_r, FAQ_lF, BC, h);
n = length(fx);
fourPts = n + C;
fx_der = real(ifft(fc_coeffs .* filter .* der_coeffs))*fourPts;
fx_der = fx_der(1 : n);

return
