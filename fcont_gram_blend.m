
% fcont_gram_blend : This routine computes the Fourier continuation and the
% FC coefficients of a sample function.
%
% Inputs: 
%   fx : (real) values of the functions f at nodes x_i
%   d : number of Gram interpolation points 
%   C : number of continuation points
%   A : matrix containing the values of the extended Gram polynomials at 
%   the continuation points
%   Ql : matrix containing the values of the Gram polynomials at the first 
%        d points of the interval. 
%   Qr : matrix containing the values of the Gram polynomials at the first 
%        d points of the interval.
%   BC : 2 by 2 matrix
%       BC(1, :) says if the left and right boundary are Dirichlet 
%       or Neumann : 0 for Dirichlet, 1 for Neumann
%       BC(2, :) gives the boundary condition to set at each Boundary in
%       the case of a Neumann BC
%       h : grid size
%
% Outputs:
%   fx_cont_coeffs : complex vector of size lenfgth(fx) + C containing the
%   FC coefficients of fx
%   fcont : (real) vector of size lenfgth(fx) + C containing the values of
%   the fx and its continuation on the extended grid
%
%   Author : Daniel Leibovici
%   Email : dleibovi@caltech.edu

function [fx_cont_coeffs, fcont, fc_l, fc_r] = fcont_gram_blend(fx, d, C, ... 
    A, Ql, Qr, BC, h)



fourPts = length(fx) + C;
h0 = 1/100;

fr = fx((length(fx) - (d-1)):length(fx));
if (BC(1, 2) == 1)
    fr(end) = BC(2, 2)*h/h0;
end

fl = fx(1:d);
if (BC(1, 1) == 1)
    fl(1) = - BC(2, 1)*h/h0;
end

fc_r = A * (Qr.' *fr);
fc_l = flipud(A * (Ql .' * flipud(fl)));

fc_r = double(fc_r);
fc_l = double(fc_l);

fcont = [fx; fc_l + fc_r];

fx_cont_coeffs = fft(fcont)/fourPts;

return
