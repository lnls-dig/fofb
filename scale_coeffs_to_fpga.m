function [sos,g] = scale_coeffs_to_fpga(sos,g,max_fpga_coeff)
% BETTER_USE_FPGA_COEFFS  Scale coefficients to better use FPGA resolution.
%
% Scale the coefficients of a given SOS system to better use FPGA resolution.
% Adjust the overall system gain accordingly.
%
% INPUTS:
%   sos:            Second-order section representation (L-by-6 matrix, where L
%                   is the number of systems).
%   g:              Overall system gain.
%   max_fpga_coeff: The maximum coefficient that can be represented in FPGA.
%
% OUTPUTS:
%   sos:            Second-order section representation (L-by-6 matrix, where L
%                   is the number of systems) after scaling.
%   g:              Overall system gain after scaling.

% NOTE: The denominator can't be changed because gateware assumes a0 = 1.

factors = max_fpga_coeff./max(abs(sos(:,1:3))');
sos(:,1:3) = diag(factors)*sos(:,1:3);
g = g/prod(factors);

end
