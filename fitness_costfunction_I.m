% This function returns the control effort (information) for a given set of
% ResolutionParams.
%
% The form of the reolution parameter vector is n1 m1 n2 m2 ..., where n is
% the number of amplitude levels and m the number of repeated measures. The
% information is then calculated to be
% 
% I = sum_i m_i log2 (n_i)
%
% DH 31.08.2015
% 

function CE_I = fitness_costfunction_I(ResolutionParams)

CE_I = 2 * sum( (1+ResolutionParams(2:2:end)) .* log2(1 + ResolutionParams(1:2:end-1)));

% the factor 2 recognizes, that each resolution pair is used twice in the
% symmetric 