% This function calculates the optimal performance for a given set of
% ResolutionParams. It optimizes the control parameters to find the minimum
% value for Perf.
%
% DH 15.08.2014 30.8.2014 30.09.2014 4.11.2014

function Perf = fitness_nonlcon(ResolutionParams,performance_fcn, PS_problem_name)

% to save time:
% check if this case has been calculated before
OPFilename = fullfile( [PS_problem_name '_results_rp'] , ['rp_' regexprep(num2str(ResolutionParams),' ','_') '.mat']);
try
    if exist(OPFilename,'file')
        load(OPFilename) % this file contails 'ResolutionParams', 'ControlParams', 'CE_I', 'Perf'
        return
    end
catch errormessage
    warning('could not open file with previous resutls')
    disp(errormessage.message)
end

Perf = performance_fcn(ResolutionParams);

%% save the result to avoid double calculation
% Additionally calculate the information to save all important results
CE_I = fitness_costfunction_I(ResolutionParams);

try
    save(OPFilename, 'ResolutionParams', 'CE_I', 'Perf')
catch errormessage
    warning('could not save resutls for later reuse')
    disp(errormessage.message)
end