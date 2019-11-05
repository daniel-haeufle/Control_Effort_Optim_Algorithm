% This is an example file for how to setup the optimization problem to
% determine control effort. For more details, please see our paper:
% currently under review, check for updates later.
%

clear

%% initialize problem
PS_problem.name       = 'nms_stim_discr';            % a name for the problem
PS_problem.nvars      = 14;                          % number of variables
PS_problem.perf_fcn   = @nms_stim_discr_performance; % performance function
PS_problem.nonlcon    = @(ResolutionParams)fitness_nonlcon(ResolutionParams, PS_problem.perf_fcn, PS_problem.name);            % constraint function
PS_problem.fitnessfcn = @fitness_costfunction_I;     % fitness function of the control effort
PS_problem.lb         = 1 * ones([1, PS_problem.nvars]);           % lower bound
PS_problem.ub         = 2^15 * ones([1, PS_problem.nvars]);        % upper bound, this is an initial guess and the starting vector of the optimization

PS_problem.MeshExpansion      = 2;                   % expansion after a sucessfull poll
PS_problem.MeshContraction    = 0.5;                 % contraction rate if no better solution was polled

PS_problem.MaxFunEval         = inf;  % Maximum number of function evaluations
PS_problem.Display            = 1;                   % 0: only display final result; 1: iterative display output to command line; 2: additionally display figure
PS_problem.memory             = 1;                   % 1: use internal memory to avoid recalcuations
PS_problem.memory_save        = 1;                   % 1: save pattern search memory to disk
PS_problem.memory_save_interv = 10;                  % interval at wich to save to disk (positive integer)
PS_problem.memory_save_name   = [PS_problem.name 'PS_memory_backup'];  % file name

%% turn off warning
warning('off','Simulink:blocks:TDelayBufferTooSmall')

%% %%%%%%%%%%%%%
% optimization %
%%%%%%%%%%%%%%%%
diary([PS_problem.name '_diary.txt'])
addpath('../')
mkdir([PS_problem.name '_results_rp'])

tic;
% start thie optimization:
[ResolutionParams, ControlEffort, PS_funEvals, PS_errors, PS_memory] = PS_CEmin(PS_problem);
optim_duration_in_seconds = toc;
save([PS_problem.name '_results_PS_optim'])

PS_problem.perf_fcn(ResolutionParams, [PS_problem.name '_results_PS_optim_detail']);
PS_problem.perf_fcn(ResolutionParams, [PS_problem.name '_results_PS_optim_detail_50s']);

diary off