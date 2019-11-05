%
%
%
% Adapted pattern search algorithm to calculate control effort.
%
% Ver 2.0 including visualization
% Ver 2.1 revised structure of the code for more transparency
% Ver 3.0 incuded internal memory to avoid recalculation
% Ver 4.0 error estimation based on stage3 search 06.02.2015
% Ver 4.1 corrected bug in display which caused false "sucess=1" dispalys
% Ver 4.2 revised error estimation 13.03.2015
% Ver 4.3 revised error estimation 31.8.2015 % used for the publication
% currently under review, check for updates later.
% 

function [x_out,f_val_out,funEvals,optim_errors,PS_memory] = PS_CEmin(PS_problem)
%% check consistency of input
if ~length(PS_problem.lb)==PS_problem.nvars; error('size of lower bound vector does not mach the number of variables'); end
if ~length(PS_problem.ub)==PS_problem.nvars; error('size of upper bound vector does not mach the number of variables'); end
% check PS_problem struct
if ~isfield(PS_problem,'memory'); PS_problem.memory = 1; end
% to be continued as required

%% initialize search

% poll initial case
x_poll                     = PS_problem.ub;  % initial starting point at the upper bound
[f_val_poll, nonlcon_poll] = evaluate_x_poll(x_poll,inf,PS_problem,'initialize');  % calculate the value for the fitness function and the constraint
funEvals                   = 1;              % count the numerber of function evaluations
if (nonlcon_poll>0)||isnan(nonlcon_poll); % check for a working upper bound solution.
    error('Nonlinear Constraint of the initial point is violated');
end

% save the sucessful parameters and results
x       = x_poll;
f_val   = f_val_poll;

% display progress and current search parameters
if PS_problem.Display
    disp('Maximum resolution results:')
    display_output(x_poll,f_val_poll,nonlcon_poll,funEvals,[],x,PS_problem,1)
    disp('Stage 1: simulatniously reducing the resolutions')
end


%% stage one:
% searches by systematically reducing the resolutions of all signals at the
% same time. This is done by a bisection algorithm starting in the middle
% between lower and upper bounds.
meshsize = PS_problem.ub - PS_problem.lb;

while max(meshsize)>=1
    % determine next x_poll depending on previous results
    meshsize = meshsize/2;
    if nonlcon_poll<=0  % then the constraint was met and a new better solution was found
        x       = x_poll;            % save the last working resolution parameters x
        f_val   = f_val_poll;
        
        x_poll  = x_poll - meshsize; % go to lower resolutions
    else %unsuccessfull poll
        x_poll  = x_poll + meshsize; % go to higher resolutions
    end
    
    % calculate performance and constraints
    [f_val_poll, nonlcon_poll] = evaluate_x_poll(x_poll,f_val,PS_problem);
    funEvals = funEvals+1;
    
    if PS_problem.Display;
        display_output(x_poll,f_val_poll,nonlcon_poll,funEvals,meshsize,x,PS_problem,(nonlcon_poll<=0 && f_val_poll<f_val));
    end
end

%% stage two
% reduce resolution in each signal individually in a pattern search manner.
% If a limit in one signal is found, also check to increase the reolution
% in this signal again while reducing the resolution in each of the other
% signals (this is tha major improvement compared to classical
% patternsearch as it resembles a linear combination of dimensions)

if PS_problem.Display; disp('Stage 2: Patternsearch'); end

% pattern vectors
pattern_vecs              = -eye(PS_problem.nvars);
meshsize                  = x / 4; %initial mesh size
last_successful_direction = 1; % has to be initialized

% the following loop is run until the meshsize is to small:
while (max(meshsize)>=0.5) && PS_problem.MaxFunEval > funEvals
    poll_successful = false;
    i1 = 1;
    
    % poll for a given mesh size
    while i1<=size(pattern_vecs,1) && PS_problem.MaxFunEval > funEvals
        x_poll          = x + meshsize.*pattern_vecs(i1,:);
        x_poll          = min(x_poll,PS_problem.ub); % respect the upper bounds
        x_poll          = max(x_poll,PS_problem.lb); % respect the lower bounds
        
        [f_val_poll, nonlcon_poll] = evaluate_x_poll(x_poll,f_val,PS_problem);
        funEvals = funEvals+1;
        if (f_val_poll<f_val) && (nonlcon_poll<=0) % check if result is better and the constraints are met
            poll_successful           = true;
            if i1<=PS_problem.nvars % in rare cases, the last successful direction is a linear combination. In this case, the last_successful_direction is not updated
                last_successful_direction = i1; % this direction is saved to generate linear combinations for the poll which allow to go back in a previously sucessful direction.
            end
            x       = x_poll;
            f_val   = f_val_poll;
            i1 = size(pattern_vecs,1) + 1; %exit while loop
        else % constraints are not met
            i1 = i1+1;
        end
        
        if PS_problem.Display;
            display_output(x_poll,f_val_poll,nonlcon_poll,funEvals,meshsize,x,PS_problem,poll_successful)
        end
        
    end
    
    % change mesh size and craete new poll pattern vectors
    if poll_successful  %increase mesh size
        meshsize        = meshsize*PS_problem.MeshExpansion;
        pattern_vecs    = -eye(PS_problem.nvars);
    else                %decrease mesh size
        meshsize        = meshsize*PS_problem.MeshContraction;
        pattern_vecs_1  = -eye(PS_problem.nvars);
        pattern_vecs_1(:,last_successful_direction) = 0.5;
        pattern_vecs    = vertcat(-eye(PS_problem.nvars),pattern_vecs_1);
        clear pattern_vecs_1
    end
end

%% stage 3
% check linear combinations with small variations near the optimum
if PS_problem.Display; disp('Stage 3: Check linear combinations with small variations near the optimum'); end
if PS_problem.Display; disp('Stage 3: meshsize = 1'); end
x = ceil(x); % from this point on forward, only integer manipulations with x

poll_successful = true;
while poll_successful
    if PS_problem.MaxFunEval <= funEvals; break; end % maximum iterations are reached
    meshsize = 1;
    [x, f_val, funEvals, poll_successful, error_x_poll] = stage3(x, f_val, funEvals, meshsize, PS_problem);
end

% if PS_problem.Display; disp('Stage 3: meshsize = 2'); end
% 
% poll_successful = true;
% while poll_successful
%     if PS_problem.MaxFunEval <= funEvals; break; end % maximum iterations are reached
%     meshsize = 2;
%     [x, f_val, funEvals, poll_successful, error_x_poll] = stage3(x, f_val, funEvals, meshsize, PS_problem);
% end

%% estimate error / certainty interval
if PS_problem.Display; disp('Stage 4: Prepare function output'); end

% The f_val_certanty gives the range of tested f_vals which were rejected
% due to the constraint-limit. It can be said that there is no better
% solution in this range for the given poll combinations. The larger the
% number, the better.
optim_errors.f_val_certanty = min(error_x_poll(:,end));

% The x_poll_certanty gives the range for each optim variable which was
% tested in stage 3
optim_errors.x_poll_certanty = x - min(error_x_poll(:,1:end-1),[],1);

% calculate f_val_error:
% the maximum difference found in f_val for the smallest vatiations in x
% The error, thus, specifies the interval (of f_val) in which a better
% solution may be found if the resolution of the poll search (meshsize) is
% refined. The smaller the number the better.
diff2_fval = nan(PS_problem.nvars,1); % preallocate
for i1 = 1:PS_problem.nvars % go through each optimization variable
    x_poll     = x;
    x_poll(i1) = x_poll(i1)-1; % change each channel by -1. This corresponds to the smallest difference tested in each channel
    x_poll     = min(x_poll,PS_problem.ub); % respect the upper bounds
    x_poll     = max(x_poll,PS_problem.lb); % respect the lower bounds
    [f_val_poll, ~] = evaluate_x_poll(x_poll,f_val,PS_problem); % retrieve f_val (from memory)
    diff2_fval(i1) = f_val-f_val_poll; % all > 0
end
optim_errors.f_val_error  = max(diff2_fval);

%% function output
x_out     = x;
f_val_out = f_val;
if PS_problem.memory % if each calculation was saved
    [~, ~, PS_memory] = evaluate_x_poll(x_poll,inf,PS_problem,'final');
    if PS_problem.memory_save % resluts are saved to disk at a certain interval
        save(PS_problem.memory_save_name,'PS_memory')
    end
else
    PS_memory = [];
end

%% display sucess
disp(['*** optimization completed after ' num2str(funEvals) 'iterations  ***'])
if funEvals==PS_problem.MaxFunEval
    disp('*** Maximum number of iterations reached ***')
    disp('*** Error estimation may not be adequate if optimization is stoped due to maximum number of iterations ***')
end
disp(['minimum f_val=' num2str(f_val_out)])
disp(x_out)
disp('***')

end

function [x, f_val, funEvals, poll_successful, error_x_poll] = stage3(x, f_val, funEvals, meshsize, PS_problem)

error_x_poll = [x f_val]; %this array will contain all tested solutions to estimate the error

poll_successful = false; %until a better solution is found, this holds

for i1 = 1:PS_problem.nvars
    for i2 = 1:PS_problem.nvars
        for i3 = 1:PS_problem.nvars
            x_poll = x;
            x_poll(i1) = x_poll(i1) + meshsize;
            x_poll(i2) = x_poll(i2) - meshsize;
            x_poll(i3) = x_poll(i3) - meshsize;
            x_poll     = min(x_poll,PS_problem.ub); % respect the upper bounds
            x_poll     = max(x_poll,PS_problem.lb); % respect the lower bounds
            
            % calculate f_val for x_poll
            [f_val_poll, nonlcon_poll] = evaluate_x_poll(x_poll,f_val,PS_problem);
            funEvals = funEvals+1; % increase the counter of function evaluations
            
            %save for error calculation
            error_x_poll = unique(vertcat(error_x_poll,[x_poll f_val_poll]),'rows');
            
            if (f_val_poll<f_val) && (nonlcon_poll<=0) % check if result is better and the constraints are met
                poll_successful           = true;
                x                         = x_poll;
                f_val                     = f_val_poll;
            end
            
            if PS_problem.Display;
                display_output(x_poll,f_val_poll,nonlcon_poll,funEvals,meshsize,x,PS_problem,poll_successful)
            end
            
            if poll_successful; return; end % start over with stage 3 if a better solution was found
            if PS_problem.MaxFunEval <= funEvals; return; end % maximum iterations are reached
        end
    end
end


end

%%
function [f_val_poll, nonlcon_poll, varargout] = evaluate_x_poll(x_poll,f_val,PS_problem,varargin)

x_poll = ceil(x_poll);

persistent PS_memory
if PS_problem.memory
    
    if ~isempty(varargin) % initialize at first run or return memory at final call
        if strcmp(varargin{1}, 'initialize')
            if isfield(PS_problem,'memory_initial');
                PS_memory = PS_problem.memory_initial;
            else %calculate first
                f_val_poll   = PS_problem.fitnessfcn(x_poll);
                nonlcon_poll = PS_problem.nonlcon(x_poll);
                PS_memory    = [x_poll f_val_poll nonlcon_poll];
                return
            end
        elseif strcmp(varargin{1}, 'final')
            f_val_poll   = [];
            nonlcon_poll = [];
            varargout = {PS_memory};
            return
        else
            error ('optional argument to function evaluate_x_poll is not valid')
        end
    end
    
    % check, if the current x_poll has been calculated before
    previously_calculated = PS_memory(:, 1) == x_poll(1);
    for i1 = 2:(size(PS_memory, 2)-2);
        previously_calculated(previously_calculated) = PS_memory(previously_calculated, i1) == x_poll(i1);
    end
    previously_calculated = find(previously_calculated);
    if isempty(previously_calculated);
        previously_calculated = 0;
    end
    
    if previously_calculated >= 1 % then retrieve previous results and return
        f_val_poll   = PS_memory(previously_calculated, end-1);
        nonlcon_poll = PS_memory(previously_calculated, end);
        return
        % else continue and calculate
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate fitness and constraint %
f_val_poll       = PS_problem.fitnessfcn(x_poll);

if f_val_poll < f_val %only check for constraints if a better f_val was found
    nonlcon_poll = PS_problem.nonlcon(x_poll);
else % if the solution is not better, don't calculate the constraint
    nonlcon_poll = inf; % meaning, that it has not been caclulated
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if PS_problem.memory % the results are stored in memory to avoid recalcultation
    PS_memory = [PS_memory ; [x_poll f_val_poll nonlcon_poll]];
    if PS_problem.memory_save && (mod(size(PS_memory,1),PS_problem.memory_save_interv)==0) % resluts are saved to disk at a certain interval
        save(PS_problem.memory_save_name,'PS_memory')
    end
end

end

%%
function display_output(x_poll,f_val_poll,nonlcon_poll,funEvals,meshsize,x,PS_problem,poll_successful)
x      = ceil(x);
x_poll = ceil(x_poll);
% Output to commandline
fprintf('Iter=%5d Success=%1d F_val=%9.2e Constr=%9.2e MaxMeshsize=%.2e x='...
    ,funEvals,poll_successful,f_val_poll,nonlcon_poll,max(meshsize));fprintf('%5d ',x_poll);fprintf('\n');

persistent fh_display
if PS_problem.Display == 2
    
    % Output to figure
    % Step 1: visualize grid
    M_grid = ones(max(PS_problem.ub)-min(PS_problem.lb)+1,size(PS_problem.lb,2));
    size_M_grid = size(M_grid);

    if isempty(fh_display)
        fh_display = figure; % initialize a figure for the plot
    end
    figure(fh_display)
    clf
    hold on
    spy(M_grid,'bo',10)
    set(gca,'YDir','normal');
    
    % Step 2: bounds
    scatter(1:size_M_grid(2),PS_problem.lb,200,[0.5 0.5 0.5],'filled')
    scatter(1:size_M_grid(2),PS_problem.ub,200,[0.5 0.5 0.5],'filled')
    
    % Step 3: tested cases
    % spy(M_all,'b',30)
    
    % Step 4: last working
    scatter(1:size_M_grid(2),x,130,[0 .8 0],'filled')
    
    % Step 5: current tested
    scatter(1:size_M_grid(2),x_poll,50,[.9 0 0],'filled')
    
    drawnow
end
end