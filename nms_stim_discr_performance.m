%
% Function to calculate the performance in the walking model for the STIM
% scenarion (see paper for more details)
% The corresponding simulink model is required
% 
function Perf = nms_stim_discr_performance(ResolutionParams, save_pathname)

% load model parameters
nms_model_MechInit
nms_model_ControlInit

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vary parameters for optimization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Signal resolution parameters for ControlEffort Minimization
MP.contr.HFL_du = 1/ResolutionParams(1);
MP.contr.HFL_dt = 10/ResolutionParams(2);%MP.T_end/ResolutionParams(2);
MP.contr.GLU_du = 1/ResolutionParams(3);
MP.contr.GLU_dt = 10/ResolutionParams(4);
MP.contr.HAM_du = 1/ResolutionParams(5);
MP.contr.HAM_dt = 10/ResolutionParams(6);
MP.contr.VAS_du = 1/ResolutionParams(7);
MP.contr.VAS_dt = 10/ResolutionParams(8);
MP.contr.GAS_du = 1/ResolutionParams(9);
MP.contr.GAS_dt = 10/ResolutionParams(10);
MP.contr.TA_du  = 1/ResolutionParams(11);
MP.contr.TA_dt  = 10/ResolutionParams(12);
MP.contr.SOL_du = 1/ResolutionParams(13);
MP.contr.SOL_dt = 10/ResolutionParams(14);

%% %%%%%%%%%%%
% Simulation %
%%%%%%%%%%%%%%
try % catch simulation errors in order to avoid optimization stop
    % simulation with optimization parameters
    SimData=sim('nms_stim_discr_model','SrcWorkspace','current','SimulationMode', 'normal');
catch SimError
    Perf = nan;
    display(SimError.message)
    return
end

%% %%%%%%%%%%%%%%%%%%%%
% calculating fitness %
%%%%%%%%%%%%%%%%%%%%%%%

SimData_Stop = SimData.get('SimData_Stop');
Perf         = MP.T_end-SimData_Stop.Time(end)-MP.ts_max;

%% %%%%%%%%%%%%%%%%%%%%%%
% save data if required %
%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('save_pathname','var')
    SimData_Kinematics = SimData.get('SimData_Kinematics');
    SimData_Torques    = SimData.get('SimData_Torques');
    SimData_Animation  = SimData.get('SimData_Animation');
    SimData_Stop       = SimData.get('SimData_Stop');
    save(save_pathname)
end
%% output for debugging
%fitness_plotscript
%fprintf('Perf= %9.2e ', Perf);
%fprintf('ResolutionParams= '); fprintf(' %5d ', ResolutionParams);
%fprintf('ControlParams= ')   ; fprintf(' %9.2e ', ControlParams);
%fprintf('\n');