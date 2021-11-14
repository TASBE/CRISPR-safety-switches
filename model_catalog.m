% This file is a list of all models in the order that we'd like to explore and plot them

addpath('generated_models');

% Each cell row is: {name, function}
models = {
    'No Delay'                              @Basic_kill_switch;
    % Singles
    'Activator'                             @Activator_Kill_Switch;
    'Repressor'                             @Repressor_Kill_Switch;
    'Cre-ON'                                @Cre_on_Kill_Switch;
    'Cre-OFF'                               @Cre_off_Kill_Switch;
    % Chains
    'Sequential Activator \rightarrow Activator' @Chain_Activator_Activator_Kill_Switch;
    'Sequential Activator \rightarrow Cre-OFF'   @Chain_Activator_Cre_off_Kill_Switch;
    'Sequential Activator \rightarrow Cre-ON'    @Chain_Activator_Cre_on_Kill_Switch;
    'Sequential Activator \rightarrow Repressor' @Chain_Activator_Repressor_Kill_Switch;
    'Sequential Cre-OFF \rightarrow Activator'   @Chain_Cre_off_Activator_Kill_Switch;
    'Sequential Cre-OFF \rightarrow Cre-OFF'     @Chain_Cre_off_Cre_off_Kill_Switch;
    'Sequential Cre-OFF \rightarrow Cre-ON'      @Chain_Cre_off_Cre_on_Kill_Switch;
    'Sequential Cre-OFF \rightarrow Repressor'   @Chain_Cre_off_Repressor_Kill_Switch;
    'Sequential Cre-ON \rightarrow Activator'    @Chain_Cre_on_Activator_Kill_Switch;
    'Sequential Cre-ON \rightarrow Cre-OFF'      @Chain_Cre_on_Cre_off_Kill_Switch;
    'Sequential Cre-ON \rightarrow Cre-ON'       @Chain_Cre_on_Cre_on_Kill_Switch;
    'Sequential Cre-ON \rightarrow Repressor'    @Chain_Cre_on_Repressor_Kill_Switch;
    'Sequential Repressor \rightarrow Activator' @Chain_Repressor_Activator_Kill_Switch;
    'Sequential Repressor \rightarrow Cre-OFF'   @Chain_Repressor_Cre_off_Kill_Switch;
    'Sequential Repressor \rightarrow Cre-ON'    @Chain_Repressor_Cre_on_Kill_Switch;
    'Sequential Repressor \rightarrow Repressor' @Chain_Repressor_Repressor_Kill_Switch;
    % Joint
    'Parallel Activator/Activator'             @Joint_Activator_Activator_Kill_Switch;
    'Parallel Activator/Cre-OFF'               @Joint_Activator_Cre_off_Kill_Switch;
    'Parallel Activator/Cre-ON'                @Joint_Activator_Cre_on_Kill_Switch;
    'Parallel Activator/Repressor'             @Joint_Activator_Repressor_Kill_Switch;
    'Parallel Cre-OFF/Cre-OFF'                 @Joint_Cre_off_Cre_off_Kill_Switch;
    'Parallel Cre-OFF/Cre-ON'                  @Joint_Cre_off_Cre_on_Kill_Switch;
    'Parallel Cre-OFF/Repressor'               @Joint_Cre_off_Repressor_Kill_Switch;
    'Parallel Cre-ON/Cre-ON'                   @Joint_Cre_on_Cre_on_Kill_Switch;
    'Parallel Cre-ON/Repressor'                @Joint_Cre_on_Repressor_Kill_Switch;
    'Parallel Repressor/Repressor'             @Joint_Repressor_Repressor_Kill_Switch;
    };
n_models = size(models,1);

MODEL_NAME = 1;
MODEL_FUN = 2;