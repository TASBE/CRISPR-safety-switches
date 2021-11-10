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
    'Chain Cre-ON \rightarrow Activator'    @Chain_Cre_TF_Activator_Cre_on_Kill_Switch;
    'Chain Cre-OFF \rightarrow Activator'   @Chain_Cre_TF_Activator_Cre_off_Kill_Switch;
    'Chain Cre-ON \rightarrow Repressor'    @Chain_Cre_TF_Repressor_Cre_on_Kill_Switch;
    'Chain Cre-OFF \rightarrow Repressor'   @Chain_Cre_TF_Repressor_Cre_off_Kill_Switch;
    'Chain Activator \rightarrow Cre-ON'    @Chain_TF_Cre_Activator_Cre_on_Kill_Switch;
    'Chain Activator \rightarrow Cre-OFF'   @Chain_TF_Cre_Activator_Cre_off_Kill_Switch;
    'Chain Repressor \rightarrow Cre-ON'    @Chain_TF_Cre_Repressor_Cre_on_Kill_Switch;
    'Chain Repressor \rightarrow Cre-OFF'   @Chain_TF_Cre_Repressor_Cre_off_Kill_Switch;
    % Joint
    'Joint Activator/Cre-ON'                @Joint_Activator_Cre_on_Kill_Switch;
    'Joint Activator/Cre-OFF'               @Joint_Activator_Cre_off_Kill_Switch;
    'Joint Repressor/Cre-ON'                @Joint_Repressor_Cre_on_Kill_Switch;
    'Joint Repressor/Cre-OFF'               @Joint_Repressor_Cre_off_Kill_Switch;
    };
n_models = size(models,1);

MODEL_NAME = 1;
MODEL_FUN = 2;