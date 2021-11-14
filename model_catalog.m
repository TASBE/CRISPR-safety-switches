% This file is a list of all models in the order that we'd like to explore and plot them

addpath('generated_models');

% Each cell row is: {name, function}
models = {
    'No Delay'                              @Basic_kill_switch                              []              []
    % Singles
    'Activator'                             @Activator_Kill_Switch                          'alpha_p_TF'    []
    'Repressor'                             @Repressor_Kill_Switch                          'alpha_p_TF'    []
    'Cre-ON'                                @Cre_on_Kill_Switch                             'alpha_p_Cre'   []
    'Cre-OFF'                               @Cre_off_Kill_Switch                            'alpha_p_Cre'   []
    % Chains
    'Sequential Activator \rightarrow Activator' @Chain_Activator_Activator_Kill_Switch     'alpha_p_TF'    'alpha_p_TF2'
    'Sequential Activator \rightarrow Cre-OFF'   @Chain_Activator_Cre_off_Kill_Switch       'alpha_p_TF'    'alpha_p_Cre'
    'Sequential Activator \rightarrow Cre-ON'    @Chain_Activator_Cre_on_Kill_Switch        'alpha_p_TF'    'alpha_p_Cre'
    'Sequential Activator \rightarrow Repressor' @Chain_Activator_Repressor_Kill_Switch     'alpha_p_TF'    'alpha_p_TF2'
    'Sequential Cre-OFF \rightarrow Activator'   @Chain_Cre_off_Activator_Kill_Switch       'alpha_p_TF'    'alpha_p_Cre'
    'Sequential Cre-OFF \rightarrow Cre-OFF'     @Chain_Cre_off_Cre_off_Kill_Switch         'alpha_p_Cre'   'alpha_p_CreH'
    'Sequential Cre-OFF \rightarrow Cre-ON'      @Chain_Cre_off_Cre_on_Kill_Switch          'alpha_p_Cre'   'alpha_p_CreH'
    'Sequential Cre-OFF \rightarrow Repressor'   @Chain_Cre_off_Repressor_Kill_Switch       'alpha_p_TF'    'alpha_p_Cre'
    'Sequential Cre-ON \rightarrow Activator'    @Chain_Cre_on_Activator_Kill_Switch        'alpha_p_TF'    'alpha_p_Cre'
    'Sequential Cre-ON \rightarrow Cre-OFF'      @Chain_Cre_on_Cre_off_Kill_Switch          'alpha_p_Cre'   'alpha_p_CreH'
    'Sequential Cre-ON \rightarrow Cre-ON'       @Chain_Cre_on_Cre_on_Kill_Switch           'alpha_p_Cre'   'alpha_p_CreH'
    'Sequential Cre-ON \rightarrow Repressor'    @Chain_Cre_on_Repressor_Kill_Switch        'alpha_p_TF'    'alpha_p_Cre'
    'Sequential Repressor \rightarrow Activator' @Chain_Repressor_Activator_Kill_Switch     'alpha_p_TF'    'alpha_p_TF2'
    'Sequential Repressor \rightarrow Cre-OFF'   @Chain_Repressor_Cre_off_Kill_Switch       'alpha_p_TF'    'alpha_p_Cre'
    'Sequential Repressor \rightarrow Cre-ON'    @Chain_Repressor_Cre_on_Kill_Switch        'alpha_p_TF'    'alpha_p_Cre'
    'Sequential Repressor \rightarrow Repressor' @Chain_Repressor_Repressor_Kill_Switch     'alpha_p_TF'    'alpha_p_TF2'
    % Joint
    'Parallel Activator/Activator'             @Joint_Activator_Activator_Kill_Switch       'alpha_p_TF'    'alpha_p_TF2'
    'Parallel Activator/Cre-OFF'               @Joint_Activator_Cre_off_Kill_Switch         'alpha_p_TF'    'alpha_p_Cre'
    'Parallel Activator/Cre-ON'                @Joint_Activator_Cre_on_Kill_Switch          'alpha_p_TF'    'alpha_p_Cre'
    'Parallel Activator/Repressor'             @Joint_Activator_Repressor_Kill_Switch       'alpha_p_TF'    'alpha_p_TF2'
    'Parallel Cre-OFF/Cre-OFF'                 @Joint_Cre_off_Cre_off_Kill_Switch           'alpha_p_Cre'   'alpha_p_CreH'
    'Parallel Cre-OFF/Cre-ON'                  @Joint_Cre_off_Cre_on_Kill_Switch            'alpha_p_Cre'   'alpha_p_CreH'
    'Parallel Cre-OFF/Repressor'               @Joint_Cre_off_Repressor_Kill_Switch         'alpha_p_TF'    'alpha_p_Cre'
    'Parallel Cre-ON/Cre-ON'                   @Joint_Cre_on_Cre_on_Kill_Switch             'alpha_p_Cre'   'alpha_p_CreH'
    'Parallel Cre-ON/Repressor'                @Joint_Cre_on_Repressor_Kill_Switch          'alpha_p_TF'    'alpha_p_Cre'
    'Parallel Repressor/Repressor'             @Joint_Repressor_Repressor_Kill_Switch       'alpha_p_TF'    'alpha_p_TF2'
    };
n_models = size(models,1);

MODEL_NAME = 1;
MODEL_FUN = 2;
MODEL_PARAM1 = 3;
MODEL_PARAM2 = 4;
