% Tests that run the examples in the folder mnt/demo

setup_test_path;
examples = {'mnt_example_network', 'mnt_example_fluxes', 'mnt_example_kinetics', ...
            'mnt_example_control_analysis', 'mnt_example_graphics', 'mnt_example_import_export'};
old_setting = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');
restore_setting = onCleanup(@() set(0, 'DefaultFigureVisible', old_setting));

%% All examples run without errors
for it = 1:length(examples)
  evalc(examples{it});
end
close all;
