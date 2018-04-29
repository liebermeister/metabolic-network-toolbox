function [network, r, r_orig, kinetic_data, r_samples, parameter_prior, r_mean, r_std, r_geom_mean, r_geom_std] = parameter_balancing_sbtab(model_file, data_file, pb_options)

%  [network, r, r_orig, kinetic_data, r_samples, parameter_prior, r_mean, r_std] = parameter_balancing_sbtab(model_file, data_file, pb_options)
%
% Wrapper function for parameter balancing 
% With no options set, the function uses directly the information from input files; it is best suited for comparisons with other PB tools
%
%  o reads input data (model, priors, kinetic and metabolic data) from SBtab files
%  o runs parameter balancing (without and with concentration data)
% 
% Function arguments:
%  o model_file:                   SBtab or SBML model filename; files with ".xml" extension are assumed to be SBML, otherwise SBtab
%  o data_file:            SBtab data file (kinetic and other constants)
%  o pb_options.parameter_prior_file: File with parameter priors (optional)
% 
% pb_options: struct with options; for default values, see 'parameter_balancing_default_options'
%
% This function calls the functions 'parameter_balancing_task', 'parameter_balancing', and 'parameter_balancing_output'

pb_options = join_struct(parameter_balancing_default_options, pb_options);

% ----------------------------------------------------------
% load model
  
if strcmp(model_file(end-3:end), '.xml'),
  network = network_sbml_import(model_file);
else,
  network = sbtab_to_network(model_file, struct('kinetic_law','cs'));
end

if isfield(network,'metabolite_is_an_enzyme'),
  if sum(network.metabolite_is_an_enzyme),
    warning('At least one of the model species appears to be an enzyme. Parameter balancing currently cannot handle this case. Please remove all enzyme species from your model and run parameter balancing again.')
  end
end

% ----------------------------------------------------------
% load parameter_prior and define relevant quantities

parameter_prior = parameter_balancing_prior([],pb_options.parameter_prior_file,1); 
parameter_prior = pb_parameter_prior_adjust(parameter_prior, pb_options); 

[model_quantities, basic_quantities, data_quantities, pseudo_quantities] = parameter_balancing_quantities(parameter_prior, network, pb_options);

% ----------------------------------------------------------
% load kinetic data

if length(data_file), display(sprintf('  Using data file %s', data_file)); end  

kinetic_data      = data_integration_load_kinetic_data(data_quantities, parameter_prior, network, data_file, 1, 0);
kinetic_data_orig = kinetic_data;
kinetic_data      = pb_kinetic_data_adjust(kinetic_data, parameter_prior, network, pb_options);

% parameter_balancing_kinetic_data_show(kinetic_data);

% -----------------------------------------------------------
% run simple parameter balancing without metabolic data

task = parameter_balancing_task(network, kinetic_data, parameter_prior, model_quantities, basic_quantities, pseudo_quantities);

pm(task.q.prior.mean,task.q.names)
pm(task.q.prior.std.^2,task.q.names)
pm(task.xdata.mean, task.xdata.names)
pm(task.xdata.std.^2, task.xdata.names)

res                         = parameter_balancing_calculation(task, parameter_prior,pb_options);
[r,r_mean,r_std,r_geom_mean,r_geom_std,r_orig,r_samples] = parameter_balancing_output(res,kinetic_data_orig,pb_options);
network.kinetics            = r;


% ===========================================================
% Test graphics

% figure(1); clf
% subplot(2,1,1); plot(task.q.prior.mean,'c'); hold on; plot(res.q_posterior.mean,'b'); ylabel('Mean'); title('Basic quantities')
% set(gca,'XTick',cumsum([0;task.q.numbers(1:end-1)]),'XTickLabel',fieldnames(task.q.indices));
% subplot(2,1,2); plot(task.q.prior.std,'c');  hold on; plot(res.q_posterior.std,'b'); ylabel('Std'); 
% set(gca,'XTick',cumsum([0;task.q.numbers(1:end-1)]),'XTickLabel',fieldnames(task.q.indices),'YScale','Log');
% 
% figure(2); clf
% subplot(2,1,1); plot(task.xdata.mean,'r'); hold on; plot(res.xdata_posterior.mean,'b'); ylabel('Mean'); title('Data quantities') 
% set(gca,'XTick',cumsum([0;task.xdata.numbers(1:end-1)]),'XTickLabel',fieldnames(task.xdata.indices));
% subplot(2,1,2); plot(task.xdata.std,'r');  hold on; plot(res.xdata_posterior.std,'b'); ylabel('Std'); 
% set(gca,'XTick',cumsum([0;task.xdata.numbers(1:end-1)]),'XTickLabel',fieldnames(task.xdata.indices),'YScale','Log');
% 
% figure(3); clf
% subplot(2,1,1); plot(res.xmodel_posterior.mean,'b');ylabel('Mean');  title('Model quantities')
% set(gca,'XTick',cumsum([0;task.xmodel.numbers(1:end-1)]),'XTickLabel',fieldnames(task.xmodel.indices));
% subplot(2,1,2); plot(res.xmodel_posterior.std,'b'); ylabel('Std'); 
% set(gca,'XTick',cumsum([0;task.xmodel.numbers(1:end-1)]),'XTickLabel',fieldnames(task.xmodel.indices),'YScale','Log');


% % --------------------------------------------------------------------------------------
% % Run single parameter balancing with metabolic data
% 
% % insert first steady state and check the resulting flux predictions
% % [nr,nm,nx,KM_indices,KA_indices,KI_indices,nKM,nKA,nKI] = network_numbers(network);
% 
% metabolic_data.v.mean = data.flux_data.DataMean;
% metabolic_data.c.mean = data.metabolite_data.DataMean;
% metabolic_data.u.mean = data.protein_data.DataMean;
% metabolic_data.v.std  = 0.1 * mean(abs(metabolic_data.v.mean(:))) * ones(size(metabolic_data.v.mean));
% metabolic_data.c.std  = 0.1 * mean(abs(metabolic_data.c.mean(:))) * ones(size(metabolic_data.c.mean));
% metabolic_data.u.std  = 0.1 * mean(abs(metabolic_data.u.mean(:))) * ones(size(metabolic_data.u.mean));
% 
% % insert the first state into the kinetic data
% 
% fn = fieldnames(metabolic_data);
% for it = 1:length(fn),
%   qq = fn{it};
%   if isfield(kinetic_data,qq), 
%     kinetic_data.(qq).mean = metabolic_data.(qq).mean(:,1);
%     kinetic_data.(qq).std  = metabolic_data.(qq).std(:,1);
%     if strcmp(kinetic_data.(qq).scaling,'Logarithmic'),
%       [kinetic_data.(qq).mean_ln, kinetic_data.(qq).std_ln] = lognormal_normal2log(kinetic_data.(qq).mean,kinetic_data.(qq).std);
%     end
%   end
% end
% 
% 
% % insert flux directions for 50 percent of all reactions (abs(v) above median)
% my_v    = metabolic_data.v.mean(:,1);
% thr     = quantile(abs(my_v),0.5)
% epsilon = 0.1;
% kinetic_data.A.lower = nan * kinetic_data.A.lower;
% kinetic_data.A.upper = nan * kinetic_data.A.upper;
% kinetic_data.A.lower(find(my_v> thr)) =  epsilon;
% kinetic_data.A.upper(find(my_v<-thr)) = -epsilon;
% 
% task = parameter_balancing_task(network, kinetic_data, parameter_prior, model_quantities, basic_quantities);
% res  = parameter_balancing(task,parameter_prior);
% 
% v_mean = metabolic_data.v.mean(:,1);
% v_std  = metabolic_data.v.std(:,1);
