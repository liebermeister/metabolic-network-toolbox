function parameter_balancing_sbtab(model_name, metabolic_data_sbtab, options)

% parameter_balancing_sbtab(model_name, metabolic_data_sbtab, options)
%
% a wrapper for parameter balancing that reads all input data 
% (model, priors, kinetic and metabolic data) from files in SBtab format

% example:
% model_name        = '/home/wolfram/projekte/carbon_fixation/matlab_carbon_fixation/models/model_MOGtranshydrogenase1/MOGtranshydrogenase1';
% file_kinetic_data = '/home/wolfram/projekte/carbon_fixation/matlab_carbon_fixation/models/model_MOGtranshydrogenase1/MOGtranshydrogenase1_Quantity.tsv';

options.kinetics          = 'cs';  % ms, ...
options.parametrisation   = 'equilibrium constant'; 
options.reaction_units    = 'concentration per time'; % 'amount per time'
options.enzymes_explicit  = 1;
options.include_metabolic = 1;

[network, kinetic_data] = sbtab_to_network(model_name,'cs');

quantity_info = data_integration_load_quantity_info; % general quantity table (table of strings)

[model_quantities, basic_quantities, data_quantities] = parameter_balancing_quantities(quantity_info, network, options);

kinetic_data = data_integration_load_kinetic_data(data_quantities, quantity_info, network, file_kinetic_data);

% change some of the prior distributions
quantity_info.PriorMedian{quantity_info.symbol_index.u} = '1';
quantity_info.PriorStd{quantity_info.symbol_index.u}    = '3';
quantity_info.PriorStd{quantity_info.symbol_index.KV}   = '3';


% -----------------------------------------------------------------------
% run simple parameter balancing without metabolic data

task = parameter_balancing_task(network, kinetic_data, quantity_info, model_quantities, basic_quantities);
res  = parameter_balancing(task, quantity_info);

figure(1); clf
subplot(2,1,1); plot(task.q.prior.mean,'c'); hold on; plot(res.q_posterior.mean,'b'); ylabel('Mean'); title('Basic quantities')
set(gca,'XTick',cumsum([0;task.q.numbers(1:end-1)]),'XTickLabel',fieldnames(task.q.indices));
subplot(2,1,2); plot(task.q.prior.std,'c');  hold on; plot(res.q_posterior.std,'b'); ylabel('Std'); 
set(gca,'XTick',cumsum([0;task.q.numbers(1:end-1)]),'XTickLabel',fieldnames(task.q.indices),'YScale','Log');

figure(2); clf
subplot(2,1,1); plot(task.xdata.mean,'r'); hold on; plot(res.xdata_posterior.mean,'b'); ylabel('Mean'); title('Data quantities') 
set(gca,'XTick',cumsum([0;task.xdata.numbers(1:end-1)]),'XTickLabel',fieldnames(task.xdata.indices));
subplot(2,1,2); plot(task.xdata.std,'r');  hold on; plot(res.xdata_posterior.std,'b'); ylabel('Std'); 
set(gca,'XTick',cumsum([0;task.xdata.numbers(1:end-1)]),'XTickLabel',fieldnames(task.xdata.indices),'YScale','Log');

figure(3); clf
subplot(2,1,1); plot(res.xmodel_posterior.mean,'b');ylabel('Mean');  title('Model quantities')
set(gca,'XTick',cumsum([0;task.xmodel.numbers(1:end-1)]),'XTickLabel',fieldnames(task.xmodel.indices));
subplot(2,1,2); plot(res.xmodel_posterior.std,'b'); ylabel('Std'); 
set(gca,'XTick',cumsum([0;task.xmodel.numbers(1:end-1)]),'XTickLabel',fieldnames(task.xmodel.indices),'YScale','Log');


% --------------------------------------------------------------------------------------
% Run single parameter balancing with metabolic data

% insert first steady state and check the resulting flux predictions
% [nr,nm,nx,KM_indices,KA_indices,KI_indices,nKM,nKA,nKI] = network_numbers(network);

metabolic_data.v.mean = data.flux_data.DataMean;
metabolic_data.c.mean = data.metabolite_data.DataMean;
metabolic_data.u.mean = data.protein_data.DataMean;
metabolic_data.v.std  = 0.1 * mean(abs(metabolic_data.v.mean(:))) * ones(size(metabolic_data.v.mean));
metabolic_data.c.std  = 0.1 * mean(abs(metabolic_data.c.mean(:))) * ones(size(metabolic_data.c.mean));
metabolic_data.u.std  = 0.1 * mean(abs(metabolic_data.u.mean(:))) * ones(size(metabolic_data.u.mean));

% insert the first state into the kinetic data

fn = fieldnames(metabolic_data);
for it = 1:length(fn),
  qq = fn{it};
  if isfield(kinetic_data,qq), 
    kinetic_data.(qq).mean = metabolic_data.(qq).mean(:,1);
    kinetic_data.(qq).std  = metabolic_data.(qq).std(:,1);
    if strcmp(kinetic_data.(qq).scaling,'Logarithmic'),
      [kinetic_data.(qq).mean_ln, kinetic_data.(qq).std_ln] = lognormal_normal2log(kinetic_data.(qq).mean,kinetic_data.(qq).std);
    end
  end
end


% insert flux directions for 50 percent of all reactions (abs(v) above median)
my_v    = metabolic_data.v.mean(:,1);
thr     = quantile(abs(my_v),0.5)
epsilon = 0.1;
kinetic_data.A.lower = nan * kinetic_data.A.lower;
kinetic_data.A.upper = nan * kinetic_data.A.upper;
kinetic_data.A.lower(find(my_v> thr)) =  epsilon;
kinetic_data.A.upper(find(my_v<-thr)) = -epsilon;

task = parameter_balancing_task(network, kinetic_data, quantity_info, model_quantities, basic_quantities);
res  = parameter_balancing(task,quantity_info);

v_mean = metabolic_data.v.mean(:,1);
v_std  = metabolic_data.v.std(:,1);
