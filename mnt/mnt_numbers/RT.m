function result = RT(temperature)

% result = RT(temperature)

eval(default('temperature','300'));  % temperature in Kelvin

constants = ncsp_constants;

result = temperature * constants.R;