function [metabolites,actions] = network_adjust_names_for_sbml_export(metabolites,actions)

% [metabolites,actions] = network_adjust_names_for_sbml_export(metabolites,actions)

eval(default('actions','[]'));

for it=1:length(metabolites),
  metabolites{it} = strrep(metabolites{it},'''','_');
  metabolites{it} = strrep(metabolites{it},'-','_');
  metabolites{it} = strrep(metabolites{it},':','_');
  metabolites{it} = strrep(metabolites{it},',','_');
  metabolites{it} = strrep(metabolites{it},'.','_');
  metabolites{it} = strrep(metabolites{it},' ','_');
  metabolites{it} = strrep(metabolites{it},'+','plus');
  metabolites{it} = strrep(metabolites{it},'(','_');
  metabolites{it} = strrep(metabolites{it},')','_');
  metabolites{it} = strrep(metabolites{it},'[','_');
  metabolites{it} = strrep(metabolites{it},']','_');
  metabolites{it} = strrep(metabolites{it},'<','_');
  metabolites{it} = strrep(metabolites{it},'/','_');  
  metabolites{it} = strrep(metabolites{it},'[','_');  
  metabolites{it} = strrep(metabolites{it},']','_');  
  
  if length(strfind('0123456789',metabolites{it}(1)))
    metabolites{it} = ['_' metabolites{it}];
  end

end

for it=1:length(actions),
  actions{it} = strrep(actions{it},':','_');
  actions{it} = strrep(actions{it},' ','_');
  actions{it} = strrep(actions{it},'=','-');
  actions{it} = strrep(actions{it},'<->','_rev_');
  actions{it} = strrep(actions{it},'->','_irrev_');
  actions{it} = strrep(actions{it},'-','_');
  actions{it} = strrep(actions{it},',','_');
  actions{it} = strrep(actions{it},'.','_');
  actions{it} = strrep(actions{it},'(','_');
  actions{it} = strrep(actions{it},')','_');
  actions{it} = strrep(actions{it},'.','_');
  actions{it} = strrep(actions{it},'+','plus');
  actions{it} = strrep(actions{it},'[','_');
  actions{it} = strrep(actions{it},']','_');

  if length(strfind('0123456789',actions{it}(1)))
    actions{it} = ['_' actions{it}];
  end

end
