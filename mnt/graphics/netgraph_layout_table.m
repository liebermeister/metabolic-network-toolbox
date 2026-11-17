function layout_table = netgraph_layout_table(networkG)

% NETGRAPH_LAYOUT_TABLE SBtab table with the layout of a network for graphics
%
% layout_table = netgraph_layout_table(networkG)
%
% SBtab table (TableID 'NetworkLayout', TableType 'Position') with the layout of networkG (network for
% graphics, e.g. from netgraph_make_layout_network), to be stored in an SBtab model file together with
% the full model (see mnt_save_model_and_state). Columns: !Element (metabolite or reaction ID, adjusted
% to SBML syntax as in network_to_sbtab), !Position:x, !Position:y, !IsFixed, !IsInvisible.
% Only the elements of networkG are listed; metabolites and reactions of the full model that are not
% listed (e.g. cofactors) are not shown in network graphics.
%
% The TableID is 'NetworkLayout' (not 'Position'), so that sbtab_to_network does not apply the layout
% to the full network. The network for graphics is rebuilt by netgraph_layout_from_table.
%
% (Until 2026-11-17 the IMF function imf_layout_table.)
%
% See also: netgraph_layout_from_table, netgraph_print_positions

% element ids adjusted as in network_to_sbtab (option adjust_metabolite_ids, default), so that they
% match the ids in the SBtab model tables
[metabolites, actions] = network_adjust_names_for_sbml_export(networkG.metabolites, networkG.actions);
names = [column(metabolites); column(actions)];
x     = networkG.graphics_par.x;
fixed = zeros(length(names),1);
if isfield(networkG.graphics_par,'fixed') && numel(networkG.graphics_par.fixed) == length(names)
  fixed = column(networkG.graphics_par.fixed);
end
invisible = zeros(length(names),1);
if isfield(networkG.graphics_par,'metinvisible') && isfield(networkG.graphics_par,'actinvisible')
  invisible = [column(networkG.graphics_par.metinvisible); column(networkG.graphics_par.actinvisible)];
end

layout_table = sbtab_table_construct(struct('TableID','NetworkLayout','TableType','Position', ...
                 'TableName','Network layout for graphics', 'Comment','Elements not listed are not shown'), ...
                 {'Element','Position:x','Position:y','IsFixed','IsInvisible'}, ...
                 {names, column(x(1,:)), column(x(2,:)), double(fixed(:) ~= 0), double(invisible(:) ~= 0)});
