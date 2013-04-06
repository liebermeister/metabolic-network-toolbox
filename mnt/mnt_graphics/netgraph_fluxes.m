% netgraph_fluxes(network,J)

function netgraph_fluxes(network,J)

netgraph_draw(network,struct('actvalues',J,'arrowvalues',J,'arrowstyle','fluxes'));