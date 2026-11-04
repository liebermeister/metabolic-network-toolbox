% NETGRAPH_ASSIGN_POSITIONS Choose metabolite positions on a picture with the mouse
%
% [names,xpos] = netgraph_assign_positions(network,picture,offsets,outfile)
%
% For each metabolite, the program asks (on the command line) whether to choose a position
% (by a mouse click on the picture), go forward, go back, or finish.
%
% picture  image matrix (as for 'image')
% offsets  [x,y] added to the positions (default [0,0])
% outfile  file to which positions are written (optional)

function [names,xpos] = netgraph_assign_positions(network,picture,offsets,outfile)

eval(default('offsets','[0,0]','outfile','[]'));

clf;
subplot('position',[0 0 1 1]);
image(picture);
axis equal; axis tight;
title('See command line');

names = network.metabolites;
ntot = length(names);

if isfield(network,'graphics_par')
  xpos = network.graphics_par.x(:,1:ntot) + repmat(offsets',1,ntot);
else
  xpos = nan * ones(2,ntot);
end

counter = 1;

while counter <= ntot

  display(sprintf('%s %f %f',names{counter},xpos(1,counter),xpos(2,counter)));
  r = input(sprintf('[1] choose position, [2] go forward, [3] go back, [4] finish:   '));
  switch r

    case 1
      [x,y,button] = ginput(1);
      xpos(:,counter) = [x,y] + offsets;
      display(sprintf(' %f %f',x,y));  counter = counter + 1;

    case 2, counter = counter + 1;
      display([' skipped']);

    case 3
      counter = max(1,counter - 1);
      display([' going back to previous metabolite']);

    case 4, counter = ntot + 1;

  end

end


if ~isempty(outfile)
  display(sprintf('Writing positions to file %s',outfile));
  mytable([names,cellstr(num2str(xpos(1,:)')),cellstr(num2str(xpos(2,:)'))],0,outfile);
end
