function colors = netgraph_movie_background_colors(colors, use_colors, n_frames)

% colors = netgraph_movie_background_colors(colors, use_colors, n_frames)
%
% Background colours for the frames of a network movie (used by the netgraph_..._movie functions)
%
% colors      matrix of RGB colours (one row per colour, or one column per colour);
%             the colours are interpolated to obtain one colour per frame
% use_colors  flag; if 0, an empty matrix is returned
% n_frames    number of frames
%
% Returns a (n_frames x 3) matrix, or [] if no background colours are used

if ~use_colors || isempty(colors)
  colors = []; return
end

if size(colors,2) ~= 3
  colors = colors';
end

nb = size(colors,1);
if nb == 1
  colors = repmat(colors, n_frames, 1);
elseif n_frames == 1
  colors = colors(1,:);
else
  colors = interp1(linspace(0,1,nb), colors, linspace(0,1,n_frames));
end
