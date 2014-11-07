function movie_save(filename,M,my_delay)

% movie_save(filename,M,resolution)
%
% argument 'resolution' is currently unused
% delay is in 1/100 seconds

eval(default('my_delay','[]','resolution','30'));

axis tight;
j=1;
for itt = 1:length(M),
  P = frame2im(M(itt));
  imwrite(P,['/tmp/' filename '_' num2str(j) '.bmp'], 'bmp');
  j=j+1;
end

% make movie

disp(['Converting frames to '  filename '.gif']);

string = '! convert ';

if length(my_delay),
  string = [string ' -delay ' num2str(my_delay)]; 
end

for it = 1:j-1,  
  string = [string ' /tmp/' filename '_' num2str(it) '.bmp']; 
end
string = [string ' ' filename '.gif'];
string = [string '; rm /tmp/' filename '_*.bmp'];
eval(string);
