% ----------------------------------------------------------------------
% Create a pink colormap for the current image
% ----------------------------------------------------------------------

function colors( min, max, flipit, jetinstead)

if ~exist('min'), min= 1; end
if ~exist('max'), max= 255; end
if ~exist('flipit'), flipit= 0; end
if ~exist('jet'), jetinstead= 0; end

if flipit==1, 
   h= ones(255,3);  
   h(min:max,:)= flipud(pink(max-min+1));
   for i=max:255, h(i,:)= h(max,:); end
else 
   h= zeros(255,3); 
   h(min:max,:)= pink(max-min+1);
   for i=max:255, h(i,:)= h(max,:); end
end

if jetinstead==1,
   h= zeros(255,3);
   dmy= jet(1.15*(max-min+1));
   h(min:max,:)= dmy(1:max-min+1,:);
   for i=1:min, h(i,:)= h(min,:); end
   for i=max:255, h(i,:)= h(max,:); end
end

colormap(h)

% show colormap in a seperate figure
%dmy= gcf
%figure
%image(-100:0,0:0,1:255)
%colormap(h)
%figure(dmy)
