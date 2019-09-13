% Load saved figures
c=hgload('no_diffusion.fig');
k=hgload('with_diffusion.fig');
% Prepare subplots
figure
h(1)=subplot(1,2,1);
h(2)=subplot(1,2,2);
% Paste figures on the subplots
copyobj(allchild(get(c,'CurrentAxes')),h(1));
copyobj(allchild(get(k,'CurrentAxes')),h(2));
% Add legends
l(1)=legend(h(1),'No Diffusion')
l(2)=legend(h(2),'With Diffusion')