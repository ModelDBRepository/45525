function pmovie()
global RAW2D 
[m,n,k]= size(RAW2D)

lo= -1.25; hi= -0.5;
S= (RAW2D<lo)*lo + (RAW2D>hi)*hi + ((RAW2D>=lo)&(RAW2D<=hi)).*RAW2D - lo;
S= round(4*S+1);

%colormap(flipud(pink(4*(hi-lo)+1)))
colors(1,128,1,1)

subplot(1,1,1);
for i=1:k
   image(30*S(:,:,i))
   set(gcf,'DefaultTextFontSize',10)
   ss=i*25;
   pp1=num2str(ss);
   pp2=['T=' pp1 ' '];
   text(3,-3,pp2)
   ii=num2str(i);
   rr=['mov' ii '.png']
   print('-dpng',rr)
end

%colormap(flipud(hot(4*(hi-lo)+1)))
%colormap(flipud(pink(4*(hi-lo)+1)))
%aa=axes('position',[0.92 0.33 0.02 0.2])
%colorbar(aa)
%colors(1,128,1,1)
%set(aa,'YLim',[50 100])
%set(aa,'YTicklabel',{'0';' ';'1'})

