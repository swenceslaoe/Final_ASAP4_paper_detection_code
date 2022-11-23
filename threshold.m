function thr = threshold(F)

nF =(F-min(F))/(max(F)-min(F));
thr = graythresh(nF);
thr = thr*(max(F)-min(F))+min(F);
