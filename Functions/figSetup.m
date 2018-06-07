% setup figure for normal printing layout
set(gcf,'Color','w');
set(gca,'FontName','Calibri','FontSize',18);
fHandles = findall(gcf,'type','line');
for ff = fHandles
    set(ff,'LineWidth',2);
end