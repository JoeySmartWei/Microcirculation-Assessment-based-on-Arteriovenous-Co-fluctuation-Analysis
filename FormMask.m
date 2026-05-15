function [RatROIs] = FormMask(CBFimg,VesslesType,ROINums)

% Vessels Position
[wid len channels] = size(CBFimg);
figure(1);set(gcf,'position',[150 150 1300 850]);
imagesc(CBFimg);axis image;
title(['Please Choose ',VesslesType,' ROIs Points'],'FontSize',30);
[x_vessel,y_vessel] = ginput(ROINums);close;
Location = round([x_vessel y_vessel]);

Regs = size(Location,1);
RatROIs = false(wid,len,Regs+1);
for reg = 1:Regs
    
    % 指定圆心位置、半径
    Rmax = 3;
    x0 = Location(reg,1);
    y0 = Location(reg,2);
    theta = 0:pi/100:2*pi;
    % 生成圆坐标位置
    xPos = x0 + Rmax * cos (theta);
    yPos = y0 + Rmax * sin (theta);
    % 创造掩膜
    BW = false(wid,len);
    addedRegion = poly2mask(xPos, yPos, wid, len);
    Inx = find(addedRegion == 1);
    addedRegion(Inx(end)) = 0;
    BW = BW | addedRegion;
    RatROIs(:,:,reg) = BW;
end

for reg = 1:Regs
    RatROIs(:,:,Regs+1) = RatROIs(:,:,Regs+1) | RatROIs(:,:,reg);
end