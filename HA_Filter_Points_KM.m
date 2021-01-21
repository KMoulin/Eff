function [HA_filter]= HA_Filter_Points_KM(HA, Mask_Depth)

    
 for cpt_ph=1:1:size(HA,2)
     
    ListHA=HA(:,cpt_ph);
    ListDist=Mask_Depth(:,cpt_ph);

    ListDist(isnan(ListHA))=nan;
    ListHA(isnan(ListDist))=nan;
    ListDist(isnan(ListDist))=0;
    ListHA(isnan(ListHA))=0;
    
    f = fittype('a*x+b'); 
    fit1 = fit(ListDist,ListHA,f,'StartPoint',[1 1]);
    fdata = feval(fit1,ListDist); 
    I = abs(fdata - ListHA) > 1.5*std(ListHA); 
    outliers = excludedata(ListDist,ListHA,'indices',I);
     
    ListOutliers=ListHA(outliers);
    ListDistOutliers=ListDist(outliers);
    ListOutliers(ListDistOutliers<0.5&ListOutliers<0)=-ListOutliers(ListDistOutliers<0.5&ListOutliers<0);
    ListOutliers(ListDistOutliers>0.5&ListOutliers>0)=-ListOutliers(ListDistOutliers>0.5&ListOutliers>0);
    ListHA(outliers)=ListOutliers;
  
    HA_filter(:,cpt_ph)=ListHA;
 end
 
end
