function [DTI_param, nodes_DTI] = Compute_DTI_nodes(Folder)

load([Folder.DTI 'Interpolation.mat']); 
load([Folder.DTI 'ROI.mat']);                     
load([Folder.DTI 'HA2.mat']);
load([Folder.DTI 'DTI.mat']);  

DTIinfo=dicominfo([enum.dcm_dir '\' enum.dataset.slc.b(1).dir(1).avg(1).filename]);


DTI_param.xres= 0.8;%enum.Pixel(1);
DTI_param.yres= 0.8;%enum.Pixel(2);
DTI_param.Thick=8;%enum.Thickness;
DTI_param.TD=enum.dataset.slc.b(1).dir(1).avg(1).tt; %+enum.TE ceil(cDTI_trigger/temporal_res); % Link to enum => enum.dataset.slc.b(1).dir(1).avg(1).tt  / divided by dense temporal resolution
DTI_param.ImOrien=DTIinfo.ImageOrientationPatient;
DTI_param.dTrans=DTIinfo.ImagePositionPatient;


DTI_param.Im_1=DTI_param.ImOrien(1:3); % attention: this is because matlab read first the colume and second the row
DTI_param.Im_2=DTI_param.ImOrien(4:6);     
DTI_param.Im_3=cross(DTI_param.Im_1,DTI_param.Im_2);
    

nodes_DTI=[];




nodes_DTI.Quaternion=[DTI_param.Im_1(1)*DTI_param.xres DTI_param.Im_2(1)*DTI_param.yres 0   DTI_param.dTrans(1);
                      DTI_param.Im_1(2)*DTI_param.xres DTI_param.Im_2(2)*DTI_param.yres 0   DTI_param.dTrans(2);
                      DTI_param.Im_1(3)*DTI_param.xres DTI_param.Im_2(3)*DTI_param.yres 0   DTI_param.dTrans(3);
                      0                      0                 0   1;];

nodes_DTI.Quaternion2=[DTI_param.Im_1(1)*DTI_param.xres DTI_param.Im_2(1)*DTI_param.yres DTI_param.Im_3(1)*DTI_param.yres  ;
                       DTI_param.Im_1(2)*DTI_param.xres DTI_param.Im_2(2)*DTI_param.yres DTI_param.Im_3(2)*DTI_param.yres  ;
                       DTI_param.Im_1(3)*DTI_param.xres DTI_param.Im_2(3)*DTI_param.yres DTI_param.Im_3(3)*DTI_param.yres  ;];

nodes_DTI.Rotation=[DTI_param.Im_1(1) DTI_param.Im_2(1) DTI_param.Im_3(1) ;
                    DTI_param.Im_1(2) DTI_param.Im_2(2) DTI_param.Im_3(2) ;
                    DTI_param.Im_1(3) DTI_param.Im_2(3) DTI_param.Im_3(3) ;];

Vect_DTI=[];
[Vect_DTI(:,2) Vect_DTI(:,1)]=find(LV_Mask);
Vect_DTI(:,3)=0;
Vect_DTI(:,4)=1;
Vect_DTI_Q=(nodes_DTI.Quaternion*Vect_DTI')';



[tmp_image.cc,tmp_image.rr,tmp_image.ll,tmp_image.wd] = Eff_Toolbox.ComputeDirections(LV_Mask, P_Epi, P_Endo);

% Just check that the fiber are correctly aligned with the circunferential
% direction 
% tmp_dot=dot(squeeze(EigVect1),tmp_image.cc,3);
% for x=1:1:size(EigVect1,1)
%     for y=1:1:size(EigVect1,2)
%         if (tmp_dot(x,y)<0)
%             EigVect1(x,y,1,:)=-EigVect1(x,y,1,:);
%         end
%         if ((tmp_image.wd(x,y)<0.4 & EigVect1(x,y,1,3)<0) | (tmp_image.wd(x,y)>0.6 & EigVect1(x,y,1,3)>0))
%             EigVect1(x,y,1,3)=-EigVect1(x,y,1,3);
%         end
%     end
% end



% Idx=find(dot(nodes_DTI.ff,nodes_DTI.cc,2)<0);
% nodes_DTI.ff(Idx,:)=-nodes_DTI.ff(Idx,:);
Mask_tmp=ones(size(Mask_AHA,1),size(Mask_AHA,2),size(Mask_AHA,3),6);
Mask_tmp=sum(Mask_AHA.*cumsum(Mask_tmp,4),4);
for cpt_v=1:1:size(Vect_DTI,1)
    
    nodes_DTI.points(cpt_v,:)=[Vect_DTI_Q(cpt_v,1),Vect_DTI_Q(cpt_v,2),Vect_DTI_Q(cpt_v,3)]; 
    nodes_DTI.AHA(cpt_v)= Mask_tmp(Vect_DTI(cpt_v,2),Vect_DTI(cpt_v,1));
    nodes_DTI.ff(cpt_v,:)= nodes_DTI.Rotation*squeeze(EigVect1(Vect_DTI(cpt_v,2),Vect_DTI(cpt_v,1),1,:));
    nodes_DTI.ff2(cpt_v,:)= nodes_DTI.Rotation*squeeze(EigVect2(Vect_DTI(cpt_v,2),Vect_DTI(cpt_v,1),1,:));
    nodes_DTI.ff3(cpt_v,:)= nodes_DTI.Rotation*squeeze(EigVect3(Vect_DTI(cpt_v,2),Vect_DTI(cpt_v,1),1,:));
    nodes_DTI.cc(cpt_v,:)=nodes_DTI.Rotation*squeeze(tmp_image.cc(Vect_DTI(cpt_v,2),Vect_DTI(cpt_v,1),:));
    nodes_DTI.rr(cpt_v,:)=nodes_DTI.Rotation*squeeze(tmp_image.rr(Vect_DTI(cpt_v,2),Vect_DTI(cpt_v,1),:));
    nodes_DTI.ll(cpt_v,:)=nodes_DTI.Rotation*squeeze(tmp_image.ll(Vect_DTI(cpt_v,2),Vect_DTI(cpt_v,1),:));
    nodes_DTI.wd(cpt_v,:)=squeeze(tmp_image.wd(Vect_DTI(cpt_v,2),Vect_DTI(cpt_v,1),:));
    
end



% Save the ROI inside the absolute coordinate system.
tmp_image.ROI_epi=[];
tmp_image.ROI_epi(:,1:2)=P_Epi;
tmp_image.ROI_epi(:,3)=0;
tmp_image.ROI_epi(:,4)=1;

tmp_image.ROI_endo=[];
tmp_image.ROI_endo(:,1:2)=P_Endo;
tmp_image.ROI_endo(:,3)=0;
tmp_image.ROI_endo(:,4)=1;
tmp_image.endo=(nodes_DTI.Quaternion*tmp_image.ROI_endo')';
tmp_image.epi=(nodes_DTI.Quaternion*tmp_image.ROI_epi')';

nodes_DTI.ROI.phase(1).epi  = tmp_image.epi(:,1:3);
nodes_DTI.ROI.phase(1).endo = tmp_image.endo(:,1:3);

end