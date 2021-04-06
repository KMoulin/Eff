function [DENSE_param, nodes_DENSE] = Compute_DENSE_nodes(Folder)

DENSE_param=[];

nodes_DENSE=[];

DirZ_vect=[0;0;1];
NPhase=1000;
for cpt_slc=1:1:length(Folder.Dense_files)
    % Loading DENSE structures
    
    DENSE_param(cpt_slc).struct_data   = load(fullfile(Folder.DENSE,Folder.Dense_files(cpt_slc).name)); % DENSE data above

    % DENSE temporal resolution
    DENSE_param(cpt_slc).xres         = DENSE_param(cpt_slc).struct_data.DENSEInfo.PixelSpacing(1); %ms
    DENSE_param(cpt_slc).yres         = DENSE_param(cpt_slc).struct_data.DENSEInfo.PixelSpacing(2); %ms
    DENSE_param(cpt_slc).temporal_res = DENSE_param(cpt_slc).struct_data.SequenceInfo(1, 1).RepetitionTime;
   
    NPhase=min(DENSE_param(cpt_slc).struct_data.DENSEInfo.Number,NPhase);
    DENSE_param(cpt_slc).NPhase= min(DENSE_param(cpt_slc).struct_data.DENSEInfo.Number,NPhase);
    DENSE_param(cpt_slc).RefPhase=1;
    
    DENSE_param(cpt_slc).ImOrien=DENSE_param(cpt_slc).struct_data.SequenceInfo(1,1).ImageOrientationPatient;
    DENSE_param(cpt_slc).dTrans=DENSE_param(cpt_slc).struct_data.SequenceInfo(1, 1).ImagePositionPatient;

    DENSE_param(cpt_slc).Im_1=DENSE_param(cpt_slc).ImOrien(1:3); % attention: this is because matlab read first the colume and second the row
    DENSE_param(cpt_slc).Im_2=DENSE_param(cpt_slc).ImOrien(4:6);     
    DENSE_param(cpt_slc).Im_3=cross(DENSE_param(cpt_slc).Im_1,DENSE_param(cpt_slc).Im_2);

    
    tmp_image.Vect=[];

    tmp_image.Vect(:,1)= DENSE_param(cpt_slc).struct_data.DisplacementInfo.X;
    tmp_image.Vect(:,2)= DENSE_param(cpt_slc).struct_data.DisplacementInfo.Y;
    tmp_image.Vect(:,3)= 0;
    tmp_image.Vect(:,4)= 1;
    
    nodes_DENSE(cpt_slc).dT=DENSE_param(cpt_slc).temporal_res*1e-3; % ms
    
   nodes_DENSE(cpt_slc).Quaternion=[DENSE_param(cpt_slc).Im_1(1)*DENSE_param(cpt_slc).xres DENSE_param(cpt_slc).Im_2(1)*DENSE_param(cpt_slc).yres 0   DENSE_param(cpt_slc).dTrans(1);
                                    DENSE_param(cpt_slc).Im_1(2)*DENSE_param(cpt_slc).xres DENSE_param(cpt_slc).Im_2(2)*DENSE_param(cpt_slc).yres 0   DENSE_param(cpt_slc).dTrans(2);
                                    DENSE_param(cpt_slc).Im_1(3)*DENSE_param(cpt_slc).xres DENSE_param(cpt_slc).Im_2(3)*DENSE_param(cpt_slc).yres 0   DENSE_param(cpt_slc).dTrans(3);
                                    0                        0                        0   1;];
            
   nodes_DENSE(cpt_slc).Quaternion2=[DENSE_param(cpt_slc).Im_1(1)*DENSE_param(cpt_slc).xres DENSE_param(cpt_slc).Im_2(1)*DENSE_param(cpt_slc).yres DENSE_param(cpt_slc).Im_3(1)*DENSE_param(cpt_slc).yres  ;
                                     DENSE_param(cpt_slc).Im_1(2)*DENSE_param(cpt_slc).xres DENSE_param(cpt_slc).Im_2(2)*DENSE_param(cpt_slc).yres DENSE_param(cpt_slc).Im_3(2)*DENSE_param(cpt_slc).yres  ;
                                     DENSE_param(cpt_slc).Im_1(3)*DENSE_param(cpt_slc).xres DENSE_param(cpt_slc).Im_2(3)*DENSE_param(cpt_slc).yres DENSE_param(cpt_slc).Im_3(3)*DENSE_param(cpt_slc).yres  ;];
    
   nodes_DENSE(cpt_slc).Rotation=[DENSE_param(cpt_slc).Im_1(1) DENSE_param(cpt_slc).Im_2(1) DENSE_param(cpt_slc).Im_3(1) ;
                                  DENSE_param(cpt_slc).Im_1(2) DENSE_param(cpt_slc).Im_2(2) DENSE_param(cpt_slc).Im_3(2) ;
                                  DENSE_param(cpt_slc).Im_1(3) DENSE_param(cpt_slc).Im_2(3) DENSE_param(cpt_slc).Im_3(3) ;];
    
    
    tmp_image.DD=[];
    tmp_image.DD(:,1,:)=DENSE_param(cpt_slc).struct_data.DisplacementInfo.dX;
    tmp_image.DD(:,2,:)=DENSE_param(cpt_slc).struct_data.DisplacementInfo.dY;
    tmp_image.DD(:,3,:)=DENSE_param(cpt_slc).struct_data.DisplacementInfo.dZ;
    tmp_image.DD(:,4,:)=0;    
    
    tmp_image.DDrot=[];
    
    % Rotate and scale the displacement field in the new coordinate system,
    % 
    for cpt_t=1:1:size(tmp_image.DD,3)
        tmp_image.DDrot(:,:,cpt_t)= (nodes_DENSE(cpt_slc).Quaternion2(1:3,1:3)*squeeze(tmp_image.DD(:,1:3,cpt_t))')';
    end
    
    tmp_image.ROI_epi=cell2mat(DENSE_param(cpt_slc).struct_data.ROIInfo.Contour(1,1));
    tmp_image.ROI_epi(:,3)=0;
    tmp_image.ROI_epi(:,4)=1;
    tmp_image.ROI_endo=cell2mat(DENSE_param(cpt_slc).struct_data.ROIInfo.Contour(1,2));
    tmp_image.ROI_endo(:,3)=0;
    tmp_image.ROI_endo(:,4)=1;
    
    % Remove the point outside the boundary considering only the first cardiac phase (usefull for LA DENSE) 
    tmp_image.In_out_epi=Collision_ToolBox.ROI_Points(tmp_image.ROI_epi,tmp_image.Vect);    % All the point that are inside the EPI boundary
    tmp_image.In_out_endo=Collision_ToolBox.ROI_Points(tmp_image.ROI_endo,tmp_image.Vect);  % All the point that are inside the ENDO boundary
    tmp_image.In_out=(~tmp_image.In_out_endo)&(tmp_image.In_out_epi);                       % All the point that are inside the EPI and not inside ENDO

    if DENSE_param(cpt_slc).struct_data.ROIInfo.ROIType=='LA'
        nodes_DENSE(cpt_slc).SA=0;
    else
        nodes_DENSE(cpt_slc).SA=1;
    end
    
   
    tmp_image.Vect_Q=(nodes_DENSE(cpt_slc).Quaternion*tmp_image.Vect')';
    
    nodes_DENSE(cpt_slc).DField        = tmp_image.DD(tmp_image.In_out,:,:);
    nodes_DENSE(cpt_slc).DField_rot    = tmp_image.DDrot(tmp_image.In_out,:,:);
    nodes_DENSE(cpt_slc).points(:,:,1) = tmp_image.Vect_Q(tmp_image.In_out,1:3);%;(Quaternion*Vect(In_out,:,1)')';
    
    
    %%% Up to this point everything has been measured from DENSE
    
    % Try to interpolate the long axis slice along the SA Z direction
    if nodes_DENSE(cpt_slc).SA
        nodes_DENSE(cpt_slc).DirZ_vect=nodes_DENSE(cpt_slc).Rotation*[0;0;1]; % projection direction wanted
        nodes_DENSE(cpt_slc).DirZ_vect=nodes_DENSE(cpt_slc).DirZ_vect./norm(nodes_DENSE(cpt_slc).DirZ_vect); % In case the rotation matrix is not scaled, unlikely
    end
    
    for cpt_ph = 1:DENSE_param(cpt_slc).NPhase 
        tmp_image.ROI_epi=cell2mat(DENSE_param(cpt_slc).struct_data.ROIInfo.Contour(cpt_ph,1));
        tmp_image.ROI_epi(:,3)=0;
        tmp_image.ROI_epi(:,4)=1;
        tmp_image.ROI_endo=cell2mat(DENSE_param(cpt_slc).struct_data.ROIInfo.Contour(cpt_ph,2));
        tmp_image.ROI_endo(:,3)=0;
        tmp_image.ROI_endo(:,4)=1;
        nodes_DENSE(cpt_slc).ROI.phase(cpt_ph).epi  = (nodes_DENSE(cpt_slc).Quaternion*tmp_image.ROI_epi')';
        nodes_DENSE(cpt_slc).ROI.phase(cpt_ph).endo = (nodes_DENSE(cpt_slc).Quaternion*tmp_image.ROI_endo')';
    end
end


% The number of phase in the analysis is equal to the min number of phases
for cpt_slc=1:1:length(nodes_DENSE)
    DENSE_param(cpt_slc).NPhase =NPhase;
    nodes_DENSE(cpt_slc).DField=nodes_DENSE(cpt_slc).DField(:,:,1:NPhase);
    nodes_DENSE(cpt_slc).DField_rot=nodes_DENSE(cpt_slc).DField_rot(:,:,1:NPhase);
end


%% Apply the displacement field to the data and Interpolate the Z direction from LA axis to SA (todo interpolate the XY direction of SA to LA)
if strcmp(DENSE_param(1).struct_data.DENSEInfo.Type,'xy')
    [nodes_DENSE]=Eff_Toolbox.Add_Field_LA(nodes_DENSE,DirZ_vect,DENSE_param(cpt_slc).xres);
end
for cpt_ph = 2:DENSE_param(1).NPhase 
    %[Fx Fy Fz]=Eff_Toolbox.Interpolate_Field_LA(nodes_DENSE,cpt_ph-1,DirZ_vect); % Interpolator for the current cardiac phase.
    
    for cpt_slc=1:1:length(nodes_DENSE) % For each DENSE measurement      
             DDproj=nodes_DENSE(cpt_slc).DField_rot(:,:,cpt_ph);  % Take the displacement field in the right coordinate system  
%              if nodes_DENSE(cpt_slc).SA   % Apply displacement field
%                  
%                 tmp_proj=zeros(size(nodes_DENSE(cpt_slc).points(:,1,cpt_ph-1),1),3);
%                 tmp_proj(:,1)= Fx(nodes_DENSE(cpt_slc).points(:,1,cpt_ph-1),nodes_DENSE(cpt_slc).points(:,2,cpt_ph-1),nodes_DENSE(cpt_slc).points(:,3,cpt_ph-1)) ;
%                 tmp_proj(:,2)= Fy(nodes_DENSE(cpt_slc).points(:,1,cpt_ph-1),nodes_DENSE(cpt_slc).points(:,2,cpt_ph-1),nodes_DENSE(cpt_slc).points(:,3,cpt_ph-1)) ;
%                 tmp_proj(:,3)= Fz(nodes_DENSE(cpt_slc).points(:,1,cpt_ph-1),nodes_DENSE(cpt_slc).points(:,2,cpt_ph-1),nodes_DENSE(cpt_slc).points(:,3,cpt_ph-1)) ;
%                 
% 
%                 [Fx Fy Fz]=Eff_Toolbox.Interpolate_Field_LA(nodes_DENSE,cpt_ph-1,DirZ_vect);
%                 DDproj=DDproj+tmp_proj;
%              end
            nodes_DENSE(cpt_slc).points(:,:,cpt_ph) = nodes_DENSE(cpt_slc).points(:,1:3,1)+DDproj; % Apply the full displacement field
    end
end

end
