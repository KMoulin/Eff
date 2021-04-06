classdef Eff_Toolbox

    methods(Static)
                 %% Generate the Depth Mask Map from ROI
                 function Mask_Depth = WallDepth(LV_mask,P_Epi,P_Endo)
                        [Xq,Yq] = meshgrid(1:size(LV_mask,2),1:size(LV_mask,1));
                        Endo_Line = zeros(size(P_Endo));
                        Epi_Line = ones(size(P_Epi));
                        PosRoi = cat(1,P_Epi,P_Endo);
                        LineRoi   = cat(1,Epi_Line,Endo_Line);
                        Mask_Depth = griddata(PosRoi(:,1),PosRoi(:,2),LineRoi(:,1),Xq,Yq);
                 end
                 
                 %% Generate circunferiential,radial and longitudinal direction and depth mask from ROI
                 function [cc,rr,ll,wd] = ComputeDirections(LV_mask, P_Epi, P_Endo)
                     yres = size(LV_mask,1);   xres = size(LV_mask,2);
                    [Xq,Yq] = meshgrid(1:xres,1:yres);

                    npts = size(P_Epi,1);
                    npts2 = size(P_Endo,1);
                    
                    cc = zeros(yres,xres,3);
                    rr = zeros(yres,xres,3);
                    ll = zeros(yres,xres,3);

                    [row,col] = find(LV_mask);
                    center = mean([row,col]);


                    for z=1:size(LV_mask,3)

                        Vec = zeros(npts,2);
                        Vec2 = zeros(npts2,2);
                        positions=[];
                        vectors=[];


                         Vec(1,:) = P_Epi(1,:,z) - P_Epi(end,:,z);
                        for y = 2:npts
                            Vec(y,:) = P_Epi(y,:,z) - P_Epi(y-1,:,z);
                        end

                        Vec2(1,:) = P_Endo(1,:,z) - P_Endo(end,:,z);
                        for y = 2:npts2
                            Vec2(y,:) = P_Endo(y,:,z) - P_Endo(y-1,:,z);
                        end

                        positions = cat(1,P_Epi(:,:,z),P_Endo(:,:,z));
                        vectors   = cat(1,Vec,Vec2);



                        Vy = griddata(positions(:,1),positions(:,2),vectors(:,2),Xq,Yq);
                        Vx = griddata(positions(:,1),positions(:,2),vectors(:,1),Xq,Yq);


                        for y = 1:yres
                            for x = 1:xres
                                if LV_mask(y,x,z) ~= 0

                                    cc(y,x,:) = [Vx(y,x) Vy(y,x) 0];
                                    ll(y,x,:) = [0 0 1];
                                    rvec = [y,x] - center;
                    %                 rvec = rvec/norm(rvec);
                                    rr(y,x,:) = [rvec(2) rvec(1) 0];%cross(cc(y,x,:),ll(y,x,:));

                                end
                            end
                        end
                    end
                    
                    wd = Eff_Toolbox.WallDepth(LV_mask,P_Epi,P_Endo);
                    
                 end
                 
                 %% Query Nb_Point per Points_query in a Points_List
                 function IDXQuery = Query(Points_List,Points_query,Nb_Point)                   
                     
                     for cpt=1:1:size(Points_query,1)
                         dd=[];
                        for cpt_dim=1:1:size(Points_query,2)                 
                            dd(:,cpt_dim)=Points_query(cpt,cpt_dim)-Points_List(:,cpt_dim);
                        end
                        Euclid_Distance=sum(dd.^2,2);
                        [~,indx]= sort(Euclid_Distance);

                        IDXQuery(cpt,:)=indx(1:Nb_Point);
                        %IDXQuery(cpt,:)=k(1:min(Nb_Point,size(k,1)));
                     end
                 end
                  
                 function [nodes_reg]=Register_nodes(nodes_ref,nodes_moving,Phase,Dir_vect)
                        
                        InterpPoint=[];
                        InterpD=[];
                        Vect_proj_nodes=[];
                        nodes_reg=nodes_moving;
                        
                        center_ref=mean(nodes_ref.points(:,:,Phase));
                        center_mov=mean(nodes_reg.points(:,:,1));
                        
                        center_diff=center_ref-center_mov;
                        for cpt_direction=1:1:size(Dir_vect,2) 
                                Shift= (Vector_ToolBox.Projection_vect(center_diff, nodes_ref.Rotation*Dir_vect(:,cpt_direction))); 
                                
                                % Shift the node points
                                for cpt_point=1:1:size(nodes_reg.points,1)
                                    nodes_reg.points(cpt_point,:)=nodes_reg.points(cpt_point,:)+Shift';
                                end
                                
                                % Shift the ROI epi points
                                for cpt_point=1:1:size(nodes_reg.ROI.phase(1).epi,1)
                                    nodes_reg.ROI.phase(1).epi(cpt_point,:)=nodes_reg.ROI.phase(1).epi(cpt_point,:)+Shift';                             
                                end
                                
                                % Shift the ROI endo points
                                for cpt_point=1:1:size(nodes_reg.ROI.phase(1).endo,1)
                                    nodes_reg.ROI.phase(1).endo(cpt_point,:)=nodes_reg.ROI.phase(1).endo(cpt_point,:)+Shift';                             
                                end
                        end
                 end
                 
                 
                 function [nodes_reg]=Register_nodes3D(nodes_ref,nodes_moving,Phase,Dir_vect)
                        
                       [nodes_reg1]=Eff_Toolbox.Register_nodes(nodes_ref(1),nodes_moving,Phase,Dir_vect);
                       [nodes_reg2]=Eff_Toolbox.Register_nodes(nodes_ref(2),nodes_moving,Phase,Dir_vect);
                        nodes_reg=nodes_reg2;
                        nodes_reg.points=(nodes_reg2.points+nodes_reg1.points)/2;
                  end
                 
                 function [nodes_reg]=Register_optim_nodes(nodes_ref,nodes_moving,Phase)
                        
                        nodes_reg=nodes_moving;
                        

                        point=(inv(nodes_moving.Rotation)*nodes_moving.points(:,:,1)')';
                        point_d=(inv(nodes_moving.Rotation)*nodes_ref.points(:,:,Phase)')';
                        
                        point_endo= (inv(nodes_moving.Rotation)*nodes_ref.ROI.phase(Phase).endo(:,1:3)')';
                        point_epi = (inv(nodes_moving.Rotation)*nodes_ref.ROI.phase(Phase).epi(:,1:3)')';
                        
                        Col_endo=Collision_ToolBox.ROI_Points(point_endo,point);
                        Col_epi=Collision_ToolBox.ROI_Points(point_epi,point);
                        
                        Col=(~Col_endo)&(Col_epi);
                       
                        tmp_error=[];
                        Error=sum(Col);
                        Shift=[0 0 0];
                        cpt=1;
                        for cpt_x=-40:1:40
                            for cpt_y=-40:1:40

                                  tmp_point(:,1)=point(:,1)+cpt_x;
                                  tmp_point(:,2)=point(:,2)+cpt_y;
                                  
                                  
                                  Col_endo=Collision_ToolBox.ROI_Points(point_endo,tmp_point);
                                  Col_epi=Collision_ToolBox.ROI_Points(point_epi,tmp_point);
                        
                                  Col=(~Col_endo)&(Col_epi);
                                   cpt=cpt+1;  
                                  tmp_error(cpt)=sum(Col);
                                  
                                  if tmp_error(cpt)>Error
                                      Error=tmp_error(cpt);
                                      Shift=[cpt_x cpt_y 0];
                                  end
                                  cpt=cpt+1;
                                      
                            end
                        end
                        

                        Shift_rot=  nodes_moving.Rotation*Shift'; 

                        % Shift the node points
                        for cpt_point=1:1:size(nodes_reg.points,1)
                            nodes_reg.points(cpt_point,:)=nodes_reg.points(cpt_point,:)+Shift_rot';
                        end

                        % Shift the ROI epi points
                        for cpt_point=1:1:size(nodes_reg.ROI.phase(1).epi,1)
                            nodes_reg.ROI.phase(1).epi(cpt_point,:)=nodes_reg.ROI.phase(1).epi(cpt_point,:)+Shift_rot';                             
                        end

                        % Shift the ROI endo points
                        for cpt_point=1:1:size(nodes_reg.ROI.phase(1).endo,1)
                            nodes_reg.ROI.phase(1).endo(cpt_point,:)=nodes_reg.ROI.phase(1).endo(cpt_point,:)+Shift_rot';                             
                        end

                 end
                 
                 function [Fx Fy Fz]=Interpolate_Field_LA(nodes,cpt_ph,Dir_vect)
                        InterpPoint=[];
                        InterpD=[];
                        Vect_proj_nodes=[];
                        for cpt_slc=1:1:size(nodes,2) 
                            if nodes(cpt_slc).SA 
                                Vect_proj_nodes=nodes(cpt_slc).Rotation*Dir_vect; %% We assume the projection direction in SA is the same. 
                            end
                        end
                        
                        
                        for cpt_slc=1:1:size(nodes,2) 
                            if ~nodes(cpt_slc).SA % If it's Long Axis View, Add it to the pool of interpolation       
                                
                                %% Add option if LA <3 just take the closest points and apply it to the SA
                                    InterpPoint=[InterpPoint; nodes(cpt_slc).points(:,:,cpt_ph)];
                                    InterpD=[InterpD; (Vector_ToolBox.Projection_vect_n(nodes(cpt_slc).DField_rot(:,:,cpt_ph),Vect_proj_nodes))];            
                            end  
                        end

                        Fx = scatteredInterpolant(InterpPoint(:,1),InterpPoint(:,2),InterpPoint(:,3),InterpD(:,1),'natural','linear');   
                        Fy = scatteredInterpolant(InterpPoint(:,1),InterpPoint(:,2),InterpPoint(:,3),InterpD(:,2),'natural','linear');
                        Fz = scatteredInterpolant(InterpPoint(:,1),InterpPoint(:,2),InterpPoint(:,3),InterpD(:,3),'natural','linear');
                 end

                 function [nodes]=Add_Field_LA(nodes,Dir_vect,Pixel_cutoff)
                        LA_Points=[];
                        LA_D=[];
                        Vect_proj_nodes=[];

                        N_max=0;
                        % First Check the dataset
                        for cpt_slc=1:1:size(nodes,2) 
                            if nodes(cpt_slc).SA 
                                Vect_proj_nodes=nodes(cpt_slc).Rotation*Dir_vect; %% We assume the projection direction in SA is the same.
                            end
                        end
                        
                        % Make sure the dataset is homogenenous and copy
                        % the diplacement data from long axis
                        for cpt_slc=1:1:size(nodes,2)                            
                             if ~nodes(cpt_slc).SA % If it's Long Axis View, Add it to the pool of interpolation       
                                
                                    %% Add option if LA <3 just take the closest points and apply it to the SA
                                    LA_Points=[LA_Points; nodes(cpt_slc).points(:,:,1)];
                                    LA_D=[LA_D; (Vector_ToolBox.Projection_vect_n_t(nodes(cpt_slc).DField_rot(:,:,:),Vect_proj_nodes))];   
                                    
                            end  
                        end
        
                        
                        % Add the long axis displacement data to the short axis one.  
                        for cpt_slc=1:1:size(nodes,2)                            
                            if nodes(cpt_slc).SA % If it's SA Axis View, Look for the closest nodes from the LA pool and add the mean diplacement field. 
                                [Idx Dist]= dsearchn(LA_Points,nodes(cpt_slc).points(:,:,1));
                                Idx2=Idx(Dist<2*Pixel_cutoff); 
                                
                                % Apply the mean displacement field, it
                                % should be probably better to came up with
                                % some interpolation method instead. 
                                
                                
                                nodes(cpt_slc).DField_rot=nodes(cpt_slc).DField_rot+repmat((mean(LA_D(Idx2,:,:))),size(nodes(cpt_slc).DField_rot,1),1,1);
                            end 
                        end    
                 end
                 
                 
                 
                 %% View
                 
                 
                 
                  %% Display points color coded by the parameter value
                 function Scatter3_value(Points,Value,Range)
                     
                     JC=jet(256);
                     
                     
                     Value(isnan(Value))=0;
                     Value=(Value - Range(1))/( Range(2)-Range(1));
                     Value(Value<0)=0;
                     Value(Value>1)=1;
                     
                     
                     h=figure;
                     set(gca,'Clipping','off');
                     for cpt_t=1:1:size(Points,3)
                        clf(h)
                        hold on
                        
                        scatter3(Points(:,1,cpt_t),Points(:,2,cpt_t),Points(:,3,cpt_t),20,JC(round(Value(:,cpt_t)*255)+1,:),'filled')

                        axis([-60 60 -60 60 -60 60]);
                       
                        view(-159, -24)
                        zoom(2)
                        pause(0.1)
                    end
                 end
                    
                 function Scatter3_vector(Points,Vector)
                                     
                     
                     Vector(isnan(Vector))=0;
                     
                     
                     h=figure;
                     
                     for cpt_t=1:1:size(Points,3)
                        clf(h)
                        hold on
                        
                        scatter3(Points(:,1,cpt_t),Points(:,2,cpt_t),Points(:,3,cpt_t),20,Vector(:,:,cpt_t),'filled')

                        axis([-60 60 -60 60 -60 60]);
                       
                        view(-159, -24)
                        zoom(2)
                        set(gca,'Clipping','off');
                        pause(0.1)
                    end
                 end
                 
                 function Scatter3_vector_rot(Points,Vector,Mat_rot)
                                     
                     
                     Vector(isnan(Vector))=0;
                     
                     
                     h=figure;
                     set(gca,'Clipping','off');
                     for cpt_t=1:1:size(Points,3)
                        clf(h)
                        hold on
                        
                        scatter3(Points(:,1,cpt_t),Points(:,2,cpt_t),Points(:,3,cpt_t),20,abs((inv(Mat_rot)*squeeze(Vector(:,:,cpt_t))')'),'filled')

                        axis([-60 60 -60 60 -60 60]);
                        
                        view(-159, -24)
                        set(gca,'Clipping','off');
                        zoom(2)
                        pause(0.1)
                    end
                 end
                 
                  %% Display points color coded by the parameter value
                 function Quiver3_vector(Points,Vector)
                     
                     JC=jet(256);
                     
                     h=figure;
                     
                     for cpt_t=1:1:size(Points,3)
                        clf(h)
                        hold on
                        
                        
                       
                        quiver3(Points(:,1,cpt_t),Points(:,2,cpt_t),Points(:,3,cpt_t),Vector(:,1,cpt_t),Vector(:,2,cpt_t),Vector(:,3,cpt_t),0.5)

                        axis([-60 60 -60 60 -60 60]);
                       
                        view(-159, -24)
                        set(gca,'Clipping','off');
                         zoom(2)
                        pause(0.1)
                    end
                 end
                    
                  %% Display points color coded by the parameter value
                 function Scatter4D(Points,VectorColor,RotMat,Range)
                     
                        if nargin<4
                            Range(1)=nanmin(nanmin(nanmin(VectorColor)));
                            Range(2)=nanmax(nanmax(nanmax(VectorColor)));
                        end
                        if nargin<3
                            RotMat=[1 0 0; 0 1 0; 0 0 1];
                        end
                        
                        JC=hot(256);
                        if(size(VectorColor,2)==3) % Vector plot
                            colorMap=abs(inv(RotMat)*VectorColor(:,:,1)')';
                        else                  % Scallar plot   
                             VectorColor(isnan(VectorColor))=0;   
                             VectorColor=(VectorColor - Range(1))/( Range(2)-Range(1));
                             VectorColor(VectorColor<0)=0;
                             VectorColor(VectorColor>1)=1;
                             colorMap=JC(round(VectorColor(:,1)*255)+1,:);
                        end
                        
                        Max_X=max(max(Points(:,1,:)));
                        Min_X=min(min(Points(:,1,:)));
                        
                        
                        Max_Y=max(max(Points(:,2,:)));
                        Min_Y=min(min(Points(:,2,:)));
                        
                        Max_Z=max(max(Points(:,3,:)));
                        Min_Z=min(min(Points(:,3,:)));
                        
                        FigH = figure('Units', 'normal', 'Position', [0.1 0.1 .8 .8]);  %not quite full screen  %[left bottom width height]

                        subgroup_im = axes('Parent', FigH, 'Units', 'normal', 'Position', [0 1/10 1 9/10]);       %3/4 of the plot
                        subgroup_button = uipanel('Parent', FigH, 'Units', 'normal', 'Position', [0 0 1 1/10]);  %1/10 of the plot
                        
                        SliderTime =  uicontrol('Style','slider'                      ,'Parent', subgroup_button,  'Units', 'normal','position', [0 0.1 0.5 0.3],'Tag','slider1','min', 1, 'max', size(Points,3),'Value',1,'SliderStep', [1/size(Points,3) , 1/size(Points,3)],'Callback',@callback_Refresh ); %,'Position',[180 260 200 20]
                        SliderScale = uicontrol('Style','slider' ,'String', 'scale'   ,'Parent', subgroup_button,  'Units', 'normal','position', [0.6 0.1 0.3 0.3],'Tag','slider1','min', 1, 'max',10,'Value',1,'SliderStep', [1/10, 1/10],'Callback',@callback_Scale ); 

                        TextTime =    uicontrol('style','text','String', 'Time Frame' ,'Parent', subgroup_button,  'Units', 'normal','position', [0 0.5 0.1 0.3]); 
       
                        shandle=scatter3(Points(:,1,1),Points(:,2,1),Points(:,3,1),20,colorMap,'filled','Parent',subgroup_im);
                        
                        
                        axis([Min_X Max_X Min_Y Max_Y Min_Z-0.001 Max_Z]);
                        view(-159, -24)
                        set(gca,'Clipping','off');
                        
                     function callback_Refresh(source, eventdata)
                        cpt_time = round(source.Value);   
                         if(size(VectorColor,2)==3) % Vector plot
                           colorMap=abs(inv(RotMat)*VectorColor(:,:,cpt_time)')';
                         else  
                            colorMap=JC(round(VectorColor(:,cpt_time)*255)+1,:);
                         end
                        TextTime.String = ['Time Frame' num2str(cpt_time) ];
                        shandle.XData=Points(:,1,cpt_time);
                        shandle.YData=Points(:,2,cpt_time);
                        shandle.ZData=Points(:,3,cpt_time);
                        shandle.CData=colorMap;
                     end
                     function callback_Scale(source, eventdata)                        
                        shandle.SizeData= 20*round(source.Value);                      
                     end
                 end
                 
                  %% Display points color coded by the parameter value
                 function Quiver4D(Points,Vector)
                    
                        % create subhandles
                        
                        Max_X=max(Points(:,1,1));
                        Min_X=min(Points(:,1,1));
                        
                        
                        Max_Y=max(Points(:,2,1));
                        Min_Y=min(Points(:,2,1));
                        
                        Max_Z=max(Points(:,3,1));
                        Min_Z=min(Points(:,3,1));
                        
                        FigH = figure('Units', 'normal', 'Position', [0.1 0.1 .8 .8]);  %not quite full screen  %[left bottom width height]

                        subgroup_im = axes('Parent', FigH, 'Units', 'normal', 'Position', [0 1/10 1 9/10]);       %3/4 of the plot
                        subgroup_button = uipanel('Parent', FigH, 'Units', 'normal', 'Position', [0 0 1 1/10]);  %1/10 of the plot
                        
                        SliderTime =  uicontrol('Style','slider'                      ,'Parent', subgroup_button,  'Units', 'normal','position', [0 0.1 0.5 0.3],'Tag','slider1','min', 1, 'max', size(Points,3),'Value',1,'SliderStep', [1/size(Points,3) , 1/size(Points,3)],'Callback',@callback_Refresh ); %,'Position',[180 260 200 20]
                        SliderScale = uicontrol('Style','slider' ,'String', 'scale'   ,'Parent', subgroup_button,  'Units', 'normal','position', [0.6 0.1 0.3 0.3],'Tag','slider1','min', 1, 'max',10,'Value',1,'SliderStep', [1/10, 1/10],'Callback',@callback_Scale ); 

                        TextTime =    uicontrol('style','text','String', 'Time Frame' ,'Parent', subgroup_button,  'Units', 'normal','position', [0 0.5 0.1 0.3]); 
       
                        qhandle=quiver3(Points(:,1,1),Points(:,2,1),Points(:,3,1),Vector(:,1,1),Vector(:,2,1),Vector(:,3,1),'Parent',subgroup_im);
                        set(qhandle,'Parent',subgroup_im);
                        axis([Min_X Max_X Min_Y Max_Y Min_Z Max_Z]);
                        view(-159, -24)
                        set(gca,'Clipping','off');
                        
                     function callback_Refresh(source, eventdata)
                        cpt_time = round(source.Value);   
                        
                        TextTime.String = ['Time Frame' num2str(cpt_time) ];
                        qhandle.XData=Points(:,1,cpt_time);
                        qhandle.YData=Points(:,2,cpt_time);
                        qhandle.ZData=Points(:,3,cpt_time);
                        qhandle.UData=Vector(:,1,cpt_time);
                        qhandle.VData=Vector(:,2,cpt_time);
                        qhandle.VData=Vector(:,3,cpt_time);
                     end
                      function callback_Scale(source, eventdata)                        
                        qhandle.AutoScaleFactor= round(source.Value);
                        
                     end
                 end
                  %% Display points color coded by the parameter value
                 function Cine4D(nodes)  
                       
                        Points=[];
                        for cpt_node=1:1:size(nodes,2)
                            Points=[Points ; nodes(cpt_node).points];
                        end
                        
                        Max_X=max(Points(:,1,1));
                        Min_X=min(Points(:,1,1));
                        
                        
                        Max_Y=max(Points(:,2,1));
                        Min_Y=min(Points(:,2,1));
                        
                        Max_Z=max(Points(:,3,1));
                        Min_Z=min(Points(:,3,1));
                        
                        FigH = figure('Units', 'normal', 'Position', [0.1 0.1 .8 .8]);  %not quite full screen  %[left bottom width height]

                        subgroup_im = axes('Parent', FigH, 'Units', 'normal', 'Position', [0 1/10 1 9/10]);       %3/4 of the plot
                        subgroup_button = uipanel('Parent', FigH, 'Units', 'normal', 'Position', [0 0 1 1/10]);  %1/10 of the plot
                        
                        SliderTime =  uicontrol('Style','slider'                      ,'Parent', subgroup_button,  'Units', 'normal','position', [0 0.1 0.5 0.3],'Tag','slider1','min', 1, 'max', size(Points,3),'Value',1,'SliderStep', [1/size(Points,3) , 1/size(Points,3)],'Callback',@callback_Refresh ); %,'Position',[180 260 200 20]
                        SliderScale = uicontrol('Style','slider' ,'String', 'scale'   ,'Parent', subgroup_button,  'Units', 'normal','position', [0.6 0.1 0.3 0.3],'Tag','slider1','min', 1, 'max',10,'Value',1,'SliderStep', [1/10, 1/10],'Callback',@callback_Scale ); 

                        TextTime =    uicontrol('style','text','String', 'Time Frame' ,'Parent', subgroup_button,  'Units', 'normal','position', [0 0.5 0.1 0.3]); 
       
                        
                        shandle=plot3(Points(:,1,1),Points(:,2,1),Points(:,3,1),'bo','Parent',subgroup_im);
                        
                        
                        axis([Min_X Max_X Min_Y Max_Y Min_Z Max_Z]);
                        view(-159, -24)
                        set(gca,'Clipping','off');
                        
                     function callback_Refresh(source, eventdata)
                        cpt_time = round(source.Value);   

                        TextTime.String = ['Time Frame' num2str(cpt_time) ];
                        shandle.XData=Points(:,1,cpt_time);
                        shandle.YData=Points(:,2,cpt_time);
                        shandle.ZData=Points(:,3,cpt_time);
                     end
                     function callback_Scale(source, eventdata)                        
                        shandle.SizeData= 20*round(source.Value);                      
                     end
                 end
                 
                  %% Display points color coded by the parameter value
                 function Strain_plot(Strain)
                    h=figure;
                    plot(nanmedian(Strain))
                    hold on
                    
                    plot(quantile(Strain,0.25),'--')
                    plot(quantile(Strain,0.75),'--')
                    %plot(nanmean(Strain)+nanstd(Strain),'--')
                    %plot(nanmean(Strain)-nanstd(Strain),'--')
                 end
                 
                  %% Display points color coded by the parameter value
                 function Strain_Wd_plot(Strain,Walldepth,Range)
                    
                 for cpt_t=1:1:size(Strain,2)
                    Endo_list=find(Walldepth(:,cpt_t)<0.3);
                    Mid_list=find(Walldepth(:,cpt_t)>=0.3&Walldepth(:,cpt_t)<=0.7);
                    Epi_list=find(Walldepth(:,cpt_t)>0.7);
                    Mean_Endo(cpt_t)=nanmedian(Strain(Endo_list,cpt_t));
                    Q_Endo(cpt_t,:)=quantile(Strain(Endo_list,cpt_t),[0.25 0.75]);
                   
                    Mean_Mid(cpt_t)=nanmedian(Strain(Mid_list,cpt_t));
                    Q_Mid(cpt_t,:)=quantile(Strain(Mid_list,cpt_t),[0.25 0.75]);
                    
                    Mean_Epi(cpt_t)=nanmedian(Strain(Epi_list,cpt_t));
                    Q_Epi(cpt_t,:)=quantile(Strain(Epi_list,cpt_t),[0.25 0.75]);
                 end  
                  h=figure;
                  hold on 
                    
                    plot(Mean_Endo,'r','Linewidth',3)
                    plot(Q_Endo(:,1),'r--','Linewidth',3)
                    plot(Q_Endo(:,2),'r--','Linewidth',3)
                    
                    plot(Mean_Mid,'g','Linewidth',3)
                    plot(Q_Mid(:,1),'g--','Linewidth',3)
                    plot(Q_Mid(:,2),'g--','Linewidth',3)
                    
                    plot(Mean_Epi,'b','Linewidth',3)
                    plot(Q_Epi(:,1),'b--','Linewidth',3)
                    plot(Q_Epi(:,2),'b--','Linewidth',3)
                    
                    ax = gca;
                    ax.FontSize=15;
                    ax.FontWeight='Bold';
                    ax.LineWidth=3;
                    
                    box off
                    set(gcf,'color','w');
                    
                    if nargin>2
                        axis([0 size(Strain,2) Range(1) Range(2)])
                    end
                 end
                 function Strain_AHA_plot(Strain,AHA,Range)
                 
                 for cpt_AHA=1:1:6    
                     for cpt_t=1:1:size(Strain,2)

                        Mean_AHA(cpt_t,cpt_AHA)=nanmedian(Strain(find(AHA==cpt_AHA),cpt_t));
                        Q_AHA(cpt_t,cpt_AHA,:)=quantile(Strain(find(AHA==cpt_AHA),cpt_t),[0.25 0.75]);

                     end
                 end
                    h=figure;
                    hold on 
                    plot(Mean_AHA,'-')
                    plot(squeeze(Q_AHA(:,:,1)),'--')
                    plot(squeeze(Q_AHA(:,:,2)),'--')
                 
                    if nargin>2
                        axis([0 size(Strain,2) Range(1) Range(2)])
                    end
                 end
                 function  HAR=HA_WD(wd,HA)
                     
                     HAR=[];
                     % PLOT HA 
                    h=figure('Units', 'normal','position', [0 0 1 1]) 
                    for cpt_t=1:1:size(wd,2) %DTI_param.CardiacPhase %
                        clf(h)
                        hold on
                        
                        
                        scatter(wd(:,cpt_t),HA(:,cpt_t),100,'r','filled')
                        
                         count = 1;
                         Med=[];
                         wd_division=linspace(0,1,19);
                        for cpt_div=1:1:19-1 
                            I = find(wd(:,cpt_t)> wd_division(cpt_div) &  wd(:,cpt_t)<= wd_division(cpt_div+1) );
                            tmp=HA(I,cpt_t);
                            tmp=tmp(~isnan(tmp));
                            Med(cpt_div) = median(tmp);
                            Q(cpt_div,:) = quantile(tmp,[0.25 0.75]);
                            plot(wd_division(cpt_div)+0.05,Med(cpt_div),'k+', 'MarkerSize',30, 'LineWidth',6)
                            plot(wd_division(cpt_div)+0.05,Q(cpt_div,1),'ks', 'MarkerSize',30, 'LineWidth',6)
                            plot(wd_division(cpt_div)+0.05,Q(cpt_div,2),'ks', 'MarkerSize',30, 'LineWidth',6)
                        end
                        
                        HAR(cpt_t)=Med(1)-Med(end);
                        
                        
                        axis([0 1 -90 90 ]);
                        pause(0.05)
                    end
                 end
                 function Save_GIF_scalar(nodes,name,scalar,Range)
                     
                     JJ=jet(256);
                                       
                      img_stack=uint8([]);
                        cm=[];
                        xl=[];
                        yl=[];
                        zl=[];
                        h=figure();
                        for cpt_t=1:1:size(nodes(1).points,3)
                            
                              Color=(scalar(:,cpt_t)-Range(1))/(Range(2)-Range(1));
                              Color(Color<0)=0;
                              Color(Color>1)=1;
                              Color=round(Color*255)+1;
                            
                            
                            clf(h)
                            hold on
                            for cpt_slc=1:1:size(nodes,2) 
                                 scatter3(nodes(cpt_slc).points(:,1,cpt_t),nodes(cpt_slc).points(:,2,cpt_t),nodes(cpt_slc).points(:,3,cpt_t),50,Color,'filled')         
                            end
                          
                            view(-140, 26)
                            set(gca,'Clipping','off');

                            set(gca,'xtick',[]);
                            set(gca,'xTickLabel',[]);
                            set(gca,'ytick',[]);
                            set(gca,'yTickLabel',[]);
                            set(gca,'ztick',[]);
                            set(gca,'zTickLabel',[]);

                            set(gca,'visible','off')

                            box off
                            set(gcf,'color','w');
                            if cpt_t==1
                                xl = xlim;
                                yl = ylim;
                                zl = zlim;
                            else
                                xlim([xl(1) xl(2)]) ;
                                ylim([yl(1) yl(2)]);
                                zlim([zl(1) zl(2)]);
                            end

                            pause(0.05)
                            img=getframe(gcf);
                            [tmpp,~]=frame2im(img);
                            img_stack(:,:,:,end+1)=tmpp;
                        end


                         %[tmp_frame,cm] = rgb2ind(squeeze(img_stack(:,:,:,1)),256);
                         for cpt=2:1:size(img_stack,4)
                             [tmp_frame,cm] = rgb2ind(squeeze(img_stack(:,:,:,cpt)),256);
                             if cpt==2
                                 imwrite( uint8(tmp_frame),cm,[name '.gif'],'gif', 'Loopcount',inf,'DelayTime',0.1);   %%%% First image, delay time = 0.1s         
                             else
                                imwrite( uint8(tmp_frame),cm,[name '.gif'],'gif','WriteMode','append','DelayTime',0.1); %%%% Following images
                             end 
                         end
                     
                 end
                 function Save_GIF(nodes,name,Color)
                     
                %% PLOT nodes in 3D
                img_stack=uint8([]);
                cm=[];
                xl=[];
                yl=[];
                zl=[];
                h=figure();
                for cpt_t=1:1:size(nodes(1).points,3)
                    clf(h)
                    hold on
                    for cpt_slc=1:1:size(nodes,2) 
                        scatter3(nodes(cpt_slc).points(:,1,cpt_t),nodes(cpt_slc).points(:,2,cpt_t),nodes(cpt_slc).points(:,3,cpt_t),32,'k','filled')
                         scatter3(nodes(cpt_slc).points(:,1,cpt_t),nodes(cpt_slc).points(:,2,cpt_t),nodes(cpt_slc).points(:,3,cpt_t),30,squeeze(Color(:,cpt_slc))','filled')
                        % plot3(nodes_DENSE(cpt_slc).ROI.phase(cpt_t).endo(:,1),nodes_DENSE(cpt_slc).ROI.phase(cpt_t).endo(:,2),nodes_DENSE(cpt_slc).ROI.phase(cpt_t).endo(:,3),'r-')  
                        % plot3(nodes_DENSE(cpt_slc).ROI.phase(cpt_t).epi(:,1) ,nodes_DENSE(cpt_slc).ROI.phase(cpt_t).epi(:,2) ,nodes_DENSE(cpt_slc).ROI.phase(cpt_t).epi(:,3),'b-')
                     
                    end
                   % scatter3(nodes_DTI2.points(:,1,cpt_t),nodes_DTI2.points(:,2,cpt_t),nodes_DTI2.points(:,3,cpt_t),20,'filled','k')
                   % plot3(nodes_DTI2.ROI.phase(cpt_t).endo(:,1),nodes_DTI2.ROI.phase(cpt_t).endo(:,2),nodes_DTI2.ROI.phase(cpt_t).endo(:,3),'r-')  
                   % plot3(nodes_DTI2.ROI.phase(cpt_t).epi(:,1) ,nodes_DTI2.ROI.phase(cpt_t).epi(:,2) ,nodes_DTI2.ROI.phase(cpt_t).epi(:,3),'b-')
                   %  axis([-60 60 -60 60 -60 60]);
                  %  view(-126, -7)
                   % view(-134, -16)
                  view( -130,-51)
                 %    view(-140, -25)
                    set(gca,'Clipping','off');

                    set(gca,'xtick',[]);
                    set(gca,'xTickLabel',[]);
                    set(gca,'ytick',[]);
                    set(gca,'yTickLabel',[]);
                    set(gca,'ztick',[]);
                    set(gca,'zTickLabel',[]);
                    
                    set(gca,'visible','off')

                    box off
                    set(gcf,'color','w');
                    if cpt_t==1
                        xl = xlim;
                        yl = ylim;
                        zl = zlim;
                    else
                        xlim([xl(1) xl(2)]) ;
                        ylim([yl(1) yl(2)]);
                        zlim([zl(1) zl(2)]);
                    end

                    pause(0.05)
                    img=getframe(gcf);
                    [tmpp,~]=frame2im(img);
                    img_stack(:,:,:,end+1)=tmpp;
                end


                 %[tmp_frame,cm] = rgb2ind(squeeze(img_stack(:,:,:,1)),256);
                 for cpt=2:1:size(img_stack,4)
                     [tmp_frame,cm] = rgb2ind(squeeze(img_stack(:,:,:,cpt)),256);
                     if cpt==2
                         imwrite( uint8(tmp_frame),cm,[name '.gif'],'gif', 'Loopcount',inf,'DelayTime',0.1);   %%%% First image, delay time = 0.1s         
                     else
                        imwrite( uint8(tmp_frame),cm,[name '.gif'],'gif','WriteMode','append','DelayTime',0.1); %%%% Following images
                     end 
                 end
              end
            
                 
                 
             function Save_GIF_ff(nodes,name,Color)
                     
                %% PLOT nodes in 3D
                img_stack=uint8([]);
                cm=[];
                xl=[];
                yl=[];
                zl=[];
                h=figure();
                for cpt_t=1:1:size(nodes(1).points,3)
                    clf(h)
                    hold on
                    for cpt_slc=1:1:size(nodes,2) 
                        CC=abs(inv(nodes.Rotation)*squeeze(Color(:,:,cpt_t)'))';
                        scatter3(nodes(cpt_slc).points(:,1,cpt_t),nodes(cpt_slc).points(:,2,cpt_t),nodes(cpt_slc).points(:,3,cpt_t),32,'k','filled')
                         scatter3(nodes(cpt_slc).points(:,1,cpt_t),nodes(cpt_slc).points(:,2,cpt_t),nodes(cpt_slc).points(:,3,cpt_t),30,CC,'filled')
                        % plot3(nodes_DENSE(cpt_slc).ROI.phase(cpt_t).endo(:,1),nodes_DENSE(cpt_slc).ROI.phase(cpt_t).endo(:,2),nodes_DENSE(cpt_slc).ROI.phase(cpt_t).endo(:,3),'r-')  
                        % plot3(nodes_DENSE(cpt_slc).ROI.phase(cpt_t).epi(:,1) ,nodes_DENSE(cpt_slc).ROI.phase(cpt_t).epi(:,2) ,nodes_DENSE(cpt_slc).ROI.phase(cpt_t).epi(:,3),'b-')
                     
                    end
                   % scatter3(nodes_DTI2.points(:,1,cpt_t),nodes_DTI2.points(:,2,cpt_t),nodes_DTI2.points(:,3,cpt_t),20,'filled','k')
                   % plot3(nodes_DTI2.ROI.phase(cpt_t).endo(:,1),nodes_DTI2.ROI.phase(cpt_t).endo(:,2),nodes_DTI2.ROI.phase(cpt_t).endo(:,3),'r-')  
                   % plot3(nodes_DTI2.ROI.phase(cpt_t).epi(:,1) ,nodes_DTI2.ROI.phase(cpt_t).epi(:,2) ,nodes_DTI2.ROI.phase(cpt_t).epi(:,3),'b-')
                   %  axis([-60 60 -60 60 -60 60]);
                  %  view(-126, -7)
                    view(-140, -25)
                   
                    set(gca,'Clipping','off');

                    set(gca,'xtick',[]);
                    set(gca,'xTickLabel',[]);
                    set(gca,'ytick',[]);
                    set(gca,'yTickLabel',[]);
                    set(gca,'ztick',[]);
                    set(gca,'zTickLabel',[]);
                    
                    set(gca,'visible','off')

                    box off
                    set(gcf,'color','w');
                    if cpt_t==1
                        xl = xlim;
                        yl = ylim;
                        zl = zlim;
                    else
                        xlim([xl(1) xl(2)]) ;
                        ylim([yl(1) yl(2)]);
                        zlim([zl(1) zl(2)]);
                    end

                    pause(0.05)
                    img=getframe(gcf);
                    [tmpp,~]=frame2im(img);
                    img_stack(:,:,:,end+1)=tmpp;
                end


                 %[tmp_frame,cm] = rgb2ind(squeeze(img_stack(:,:,:,1)),256);
                 for cpt=2:1:size(img_stack,4)
                     [tmp_frame,cm] = rgb2ind(squeeze(img_stack(:,:,:,cpt)),256);
                     if cpt==2
                         imwrite( uint8(tmp_frame),cm,[name '.gif'],'gif', 'Loopcount',inf,'DelayTime',0.1);   %%%% First image, delay time = 0.1s         
                     else
                        imwrite( uint8(tmp_frame),cm,[name '.gif'],'gif','WriteMode','append','DelayTime',0.1); %%%% Following images
                     end 
                 end
             end
             function [img]=fiber2im(nodes,sz)
                        
                 img=[];
                 tmp_point=[];
                   for cpt=1:1:size(nodes.points,3)
                       tmp_point(:,:,cpt)=(inv(nodes.Rotation)*nodes.points(:,:,cpt)')';   
                       fiberRot(:,:,cpt) =(inv(nodes.Rotation)*nodes.ff(:,:,cpt)')';
                       
                   end
                   tmp_point=tmp_point(:,1:2,:); % Keep only the 2d data;   
                   
                   Bary_time=squeeze(mean(mean(tmp_point),3));
                   Max_time= max(max(abs(tmp_point)),[],2);
             
                   [X Y]=meshgrid(linspace(-Max_time(:,1),Max_time(:,1),sz),linspace(-Max_time(:,2),Max_time(:,2),sz));
                   
                   X=X+Bary_time(1);
                   Y=Y+Bary_time(2);
                   
                 for phase=1:1:size(nodes.points,3) 
                   point_endo= (inv(nodes.Rotation)*nodes.ROI.phase(phase).endo(:,1:3)')';
                   point_epi = (inv(nodes.Rotation)*nodes.ROI.phase(phase).epi(:,1:3)')';
                        
                   point_endo=point_endo(:,1:2);
                   point_epi=point_epi(:,1:2);
                   
                   CoorGrid(:,1)=X(:);
                   CoorGrid(:,2)=Y(:);
                   
                   Col_endo=Collision_ToolBox.ROI_Points(point_endo,CoorGrid);
                   Col_epi=Collision_ToolBox.ROI_Points(point_epi,CoorGrid);
                        
                   Col=(~Col_endo)&(Col_epi);
                     
                   Fx = scatteredInterpolant(tmp_point(:,1,phase),tmp_point(:,2,phase),abs(fiberRot(:,1,phase)),'natural','linear');
                   Fy = scatteredInterpolant(tmp_point(:,1,phase),tmp_point(:,2,phase),abs(fiberRot(:,2,phase)),'natural','linear');          
                   Fz = scatteredInterpolant(tmp_point(:,1,phase),tmp_point(:,2,phase),abs(fiberRot(:,3,phase)),'natural','linear');
                   
                   GridX=zeros(size(CoorGrid,1),3);

                   for cpt_p=1:1:size(CoorGrid,1)
                        if Col(cpt_p)
                            GridX(cpt_p,1)=Fx(CoorGrid(cpt_p,1), CoorGrid(cpt_p,2));
                            GridX(cpt_p,2)=Fy(CoorGrid(cpt_p,1), CoorGrid(cpt_p,2));
                            GridX(cpt_p,3)=Fz(CoorGrid(cpt_p,1), CoorGrid(cpt_p,2));
                        end
                   end
                   img(:,:,:,phase)=reshape(GridX,sz,sz,3);
  %                  figure,scatter(CoorGrid(:,1),CoorGrid(:,2),40,Col)
 %                   hold on
 %                   scatter(tmp_point(:,1,phase),tmp_point(:,2,phase),40,'r')
                end
               end
             function [img]=nodes2im(nodes,value,sz)
                        
                 img=[];
                 tmp_point=[];
                   for cpt=1:1:size(nodes.points,3)
                       tmp_point(:,:,cpt)=(inv(nodes.Rotation)*nodes.points(:,:,cpt)')';                      
                   end
                   tmp_point=tmp_point(:,1:2,:); % Keep only the 2d data;   
                   
                   Bary_time=squeeze(mean(mean(tmp_point),3));
                   Max_time= max(max(abs(tmp_point)),[],2);
             
                   [X Y]=meshgrid(linspace(-Max_time(:,1),Max_time(:,1),sz),linspace(-Max_time(:,2),Max_time(:,2),sz));
                   
                   X=X+Bary_time(1);
                   Y=Y+Bary_time(2);
                   
                 for phase=1:1:size(nodes.points,3) 
                   point_endo= (inv(nodes.Rotation)*nodes.ROI.phase(phase).endo(:,1:3)')';
                   point_epi = (inv(nodes.Rotation)*nodes.ROI.phase(phase).epi(:,1:3)')';
                        
                   point_endo=point_endo(:,1:2);
                   point_epi=point_epi(:,1:2);
                   
                   CoorGrid(:,1)=X(:);
                   CoorGrid(:,2)=Y(:);
                   
                   Col_endo=Collision_ToolBox.ROI_Points(point_endo,CoorGrid);
                   Col_epi=Collision_ToolBox.ROI_Points(point_epi,CoorGrid);
                        
                   Col=(~Col_endo)&(Col_epi);
                   
                  
                   
                   F = scatteredInterpolant(tmp_point(:,1,phase),tmp_point(:,2,phase),value(:,phase));
                   GridEff=zeros(size(CoorGrid,1),1);

                   for cpt_p=1:1:size(CoorGrid,1)
                        if Col(cpt_p)
                            GridEff(cpt_p)=F(CoorGrid(cpt_p,1), CoorGrid(cpt_p,2));
                        end
                   end
                   img(:,:,phase)=reshape(GridEff,sz,sz);
%                    figure,scatter(CoorGrid(:,1),CoorGrid(:,2),40,Col)
%                    hold on
%                    scatter(tmp_point(:,1,phase),tmp_point(:,2,phase),40,'r')
                end
               end
            
        
    end
end
