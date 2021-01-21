function [HA E2A TRA]= HA_E2A_nodes_KM2(nodes_DTI)

% Calculates the angle (degrees) between the primary eigenvector and the SA
% plane; epicardium and endocardium should be a series of points representing the
% boundaries of the myocardium.
%
% SYNTAX:  HA_KM(  EigVect1, Mask, P_Epi, P_Endo)
%
% INPUTS:   EigVect1 - First EigVector image matrix
%                 [y x slices coordinates]
%           
%           Mask -  Mask matrix
%                 [y x slices]
%
%           P_Endo - List of Coordinates of the Endocardium ROI
%                 
%           P_Epi - List of Coordinates of the Endocardium ROI
%
% OUTPUTS:  HA - HA image matrix (units [- pi pi])
%                 [y x slices]
%
% ???? 08.14.2017
% ?????
% Ennis Lab @ UCLA; http://mrrl.ucla.edu





%% using ellipses

%%

disp('Generate HA') 
h = waitbar(0,'Generate HA...');
%HA = zeros(size(EigVect1,1),size(C_Epi,3));
%TRA = zeros(size(EigVect1,1),size(C_Epi,3));

for cpt_t=1:1:size(nodes_DTI.points,3)
    
    
    % Approximate everything in 2D
    P_Epi=(inv(nodes_DTI.Rotation)*squeeze(squeeze(nodes_DTI.ROI.phase(cpt_t).epi))')';
    P_Endo=(inv(nodes_DTI.Rotation)*squeeze(squeeze(nodes_DTI.ROI.phase(cpt_t).endo))')';
    Points_Rot=(inv(nodes_DTI.Rotation)*squeeze(nodes_DTI.points(:,:,cpt_t))')';  
    
    Vec(1,:) = P_Epi(1,:) - P_Epi(end,:);
    for y = 2:size(P_Epi,1)
        Vec(y,:) = P_Epi(y,:) - P_Epi(y-1,:);
    end

    Vec2(1,:) = P_Endo(1,:) - P_Endo(end,:);
    for y = 2:size(P_Endo,1)
        Vec2(y,:) = P_Endo(y,:) - P_Endo(y-1,:);
    end
    
    positions = cat(1,P_Epi(:,:),P_Endo(:,:));
    vectors   = cat(1,Vec,Vec2);  
    
    Vy = griddata(positions(:,1),positions(:,2),vectors(:,2),Points_Rot(:,1),Points_Rot(:,2));
    Vx = griddata(positions(:,1),positions(:,2),vectors(:,1),Points_Rot(:,1),Points_Rot(:,2));
    
    for cpt_p = 1:size(nodes_DTI.points,1)                         
                                                
                E1 = inv(nodes_DTI.Rotation)*squeeze(nodes_DTI.ff(cpt_p,:,cpt_t))';
                E2 = inv(nodes_DTI.Rotation)*squeeze(nodes_DTI.ff2(cpt_p,:,cpt_t))';
                
                 % Circunferiential Vector definition
                Circ = [Vx(cpt_p) Vy(cpt_p) 0];
                Circ= Circ/norm(Circ);
                % Longitudinal Vector definition
                Long= [0 0 1];
                
                % Radial Vector definition
                Rad = cross(Circ/norm(Circ),Long/norm(Long));
                
                % Projection of the Fiber Vector onto the Circunferential
                % direction 
                E1proj_circ = Vector_ToolBox.Projection_vect(E1, Circ);
                E1proj_circ(3)=E1(3);
                %E1proj=dot(E1,Circ)*Circ/(norm(Circ)^2);
%                 Fiber_proj=[E1proj(1) E1proj(2) E1(3)];
%            
%                 Fiber_proj=Fiber_proj./norm(Fiber_proj);   
                E1proj_long = Vector_ToolBox.Projection_vect(E1, Long); % Z component of the vector
                
                %HA(cpt_p,cpt_t) = atan2(norm(E1proj_long),norm(E1proj_circ))*180/(pi);%
                 HA(cpt_p,cpt_t) =asin(E1proj_circ(3)/norm(E1proj_circ))*180/(pi);
                                
                MidFiber=cross(E1proj_circ/norm(E1proj_circ),Rad/norm(Rad));
                %MidFiber=cross(E1/norm(E1),Rad/norm(Rad));
                
                % Projection of the Sheet Vector onto the Radial
                % direction 
                E2proj = Vector_ToolBox.Projection_vect(E2, Rad);
                
                % Projection of the Sheet Vector onto the MidFiber
                % direction 
                E2proj_mid = Vector_ToolBox.Projection_vect(E2, MidFiber);
                
                % E2A represent the angle between these two projections              
                E2A(cpt_p,cpt_t) = atan2(norm(E2proj),norm(E2proj_mid))*180/(pi);
                
                
                
                E1proj_circ = Vector_ToolBox.Projection_vect(E1, Circ);
                TRA(cpt_p,cpt_t)= acos(dot(E1,Circ)/(norm(Circ)*norm(E1)))*180/(pi);
                
                
            
                if dot(Circ,E1) < 0 %&& HA(j,k) > 0  %dot(tang,tProj) > 0
                    HA(cpt_p,cpt_t) = -HA(cpt_p,cpt_t);
                end
                
                if dot(E1,Circ) < 0 %&& HA(j,k) > 0
                    TRA(cpt_p,cpt_t) = TRA(cpt_p,cpt_t)-180;
                end
        end
        waitbar(cpt_t/size(nodes_DTI.points,3),h);
    end

close(h)


%HA(MASK==0) = nan;

end

