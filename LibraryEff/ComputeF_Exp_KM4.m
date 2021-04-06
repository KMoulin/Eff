% function [Eff, Err, Ecc, Ell, Jac, V_t, V_t_nodes] = ComputeF_Exp(phi_all, phi_up, phi_down, query_up, query_down, dN_all, QPmid_im, cc, rr, ll, fail_index, f, phase)
function [nodes_dti2, Fscat_save] = ComputeF_Exp_KM4(nodes_dense, nodes_dti, ref_phase)



phi_up=nodes_dense(1).points;
phi_mid=nodes_dense(2).points;
phi_down=nodes_dense(3).points;
phi_all = [nodes_dense(1).points;nodes_dense(3).points];



nodes_dti2=nodes_dti;
nodes_dti2.dT=nodes_dense(1).dT;
% Select the DTI points we are going to move.  
QDTI_Points=[];
QROI_Points_endo=[];
QROI_Points_epi=[];
ref_phase_DTI=ref_phase;

% If there is only one cardiac phase for DTI
if size(nodes_dti.points,3)<ref_phase
    ref_phase_DTI=1;
end

QDTI_Points=      nodes_dti.points(:,1:3,ref_phase_DTI);
QROI_Points_endo= nodes_dti.ROI.phase(ref_phase_DTI).endo(:,:);
QROI_Points_epi=  nodes_dti.ROI.phase(ref_phase_DTI).epi(:,:);



for cpt_t = 1:size(phi_all,3)
    %DENSE
    
     % Intial config at the given phase from DTI (reference config)
    X_radius = [phi_up(:,:,ref_phase);phi_mid(:,:,ref_phase);phi_down(:,:,ref_phase)];
        
    % Position throught time
    x_radius = [phi_up(:,:,cpt_t);phi_mid(:,:,cpt_t);phi_down(:,:,cpt_t)];   
    
    % Interpolant of position in X Y Z as a regard of initial config. 
    Phi_x = scatteredInterpolant(X_radius(:,1), X_radius(:,2), X_radius(:,3), x_radius(:,1),'natural','linear');
    Phi_y = scatteredInterpolant(X_radius(:,1), X_radius(:,2), X_radius(:,3), x_radius(:,2),'natural','linear');
    Phi_z = scatteredInterpolant(X_radius(:,1), X_radius(:,2), X_radius(:,3), x_radius(:,3),'natural','linear');
    nodes_dti2.Phi_x(cpt_t).struct=Phi_x;
    nodes_dti2.Phi_y(cpt_t).struct=Phi_y;
    nodes_dti2.Phi_z(cpt_t).struct=Phi_z;
    % Apply deformation for each dTI point
    for cpt_p = 1:size(QDTI_Points,1)
        
        % Current DTI point
        DTI_point = squeeze(QDTI_Points(cpt_p,:));
        
        % Local definition of the orientation vectors
        ff = nodes_dti.ff(cpt_p,:,ref_phase_DTI) /norm(nodes_dti.ff(cpt_p,:,ref_phase_DTI));
        ff2= nodes_dti.ff2(cpt_p,:,ref_phase_DTI)/norm(nodes_dti.ff2(cpt_p,:,ref_phase_DTI));
        ff3= nodes_dti.ff3(cpt_p,:,ref_phase_DTI)/norm(nodes_dti.ff3(cpt_p,:,ref_phase_DTI));
        cc = nodes_dti.cc(cpt_p,:,ref_phase_DTI) /norm(nodes_dti.cc(cpt_p,:,ref_phase_DTI));
        rr = nodes_dti.rr(cpt_p,:,ref_phase_DTI) /norm(nodes_dti.rr(cpt_p,:,ref_phase_DTI));
        ll = nodes_dti.ll(cpt_p,:,ref_phase_DTI) /norm(nodes_dti.ll(cpt_p,:,ref_phase_DTI));       
        wd = nodes_dti.wd(cpt_p,ref_phase_DTI) ;
        
        Fscat = Compute_F_local (DTI_point,Phi_x,Phi_y,Phi_z);
        Fscat_save(:,:,cpt_p,cpt_t)=Fscat;
        
        %% Apply the deformation to the vectors

        tmp_ff = Fscat*ff';
        nodes_dti2.ff(cpt_p,:,cpt_t) = tmp_ff/norm(tmp_ff);

        tmp_ff2 = Fscat*ff2';
        nodes_dti2.ff2(cpt_p,:,cpt_t) = tmp_ff2/norm(tmp_ff2);

        
        tmp_ff3 = Fscat*ff3';
        nodes_dti2.ff2(cpt_p,:,cpt_t) = tmp_ff3/norm(tmp_ff3);
        
        tmp_rr = Fscat*rr';
        nodes_dti2.rr(cpt_p,:,cpt_t) = tmp_rr/norm(tmp_rr);

        tmp_cc = Fscat*cc';           
        nodes_dti2.cc(cpt_p,:,cpt_t) = tmp_cc/norm(tmp_cc);

        tmp_ll = Fscat*ll';
        nodes_dti2.ll(cpt_p,:,cpt_t) = tmp_ll/norm(tmp_ll);

        
        nodes_dti2.Jac(cpt_p,cpt_t) = det(Fscat);

        
        nodes_dti2.points(cpt_p,:,cpt_t) = [Phi_x(DTI_point); Phi_y(DTI_point); Phi_z(DTI_point)];
        
        nodes_dti2.wd(cpt_p,cpt_t)=wd;
        
        % For strain C= F' F 
        C = transpose(Fscat) * Fscat;
        
        

        E= 0.5 *(C - eye(3));
        
        % cc
        nodes_dti2.Ecc(cpt_p,cpt_t) = cc * E * transpose(cc);

        % rr
        nodes_dti2.Err(cpt_p,cpt_t) = rr * E * transpose(rr);

        % ll
        nodes_dti2.Ell(cpt_p,cpt_t) = ll * E * transpose(ll);

        % ff
        nodes_dti2.Eff(cpt_p,cpt_t) = ff * E * transpose(ff);
        
        % ff2
        nodes_dti2.Eff2(cpt_p,cpt_t) = ff2 * E * transpose(ff2);
        
         % ff3
        nodes_dti2.Eff3(cpt_p,cpt_t) = ff3 * E * transpose(ff3);
        
        
         % ff3
        nodes_dti2.Ecl(cpt_p,cpt_t) = cc * E * transpose(ll);
    end
    
    % Apply deformation for ROI endo 
    for cpt_p = 1:size(QROI_Points_endo,1)
        nodes_dti2.ROI.phase(cpt_t).endo(cpt_p,:) = [Phi_x(QROI_Points_endo(cpt_p,1:3)); Phi_y(QROI_Points_endo(cpt_p,1:3)); Phi_z(QROI_Points_endo(cpt_p,1:3))];
    end
    
    % Apply deformation for ROI epi 
    for cpt_p = 1:size(QROI_Points_epi,1)
        nodes_dti2.ROI.phase(cpt_t).epi(cpt_p,:) = [Phi_x(QROI_Points_epi(cpt_p,1:3)); Phi_y(QROI_Points_epi(cpt_p,1:3)); Phi_z(QROI_Points_epi(cpt_p,1:3))];
    end
       
    %% Calculate new Wd    
    
    
    nodes_dti2.wd(:,cpt_t) = Compute_WD_local(nodes_dti2.points(:,:,cpt_t),nodes_dti2.ROI.phase(cpt_t).endo,nodes_dti2.ROI.phase(cpt_t).epi,nodes_dti2.Rotation);

   % disp([ num2str(cpt_t/size(phi_all,3)) ' Computing Strains... ' ]);
end
end


function F = Compute_F_local (point,Phi_x,Phi_y,Phi_z)
       %Defining Fscat
        Delta = 0.1; % was in pixel should be in mm (0.1 originally)
       % 11 21 31
        Xq_plus  = point; 
        Xq_plus(1)  = Xq_plus(1)  + Delta;
        Xq_minus = point; 
        Xq_minus(1) = Xq_minus(1) - Delta;

        % J is the colum coordinate in the reference system in which we are
        % taking the derivative
        % i is the line the deformation mapping conponement we are
        % differentiate. 
        
        F(1,1) = (Phi_x(Xq_plus) - Phi_x(Xq_minus))/(2*Delta);
        F(2,1) = (Phi_y(Xq_plus) - Phi_y(Xq_minus))/(2*Delta);
        F(3,1) = (Phi_z(Xq_plus) - Phi_z(Xq_minus))/(2*Delta);

        % 12 22 32
        Xq_plus  = point; 
        Xq_plus(2)  = Xq_plus(2)  + Delta;
        Xq_minus = point; 
        Xq_minus(2) = Xq_minus(2) - Delta;

        F(1,2) = (Phi_x(Xq_plus) - Phi_x(Xq_minus))/(2*Delta);
        F(2,2) = (Phi_y(Xq_plus) - Phi_y(Xq_minus))/(2*Delta);
        F(3,2) = (Phi_z(Xq_plus) - Phi_z(Xq_minus))/(2*Delta);

        % 13 23 33
        Xq_plus  = point; 
        Xq_plus(3)  = Xq_plus(3)  + Delta;
        Xq_minus = point; 
        Xq_minus(3) = Xq_minus(3) - Delta;

        % F for Phi
        F(1,3) = (Phi_x(Xq_plus) - Phi_x(Xq_minus))/(2*Delta);
        F(2,3) = (Phi_y(Xq_plus) - Phi_y(Xq_minus))/(2*Delta);
        F(3,3) = (Phi_z(Xq_plus) - Phi_z(Xq_minus))/(2*Delta);


end

function WD= Compute_WD_local (Points,P_Endo,P_Epi,Rot)

    Endo_Line = zeros(size(P_Endo,1),1);
    Epi_Line  = ones(size(P_Epi,1),1);
    PosRoi    = cat(1,P_Endo,P_Epi);
    LineRoi   = cat(1,Endo_Line,Epi_Line);
    
    
    %% Approximate in 2D 
    PosRoi_Rot=(inv(Rot)*(PosRoi)')';
    Points_Rot=(inv(Rot)*(Points)')';
    % griddata(X,Y,V, xq, yq)
    WD=griddata(PosRoi_Rot(:,1),PosRoi_Rot(:,2),LineRoi,Points_Rot(:,1),Points_Rot(:,2));
    
%     Idx_exist=find(~isnan(WD));
%     Idx_nan=find(isnan(WD));
%     Scatter_WD = scatteredInterpolant(Points(Idx_exist,1), Points(Idx_exist,2), Points(Idx_exist,3), WD(Idx_exist),'linear','linear');
%     WD(Idx_nan)=Scatter_WD(Points(Idx_nan,1), Points(Idx_nan,2), Points(Idx_nan,3));
    
end
