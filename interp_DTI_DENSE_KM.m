function interp_DTI_DENSE_KM(dcm_dir)
warning off

Folder=[];
Folder.dcm_dir = dcm_dir;
cd(Folder.dcm_dir);


Folder.Dense_files=dir ([Folder.dcm_dir '\DENSE\*.mat']);

Folder.DENSE = [Folder.dcm_dir '\DENSE\'];
Folder.DTI   = [Folder.dcm_dir '\DTI\'];
Folder.Result= [Folder.dcm_dir '\Result\'];

mkdir(Folder.Result);

DENSE_param=[];

DirZ_vect=[0;0;1];


%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DENSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[DENSE_param, nodes_DENSE] = Compute_DENSE_nodes(Folder);
%[DENSE_param, nodes_DENSE] = Compute_DENSE_nodes_DNS('Seg2_test_volume_2.dns');
%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DTI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[DTI_param, nodes_DTI] = Compute_DTI_nodes(Folder);
DTI_param.CardiacPhase = round(DTI_param.TD/DENSE_param(1).temporal_res ) ;
%%

%% Register DTI on DENSE (red before, green after) 
[nodes_DTI]=Eff_Toolbox.Register_nodes(nodes_DENSE(2),nodes_DTI,DTI_param.CardiacPhase,[1,0,0;0,1,0]'); % Register the X and Y direction of the DTI slice and the middle DENSE slice. The selection of the DENSE slice is still hardcoded. TODO. 
[nodes_DTI]=Eff_Toolbox.Register_optim_nodes(nodes_DENSE(2),nodes_DTI,DTI_param.CardiacPhase);


%% Query 10 point in the top slice and 10 point in the bottom slice
[nodes_DTI2, Fscat_save2] = ComputeF_Exp_KM4(nodes_DENSE, nodes_DTI, DTI_param.CardiacPhase);


%% Compute strain correct from phase 1 GOOD
[nodes_DTI3, Fscat_save3] = ComputeF_Exp_KM4(nodes_DENSE, nodes_DTI2, DENSE_param(1).RefPhase);

%%
[nodes_DTI3.HA_nodes2 nodes_DTI3.E2A_nodes2 nodes_DTI3.TRA_nodes2]= HA_E2A_nodes_KM2(nodes_DTI3);
[nodes_DTI3.HA_nodes_f]=HA_Filter_Points_KM(nodes_DTI3.HA_nodes2,nodes_DTI3.wd);
nodes_DTI3.HAR=Eff_Toolbox.HA_WD(nodes_DTI3.wd,nodes_DTI3.HA_nodes_f);


%%
save([Folder.Result 'Eff_n01.mat'],'DENSE_param','DTI_param','Fscat_save2','Fscat_save3','nodes_DENSE','nodes_DTI','nodes_DTI2','nodes_DTI3');

end