classdef Vector_ToolBox

    methods(Static)
        
                 %% Make the projection of a list of Vector [n x coordinates] onto a vector [x y z]
                 function Vector_proj= Projection_vect_n(Vector_list, Vector_dir)
                      Vector_proj=zeros(size(Vector_list));
                      for cpt_vec=1:1:size(Vector_list,1)                         
                            Vector_proj(cpt_vec,:)=Vector_ToolBox.Projection_vect(Vector_list(cpt_vec,:),Vector_dir);
                      end    
                 end
        
                 %% Make the projection of a list of Vector [n x coordinates] onto a vector [x y z]
                 function Vector_proj= Projection_vect_n_t(Vector_list, Vector_dir)
                      Vector_proj=zeros(size(Vector_list));
                      for cpt_t=1:1:size(Vector_list,3)                         
                            Vector_proj(:,:,cpt_t)=Vector_ToolBox.Projection_vect_n(Vector_list(:,:,cpt_t),Vector_dir);
                      end    
                 end
                 
                %% Make the projection of a  Vector [n x coordinates] onto a vector [x y z]
                 function Vector_proj = Projection_vect(Vector, Vector_dir)                 
%                       A = [-10,10,0];
%                       B = [0,0,1];
%                       C = (dot(A,B)/norm(B)^2)*B;
                      Vector_proj = (dot(Vector,Vector_dir)/norm(Vector_dir)^2)*Vector_dir;
                 end
                 
                 %% Make the projection of a list of Vector [n x coordinates] onto a vector [x y z]
                 function Vnorm= Norm_vect_n(Vector_list)
                      Vnorm=zeros(size(Vector_list,1),1);
                      for cpt_vec=1:1:size(Vector_list,1)      
                        Vnorm(cpt_vec)=  Vector_ToolBox.Norm_vect(Vector_list(cpt_vec,:));
                      end
                 end
                 %% Make the projection of a list of Vector [n x coordinates] onto a vector [x y z]
                 function Vnorm= Norm_vect(Vector)
                      Vnorm=  sqrt(Vector(1)*Vector(1)+Vector(2)*Vector(2)+Vector(3)*Vector(3));
                 end
    end
    
end