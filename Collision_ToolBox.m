classdef Collision_ToolBox

    methods(Static)
        
                %% General function to manage Collision between Objects and Points
                 function In_out_List = Collision_Detection(Objects, Points)
                      In_out_List=false(size(Points,1),1);
                      In_out=false(size(Points,1),1);

                      for cpt_cell=1:1:size(Objects,1)
                          if size(Objects,2)>7 % Cell is a polygon 
                            In_out= Collision_ToolBox.Poly_Points_Bounding_Box(squeeze(Objects(cpt_cell,:)),Points); % 
                          else               % Cell is a circle
                            In_out= Collision_ToolBox.Circle_Points(squeeze(Objects(cpt_cell,:)),Points);  
                          end 
                          In_out_List=In_out_List+In_out; % List of the molecule associate to the corresponding cells  
                      end
                      
                 end
                %% General function to manage Collision between Objects and Points
                 function In_out_List = Collision_Detection_Mask(Mask, Points,Resolution)
                      In_out_List=false(size(Points,1),1);
                      
                      % From absolute coordinate to pixel coordinate  from
                      % 1 to 500
                      for cpt_dim=1:1:3
                          Points(:,cpt_dim)=round(Points(:,cpt_dim)./Resolution(cpt_dim))+1;
                          Idx=find(Points(:,cpt_dim)>size(Mask,cpt_dim));
                          Points(Idx,cpt_dim)=size(Mask,cpt_dim);
                          
                          Idx=find(Points(:,cpt_dim)<1);
                          Points(Idx,cpt_dim)=1;
                      end 
                      Idx=sub2ind(size(Mask),Points(:,1),Points(:,2),Points(:,3));
                      In_out_List=Mask(Idx);
                 end
                %% Function which manage the permeability 
                function [In_out] = Permeability(In_out_before,In_out_after,Perma)
                           List_diff=In_out_before-In_out_after; % 0 nothing ; -1 From cells to extra ; 1 From extra to cells
                           List_perma=find(List_diff~=0);
                           Roll_perma=(rand(length(List_perma),1));
                           In_out=List_perma(find(Roll_perma(Roll_perma>Perma)));
                end
                
                %% Collision between Circle and a list of points
                function   In_out=Circle_Points(Circle,Points)
                      tmp_dist =  sqrt( (Circle(3)-Points(:,1)).^2 + (Circle(4)-Points(:,2)).^2);
                      in_xy = ( tmp_dist<=Circle(1) );
                      in_z =  ( Points(:,3)>= ( Circle(5)  ) ) & ( Points(:,3)< ( Circle(5) + Circle(2) ) );       
                      In_out= in_xy & in_z ;
                end
                
                  %% Collision between Circle and a list of points
                function   In_out=Z_Plane(Object1,Object2)
                      In_out =  ( Object2(5)>= ( Object1(5)  ) ) & ( Object2(5)< ( Object1(5) + Object1(2) ) );       
                end
                
                %% Collision between Polynome and a list of points
                function   In_out=Poly_Points(Poly,Points)     
                     % Poly [Radius Length Pos_X Pos_Y Pos_Z %Surface %Volume Nb_Poly p1 p2 p3 p4 ..]
                       in_xy =false(1,size(Points,1));
                      for cpt_cell=1:1:Poly(8)
                           p1(:,1)=Poly(3)+Poly(8+cpt_cell)*(cos((2*pi*(cpt_cell-1))/(Poly(8))));
                           p1(:,2)=Poly(4)+Poly(8+cpt_cell)*(sin((2*pi*(cpt_cell-1))/(Poly(8))));
                           if cpt_cell==Poly(8)
                                p2(1)=Poly(3)+Poly(8+1)*cos(0);
                                p2(2)=Poly(4)+Poly(8+1)*sin(0);
                           else
                                p2(1)=Poly(3)+Poly(8+cpt_cell+1)*(cos((2*pi*(cpt_cell))/(Poly(8))));
                                p2(2)=Poly(4)+Poly(8+cpt_cell+1)*(sin((2*pi*(cpt_cell))/(Poly(8))));
                           end
                           tmp_col=(((p1(2) >= Points(:,2) & p2(2) < Points(:,2)) | (p1(2) < Points(:,2) & p2(2) >= Points(:,2))) & ( Points(:,1) < (p2(1)-p1(1))*(Points(:,2)-p1(2)) / (p2(2)-p1(2))+p1(1))); 
                           in_xy(tmp_col)=~in_xy(tmp_col);
                       end
                      % in_z =  ( Points(:,3)>= ( Circle(5)  ) ) & ( Points(:,3)< ( Cells(5) + Cells(2) ) );   
                      in_z =  ( Points(:,3)>= ( Poly(5)  ) ) & ( Points(:,3)< ( Poly(5) + Poly(2) ) );       
                      In_out= in_xy' & in_z ;
                end
                 
                 
                %% Collision between Polynome and a list of points
                function   In_out=Poly_Points_Bounding_Box(Poly,Points)     
                     % Poly [Radius Length Pos_X Pos_Y Pos_Z %Surface %Volume Nb_Poly p1 p2 p3 p4 ..]
                        In_out =false(size(Points,1),1); 
                        In_out_z=false(size(Points,1),1);  
                        In_out_ob=false(size(Points,1),1);  
                        In_out_ib=false(size(Points,1),1); 
                        
                        tmp_Poly1=Poly;
                       
                        %%  Z solver 
                        In_out_z =  ( Points(:,3)>= ( Poly(5)  ) ) & ( Points(:,3)< ( Poly(5) + Poly(2) ) ); 
                        
                        
                        % Outerbounding box check
                        tmp_Poly1(1)=max(tmp_Poly1(9:end));
                        In_out_ob(In_out_z,:)=Collision_ToolBox.Circle_Points(tmp_Poly1,Points(In_out_z,:)); 
                        
                        % InnerBounding box check
                        tmp_Poly1(1)=min(tmp_Poly1(9:end));       
                        In_out_ib(In_out_ob)=Collision_ToolBox.Circle_Points(tmp_Poly1,Points(In_out_ob,:)); 
                        
                        In_out(In_out_ib)=true; % We are sure that these ones are in;
                        
                        % We test which remains
                        Idx_poly=In_out_ob&~In_out_ib;
                        
                        % Points that are outside the max box will never be in touching, check if it's a true collision
                        In_out(Idx_poly)=Collision_ToolBox.Poly_Points(Poly,Points(Idx_poly,:));
                end
                 
                
                 %% Collision between two Polynomes
                 function   In_out=Poly_Poly(Poly1,Poly2)            
                   % Poly [Radius Length Pos_X Pos_Y Pos_Z %Surface %Volume Nb_Poly p1 p2 p3 p4 ..]
                   In_out=false;
                   pp1=[];
                   pp3=[];
                   tmp_Poly1=Poly1;
                   tmp_Poly1(1)=max(tmp_Poly1(9:end));
                   tmp_Poly2=Poly2;
                   tmp_Poly2(1)=max(tmp_Poly2(9:end));
                   
                 % if in_z =  ( Points(:,3)>= ( Poly(5)  ) ) & ( Points(:,3)< ( Poly(5) + Poly(2) ) );   
                   if Collision_ToolBox.Circle_Circle(tmp_Poly1,tmp_Poly2) % Outerbounding box check
                        % OuterBounding box are touching, check if it's a true collision
                       tmp_Poly1(1)=min(tmp_Poly1(9:end));
                       tmp_Poly2(1)=min(tmp_Poly2(9:end));
                       if ~Collision_ToolBox.Circle_Circle(tmp_Poly1,tmp_Poly2) % InnerBounding box check
                           % InnerBounding box are not touching, check if there is a vertice collision 
                           for cpt_pol=1:1:Poly1(8)
                               p1(:,1)=Poly1(3)+Poly1(8+cpt_pol)*(cos((2*pi*(cpt_pol-1))/(Poly1(8))));
                               p1(:,2)=Poly1(4)+Poly1(8+cpt_pol)*(sin((2*pi*(cpt_pol-1))/(Poly1(8))));
                               p3(:,1)=Poly2(3)+Poly2(8+cpt_pol)*(cos((2*pi*(cpt_pol-1))/(Poly2(8))));
                               p3(:,2)=Poly2(4)+Poly2(8+cpt_pol)*(sin((2*pi*(cpt_pol-1))/(Poly2(8))));
                               if cpt_pol==Poly1(8)
                                    p2(1)=Poly1(3)+Poly1(8+1)*cos(0);
                                    p2(2)=Poly1(4)+Poly1(8+1)*sin(0);
                                    p4(1)=Poly2(3)+Poly2(8+1)*cos(0);
                                    p4(2)=Poly2(4)+Poly2(8+1)*sin(0);
                               else
                                    p2(1)=Poly1(3)+Poly1(8+cpt_pol+1)*(cos((2*pi*(cpt_pol))/(Poly1(8))));
                                    p2(2)=Poly1(4)+Poly1(8+cpt_pol+1)*(sin((2*pi*(cpt_pol))/(Poly1(8))));
                                    p4(1)=Poly2(3)+Poly2(8+cpt_pol+1)*(cos((2*pi*(cpt_pol))/(Poly2(8))));
                                    p4(2)=Poly2(4)+Poly2(8+cpt_pol+1)*(sin((2*pi*(cpt_pol))/(Poly2(8))));
                               end
                               
                               if Collision_ToolBox.Poly_Line(Poly2,p1,p2)
                                   In_out=true; % Vertice collision
                                   break;
                               end
                           end
                       else
                           In_out=true; % Innerbounding box are Touching, this is a collision collide
                       end
                   else
                       In_out=false; % Outerbounding box not Touching, they will never collide
                   end   
                 end

                %% Collision between Polynome and a Line
                function   In_out=Poly_Line(Poly,p3,p4)
                     % Poly [Radius Length Pos_X Pos_Y Pos_Z %Surface %Volume Nb_Poly p1 p2 p3 p4 ..]
                        In_out=false;
                           for cpt_pol=1:1:Poly(8)
                               p1(:,1)=Poly(3)+Poly(8+cpt_pol)*(cos((2*pi*(cpt_pol-1))/(Poly(8))));
                               p1(:,2)=Poly(4)+Poly(8+cpt_pol)*(sin((2*pi*(cpt_pol-1))/(Poly(8))));
                               if cpt_pol==Poly(8)
                                    p2(1)=Poly(3)+Poly(8+1)*cos(0);
                                    p2(2)=Poly(4)+Poly(8+1)*sin(0);
                               else
                                    p2(1)=Poly(3)+Poly(8+cpt_pol+1)*(cos((2*pi*(cpt_pol))/(Poly(8))));
                                    p2(2)=Poly(4)+Poly(8+cpt_pol+1)*(sin((2*pi*(cpt_pol))/(Poly(8))));
                               end
                               if Collision_ToolBox.Line_Line(p1,p2,p3,p4)
                                  In_out=true;
                                  break
                               end
                           end
                end

                %% Collision between a Line and a Line
                function   In_out=Line_Line(p1,p2,p3,p4)
                     uA = ((p4(1)-p3(1))*(p1(2)-p3(2)) - (p4(2)-p3(2))*(p1(1)-p3(1))) / ((p4(2)-p3(2))*(p2(1)-p1(1)) - (p4(1)-p3(1))*(p2(2)-p1(2)));
                     uB = ((p2(1)-p1(1))*(p1(2)-p3(2)) - (p2(2)-p1(2))*(p1(1)-p3(1))) / ((p4(2)-p3(2))*(p2(1)-p1(1)) - (p4(1)-p3(1))*(p2(2)-p1(2))); 
                     %   // if uA and uB are between 0-1, lines are colliding
                    if (uA >= 0 && uA <= 1 && uB >= 0 && uB <= 1)    
                          In_out=true;
                    else
                          In_out=false;
                    end
                end

                %% Collision between two circles
                function   In_out=Circle_Circle(Circle1,Circle2)
                    In_out=false;
                    a = Circle1(3) - Circle2(3);
                    b = Circle1(4) - Circle2(4);
                    dist=sqrt(a*a+b*b);
                    if dist < (Circle1(1) + Circle2(1)) % distance too small = collision
                        In_out=true;
                    end
                end
               
                %% Collision inside a ROI for a list of points
                function   In_out=ROI_Points(ROI,Points)     
                     % Poly [Radius Length Pos_X Pos_Y Pos_Z %Surface %Volume Nb_Poly p1 p2 p3 p4 ..]
                       in_xy =false(1,size(Points,1));
                      for cpt_cell=1:1:size(ROI,1)
                           p1(:,1)=ROI(cpt_cell,1);
                           p1(:,2)=ROI(cpt_cell,2);
                           if cpt_cell<size(ROI,1)
                                p2(:,1)=ROI(cpt_cell+1,1);
                                p2(:,2)=ROI(cpt_cell+1,2);
                           else
                               p2(:,1)=ROI(1,1);
                               p2(:,2)=ROI(1,2);
                           end
                           tmp_col=(((p1(2) >= Points(:,2) & p2(2) < Points(:,2)) | (p1(2) < Points(:,2) & p2(2) >= Points(:,2))) & ( Points(:,1) < (p2(1)-p1(1))*(Points(:,2)-p1(2)) / (p2(2)-p1(2))+p1(1))); 
                           in_xy(tmp_col)=~in_xy(tmp_col);
                       end 
                      In_out= in_xy ;
                 end
    end
end