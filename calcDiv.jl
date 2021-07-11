

@everywhere function nodesDivergenceReconstructionSA( beginCell::Int64, endCell::Int64, testMesh::mesh2d_shared, 
  gradX::SharedArray{Float64,1}, gradY::SharedArray{Float64,1}, divergence::SharedArray{Float64,1})
  
  
  phiLeftX::Float64 = 0.0;
  phiLeftY::Float64 = 0.0;

  phiRightX::Float64 = 0.0;
  phiRightY::Float64 = 0.0;

  phiFaceX = zeros(Float64,4);
  phiFaceY  = zeros(Float64,4);


  for C = beginCell:endCell
  
	   

        phiFaceX[1] = phiFaceX[2] =  phiFaceX[3] = phiFaceX[4] = 0.0;	
		phiFaceY[1] = phiFaceY[2] =  phiFaceY[3] = phiFaceY[4] = 0.0;

		numNodesInCell::Int64 = testMesh.mesh_connectivity[C,3]; ## CMatrix mesh_connectivity - first index == 1

        phiLeftX = phiLeftY = phiRightX = phiRightY = 0.0;

        phiFaceX[1] = phiFaceX[2] = phiFaceX[3] = phiFaceX[4] = 0.0;
        phiFaceY[1] = phiFaceY[2] = phiFaceY[3] = phiFaceY[4] = 0.0;


       for T=1:numNodesInCell

           if (T == 1)
           
               phiLeftX = gradX[ testMesh.cells2nodes[C,1] ];
               phiLeftY = gradY[ testMesh.cells2nodes[C,1] ];

               phiRightX = gradX[ testMesh.cells2nodes[C,2] ];
               phiRightY = gradY[ testMesh.cells2nodes[C,2] ];

           
           elseif (T == 2)
           
               phiLeftX = gradX[ testMesh.cells2nodes[C,3] ];
               phiLeftY = gradY[ testMesh.cells2nodes[C,3] ];

               phiRightX = gradX[ testMesh.cells2nodes[C,4] ];
               phiRightY = gradY[ testMesh.cells2nodes[C,4] ];

           
           elseif (T == 3)
           

               phiLeftX = gradX[ testMesh.cells2nodes[C,5] ];
               phiLeftY = gradY[ testMesh.cells2nodes[C,5] ];

               phiRightX = gradX[ testMesh.cells2nodes[C,6] ];
               phiRightY = gradY[ testMesh.cells2nodes[C,6] ];

           
           elseif (T == 4)
           
               phiLeftX = gradX[ testMesh.cells2nodes[C,7] ];
               phiLeftY = gradY[ testMesh.cells2nodes[C,7] ];

               phiRightX = gradX[ testMesh.cells2nodes[C,8] ];
               phiRightY = gradY[ testMesh.cells2nodes[C,8] ];

           end


		   side::Float64 = testMesh.cell_edges_length[C,T];
           nx::Float64 = testMesh.cell_edges_Nx[C,T];
           ny::Float64 = testMesh.cell_edges_Ny[C,T];	

           phiFaceY[T] = 0.5*(phiLeftY + phiRightY)*-ny*side;
           phiFaceX[T] = 0.5*(phiLeftX + phiRightX)*-nx*side;

           # if (numNodesInCell == 3)
           
               # double divX = (phiFaceX[1] + phiFaceX[2] + phiFaceX[3])/areaCell;
               # double divY = (phiFaceY[1] + phiFaceY[2] + phiFaceY[3])/areaCell;

               # divergence[C] = divX + divY;
           
           # elseif (numNodes == 4)
           
               # double divX = (phiFaceX[1] + phiFaceX[2] + phiFaceX[3] + phiFaceX[4])/areaCell;
               # double divY = (phiFaceY[1] + phiFaceY[2] + phiFaceY[3] + phiFaceY[4])/areaCell;

               # divergence[C] = divX + divY;

           # end

       end ## end of nodes


		areaCell::Float64 = testMesh.cell_areas[C]; 

		divX::Float64 = (phiFaceX[1] + phiFaceX[2] + phiFaceX[3] + phiFaceX[4])/areaCell;
        divY::Float64 = (phiFaceY[1] + phiFaceY[2] + phiFaceY[3] + phiFaceY[4])/areaCell;


		divergence[C] = divX + divY;


  end ## end of cells

  

end ## end of function



@everywhere function nodesDivergenceReconstructionFastSA( beginCell::Int64, endCell::Int64, testMesh::mesh2d_shared, 
  gradX::SharedArray{Float64,1}, gradY::SharedArray{Float64,1}, divergence::SharedArray{Float64,1})
  
  
 

  for C = beginCell:endCell
  
		phiLeftX = zeros(Float64,4);
		phiLeftY = zeros(Float64,4);

		phiRightX = zeros(Float64,4);
		phiRightY = zeros(Float64,4);

		phiFaceX = zeros(Float64,4);
		phiFaceY  = zeros(Float64,4);

		side = zeros(Float64,4);
		nx = zeros(Float64,4);
		ny = zeros(Float64,4);

		numNodesInCell::Int64 = testMesh.mesh_connectivity[C,3]; ## CMatrix mesh_connectivity - first index == 1
		
		
		T1::Int64 = testMesh.mesh_connectivity[C,4];
		T2::Int64 = testMesh.mesh_connectivity[C,5];
		T3::Int64 = testMesh.mesh_connectivity[C,6];
	  
		side[1] = testMesh.cell_edges_length[C,1];
		side[2] = testMesh.cell_edges_length[C,2];
		side[3] = testMesh.cell_edges_length[C,3];
	  
		nx[1] = testMesh.cell_edges_Nx[C,1];
		nx[2] = testMesh.cell_edges_Nx[C,2];
		nx[3] = testMesh.cell_edges_Nx[C,3];
	  
		ny[1] = testMesh.cell_edges_Nx[C,1];
		ny[2] = testMesh.cell_edges_Nx[C,2];
		ny[3] = testMesh.cell_edges_Nx[C,3];
		

		phiLeftX[1] =  gradX[ testMesh.cells2nodes[C,1] ];
		phiLeftY[1] =  gradY[ testMesh.cells2nodes[C,1] ];
	
		phiRightX[1] = gradX[ testMesh.cells2nodes[C,2] ];
		phiRightY[1] = gradY[ testMesh.cells2nodes[C,2] ];

		phiLeftX[2] =  gradX[ testMesh.cells2nodes[C,3] ];
		phiLeftY[2] =  gradY[ testMesh.cells2nodes[C,3] ];
	
		phiRightX[2] = gradX[ testMesh.cells2nodes[C,4] ];
		phiRightY[2] = gradY[ testMesh.cells2nodes[C,4] ];
	
		phiLeftX[3] =  gradX[ testMesh.cells2nodes[C,5] ];
		phiLeftY[3] =  gradY[ testMesh.cells2nodes[C,5] ];
	
		phiRightX[3] = gradX[ testMesh.cells2nodes[C,6] ];
		phiRightY[3] = gradY[ testMesh.cells2nodes[C,6] ];


		if (numNodesInCell == 4)
	  
			T4::Int64 = testMesh.mesh_connectivity[C,7];
			side[4] = testMesh.cell_edges_length[C,4];
			nx[4] = testMesh.cell_edges_Nx[C,4];
			ny[4] = testMesh.cell_edges_Nx[C,4];


			phiLeftX[4] =  gradX[ testMesh.cells2nodes[C,7] ];
			phiLeftY[4] =  gradY[ testMesh.cells2nodes[C,7] ];
	
			phiRightX[4] = gradX[ testMesh.cells2nodes[C,8] ];
			phiRightY[4] = gradY[ testMesh.cells2nodes[C,8] ];

		
		end


		#phiFaceY[T] = 0.5*(phiLeftY + phiRightY)*-ny*side;
        #phiFaceX[T] = 0.5*(phiLeftX + phiRightX)*-nx*side;

		phiFaceX[1] = 0.5*(phiLeftX[1] + phiRightX[1])*-nx[1]*side[1];		
		phiFaceY[1] = 0.5*(phiLeftY[1] + phiRightY[1])*-ny[1]*side[1];

		phiFaceX[2] = 0.5*(phiLeftX[2] + phiRightX[2])*-nx[2]*side[2];		
		phiFaceY[2] = 0.5*(phiLeftY[2] + phiRightY[2])*-ny[2]*side[2];
	  
		phiFaceX[3] = 0.5*(phiLeftX[3] + phiRightX[3])*-nx[3]*side[3];		
		phiFaceY[3] = 0.5*(phiLeftY[3] + phiRightY[3])*-ny[3]*side[3];
	  
		phiFaceX[4] = 0.5*(phiLeftX[4] + phiRightX[4])*-nx[4]*side[4];		
		phiFaceY[4] = 0.5*(phiLeftY[4] + phiRightY[4])*-ny[4]*side[4]; 	


		areaCell::Float64 = testMesh.cell_areas[C]; 

		divX::Float64 = (phiFaceX[1] + phiFaceX[2] + phiFaceX[3] + phiFaceX[4])/areaCell;
        divY::Float64 = (phiFaceY[1] + phiFaceY[2] + phiFaceY[3] + phiFaceY[4])/areaCell;


		divergence[C] = divX + divY;


  end ## end of cells

  

end ## end of function