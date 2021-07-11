

@everywhere function nodesGradientReconstructionFastSA(beginCell::Int64, endCell::Int64, testMesh::mesh2d_shared, 
  scalarField::SharedArray{Float64,1}, gradX::SharedArray{Float64,1}, gradY::SharedArray{Float64,1}, gradMag::SharedArray{Float64,1})


  for C =  beginCell:endCell
  
  
	phiFaceX = zeros(Float64,4);
	phiFaceY = zeros(Float64,4);
	
	phiLeft = zeros(Float64,4);
	phiRight = zeros(Float64,4);
	
	side = zeros(Float64,4);
 	nx = zeros(Float64,4);
    ny = zeros(Float64,4);


      numNodesInCell::Int64 = testMesh.mesh_connectivity[C,3]; ## CMatrix mesh_connectivity - first index == 1
	  
	  areaCell::Float64 = testMesh.cell_areas[C]; 
	  
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
	  
	  
	  
	  

	  phiLeft[1] =  scalarField[ testMesh.cells2nodes[C,1] ];
      phiRight[1] = scalarField[ testMesh.cells2nodes[C,2] ];

	  phiLeft[2] =  scalarField[ testMesh.cells2nodes[C,3] ];
      phiRight[2] = scalarField[ testMesh.cells2nodes[C,4] ];

	  phiLeft[3] =  scalarField[ testMesh.cells2nodes[C,5] ];
      phiRight[3] = scalarField[ testMesh.cells2nodes[C,6] ];
	  

	  if (numNodesInCell == 4)
	  
	    T4::Int64 = testMesh.mesh_connectivity[C,7];
		side[4] = testMesh.cell_edges_length[C,4];
		nx[4] = testMesh.cell_edges_Nx[C,4];
		ny[4] = testMesh.cell_edges_Nx[C,4];

		phiLeft[4] =  scalarField[ testMesh.cells2nodes[C,7] ];
		phiRight[4] = scalarField[ testMesh.cells2nodes[C,8] ];
		
	  end
	  
	  
	  phiFaceX[1] = 0.5*(phiLeft[1] + phiRight[1])*-nx[1]*side[1];		
      phiFaceY[1] = 0.5*(phiLeft[1] + phiRight[1])*-ny[1]*side[1];

	  phiFaceX[2] = 0.5*(phiLeft[2] + phiRight[2])*-nx[2]*side[2];		
      phiFaceY[2] = 0.5*(phiLeft[2] + phiRight[2])*-ny[2]*side[2];
	  
	  phiFaceX[3] = 0.5*(phiLeft[3] + phiRight[3])*-nx[3]*side[3];		
      phiFaceY[3] = 0.5*(phiLeft[3] + phiRight[3])*-ny[3]*side[3];
	  
	  phiFaceX[4] = 0.5*(phiLeft[4] + phiRight[4])*-nx[4]*side[4];		
      phiFaceY[4] = 0.5*(phiLeft[4] + phiRight[4])*-ny[4]*side[4];


      gradX[C] =  (phiFaceX[1] + phiFaceX[2] + phiFaceX[3] + phiFaceX[4]) / areaCell;
      gradY[C] =  (phiFaceY[1] + phiFaceY[2] + phiFaceY[3] + phiFaceY[4]) / areaCell;
	  
	  gradMag[C] = sqrt(gradX[C]*gradX[C] + gradY[C]*gradY[C]);



  end ## end global loop for cells



end



# @everywhere @inline function calcVelocityDerivativesSA( cellsThreadsX::SharedArray{Int64,2}, testMesh::mesh2d_shared, 
	# testFields::fields2d_shared, viscousFields2dX::viscousFields2d_shared)
	
	
	# @sync @distributed for p in workers()	
	
		# beginCell::Int64 = cellsThreadsX[p-1,1];
		# endCell::Int64 = cellsThreadsX[p-1,2];

		# nodesGradientReconstructionFastPerThread(beginCell, endCell, testMesh, testfields2dX.UxNodes, viscfields2dX.dUdxCells,viscfields2dX.dUdyCells);
		# nodesGradientReconstructionFastPerThread(beginCell, endCell, testMesh, testfields2dX.UyNodes, viscfields2dX.dVdxCells,viscfields2dX.dVdyCells);
					
	# end
	
	
# end


@everywhere function nodesGradientReconstructionFastPerThread(beginCell::Int64, endCell::Int64, testMesh::mesh2d_shared, 
  scalarField::SharedArray{Float64,1}, gradX::SharedArray{Float64,1}, gradY::SharedArray{Float64,1})


  for C =  beginCell:endCell
  
  
	phiFaceX = zeros(Float64,4);
	phiFaceY = zeros(Float64,4);
	
	phiLeft = zeros(Float64,4);
	phiRight = zeros(Float64,4);
	
	side = zeros(Float64,4);
 	nx = zeros(Float64,4);
    ny = zeros(Float64,4);


      numNodesInCell::Int64 = testMesh.mesh_connectivity[C,3]; ## CMatrix mesh_connectivity - first index == 1
	  
	  areaCell::Float64 = testMesh.cell_areas[C]; 
	  
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
	  
	  
	  
	  

	  phiLeft[1] =  scalarField[ testMesh.cells2nodes[C,1] ];
      phiRight[1] = scalarField[ testMesh.cells2nodes[C,2] ];

	  phiLeft[2] =  scalarField[ testMesh.cells2nodes[C,3] ];
      phiRight[2] = scalarField[ testMesh.cells2nodes[C,4] ];

	  phiLeft[3] =  scalarField[ testMesh.cells2nodes[C,5] ];
      phiRight[3] = scalarField[ testMesh.cells2nodes[C,6] ];
	  

	  if (numNodesInCell == 4)
	  
	    T4::Int64 = testMesh.mesh_connectivity[C,7];
		side[4] = testMesh.cell_edges_length[C,4];
		nx[4] = testMesh.cell_edges_Nx[C,4];
		ny[4] = testMesh.cell_edges_Nx[C,4];

		phiLeft[4] =  scalarField[ testMesh.cells2nodes[C,7] ];
		phiRight[4] = scalarField[ testMesh.cells2nodes[C,8] ];
		
	  end
	  
	  
	  phiFaceX[1] = 0.5*(phiLeft[1] + phiRight[1])*-nx[1]*side[1];		
      phiFaceY[1] = 0.5*(phiLeft[1] + phiRight[1])*-ny[1]*side[1];

	  phiFaceX[2] = 0.5*(phiLeft[2] + phiRight[2])*-nx[2]*side[2];		
      phiFaceY[2] = 0.5*(phiLeft[2] + phiRight[2])*-ny[2]*side[2];
	  
	  phiFaceX[3] = 0.5*(phiLeft[3] + phiRight[3])*-nx[3]*side[3];		
      phiFaceY[3] = 0.5*(phiLeft[3] + phiRight[3])*-ny[3]*side[3];
	  
	  phiFaceX[4] = 0.5*(phiLeft[4] + phiRight[4])*-nx[4]*side[4];		
      phiFaceY[4] = 0.5*(phiLeft[4] + phiRight[4])*-ny[4]*side[4];


      gradX[C] =  (phiFaceX[1] + phiFaceX[2] + phiFaceX[3] + phiFaceX[4]) / areaCell;
      gradY[C] =  (phiFaceY[1] + phiFaceY[2] + phiFaceY[3] + phiFaceY[4]) / areaCell;
	  
	 ## gradMag[C] = sqrt(gradX[C]*gradX[C] + gradY[C]*gradY[C]);



  end ## end global loop for cells



end



@everywhere function nodesGradientReconstructionUconsFastSA(beginCell::Int64, endCell::Int64, testMesh::mesh2d_shared, 
  ucons::SharedArray{Float64,2}, viscousFields2dX::viscousFields2d_shared)


  for V = 2:4


	  for C =  beginCell:endCell
	  
	  
		phiFaceX = zeros(Float64,4);
		phiFaceY = zeros(Float64,4);
		
		phiLeft = zeros(Float64,4);
		phiRight = zeros(Float64,4);
		
		side = zeros(Float64,4);
		nx = zeros(Float64,4);
		ny = zeros(Float64,4);


		  numNodesInCell::Int64 = testMesh.mesh_connectivity[C,3]; ## CMatrix mesh_connectivity - first index == 1
		  
		  areaCell::Float64 = testMesh.cell_areas[C]; 
		  
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
		  
		  
		  
		  

		  phiLeft[1] =  ucons[ testMesh.cells2nodes[C,1] ,V ];
		  phiRight[1] = ucons[ testMesh.cells2nodes[C,2] ,V ];

		  phiLeft[2] =  ucons[ testMesh.cells2nodes[C,3] ,V ];
		  phiRight[2] = ucons[ testMesh.cells2nodes[C,4] ,V ];

		  phiLeft[3] =  ucons[ testMesh.cells2nodes[C,5] ,V ];
		  phiRight[3] = ucons[ testMesh.cells2nodes[C,6] ,V ];
		  

		  if (numNodesInCell == 4)
		  
			T4::Int64 = testMesh.mesh_connectivity[C,7];
			side[4] = testMesh.cell_edges_length[C,4];
			nx[4] = testMesh.cell_edges_Nx[C,4];
			ny[4] = testMesh.cell_edges_Nx[C,4];

			phiLeft[4] =  ucons[ testMesh.cells2nodes[C,7] ,V  ];
			phiRight[4] = ucons[ testMesh.cells2nodes[C,8] ,V ];
			
		  end
		  
		  
		  phiFaceX[1] = 0.5*(phiLeft[1] + phiRight[1])*-nx[1]*side[1];		
		  phiFaceY[1] = 0.5*(phiLeft[1] + phiRight[1])*-ny[1]*side[1];

		  phiFaceX[2] = 0.5*(phiLeft[2] + phiRight[2])*-nx[2]*side[2];		
		  phiFaceY[2] = 0.5*(phiLeft[2] + phiRight[2])*-ny[2]*side[2];
		  
		  phiFaceX[3] = 0.5*(phiLeft[3] + phiRight[3])*-nx[3]*side[3];		
		  phiFaceY[3] = 0.5*(phiLeft[3] + phiRight[3])*-ny[3]*side[3];
		  
		  phiFaceX[4] = 0.5*(phiLeft[4] + phiRight[4])*-nx[4]*side[4];		
		  phiFaceY[4] = 0.5*(phiLeft[4] + phiRight[4])*-ny[4]*side[4];


		if (V == 2)
		  viscousFields2dX.cdUdxCells[C] =  (phiFaceX[1] + phiFaceX[2] + phiFaceX[3] + phiFaceX[4]) / areaCell;
		  viscousFields2dX.cdUdyCells[C] =  (phiFaceY[1] + phiFaceY[2] + phiFaceY[3] + phiFaceY[4]) / areaCell;
		elseif (V == 3)
		  viscousFields2dX.cdVdxCells[C] =  (phiFaceX[1] + phiFaceX[2] + phiFaceX[3] + phiFaceX[4]) / areaCell;
		  viscousFields2dX.cdVdyCells[C] =  (phiFaceY[1] + phiFaceY[2] + phiFaceY[3] + phiFaceY[4]) / areaCell;
		elseif (V == 4)
		
		  viscousFields2dX.cdEdxCells[C] =  (phiFaceX[1] + phiFaceX[2] + phiFaceX[3] + phiFaceX[4]) / areaCell;
		  viscousFields2dX.cdEdyCells[C] =  (phiFaceY[1] + phiFaceY[2] + phiFaceY[3] + phiFaceY[4]) / areaCell;
		
		  
		end


	  end ## end global loop for cells





	end ## for UCons[:,2], UCons[:,3], UCons[:,4]; skipped for UCons[:,1]!!!  


end

