

@everywhere function SecondOrderUpwindM2(
	beginCell::Int32,endCell::Int32, bettaKJ::Float64, dt::Float64, flowTime::Float64, 
	testMesh::mesh2d_shared, testFields::fields2d_shared, thermo::THERMOPHYSICS, 
	UConsCellsOld::SharedArray{Float64,2}, FLUXES::SharedArray{Float64,2}, UconsDiffTerm::SharedArray{Float64,2}, UconsCellsNew::SharedArray{Float64,2})

	#nCells = size(UConsCellsOld,1);
	
	uLeftp = zeros(Float64,4);
	#FLUXES = zeros(Float64,nCells);
	
	for i = beginCell:endCell
    
		num_nodes::Int64 = testMesh.mesh_connectivity[i,3];
		
		
		uLeftp[1] = testFields.densityCells[i];
		uLeftp[2] = testFields.UxCells[i];
		uLeftp[3] = testFields.UyCells[i];
		uLeftp[4] = testFields.pressureCells[i];
		
	   
		
		if (num_nodes == 3)
		
			edge_flux1 = ( computeInterfaceSlope(i, Int32(1), testMesh, testFields, thermo, uLeftp, flowTime) );
			edge_flux2 = ( computeInterfaceSlope(i, Int32(2), testMesh, testFields, thermo, uLeftp, flowTime) );
			edge_flux3 = ( computeInterfaceSlope(i, Int32(3), testMesh, testFields, thermo, uLeftp, flowTime) );
				
			FLUXES[i,1] = edge_flux1[1] + edge_flux2[1] + edge_flux3[1];
			FLUXES[i,2] = edge_flux1[2] + edge_flux2[2] + edge_flux3[2];
			FLUXES[i,3] = edge_flux1[3] + edge_flux2[3] + edge_flux3[3];
			FLUXES[i,4] = edge_flux1[4] + edge_flux2[4] + edge_flux3[4];
			

		elseif (num_nodes == 4)
			
			edge_flux1 = ( computeInterfaceSlope(i, Int32(1), testMesh, testFields, uLeftp, flowTime) );
			edge_flux2 = ( computeInterfaceSlope(i, Int32(2), testMesh, testFields, uLeftp, flowTime) );
			edge_flux3 = ( computeInterfaceSlope(i, Int32(3), testMesh, testFields, uLeftp, flowTime) );
			edge_flux4 = ( computeInterfaceSlope(i, Int32(4), testMesh, testFields, uLeftp, flowTime) );
				
			FLUXES[i,1] = edge_flux1[1] + edge_flux2[1] + edge_flux3[1] + edge_flux4[1];
			FLUXES[i,2] = edge_flux1[2] + edge_flux2[2] + edge_flux3[2] + edge_flux4[2];
			FLUXES[i,3] = edge_flux1[3] + edge_flux2[3] + edge_flux3[3] + edge_flux4[3];
			FLUXES[i,4] = edge_flux1[4] + edge_flux2[4] + edge_flux3[4] + edge_flux4[4];

		else
			
			display("something wrong in flux calculations ... ")
			
		end
		
		
		UconsCellsNew[i,1] = ( UConsCellsOld[i,1] - FLUXES[i,1]*bettaKJ*dt*testMesh.Z[i] + bettaKJ*dt*UconsDiffTerm[i,1] );
		UconsCellsNew[i,2] = ( UConsCellsOld[i,2] - FLUXES[i,2]*bettaKJ*dt*testMesh.Z[i] + bettaKJ*dt*UconsDiffTerm[i,2] );
		UconsCellsNew[i,3] = ( UConsCellsOld[i,3] - FLUXES[i,3]*bettaKJ*dt*testMesh.Z[i] + bettaKJ*dt*UconsDiffTerm[i,3] );
		UconsCellsNew[i,4] = ( UConsCellsOld[i,4] - FLUXES[i,4]*bettaKJ*dt*testMesh.Z[i] + bettaKJ*dt*UconsDiffTerm[i,4] );
      
   
	end # i - loop for all cells


end


	


