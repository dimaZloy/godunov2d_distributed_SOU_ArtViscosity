

function calcDiffTerm(cellsThreadsX::SharedArray{Int64,2}, testMeshDistrX::mesh2d_shared, 
		testfields2dX::fields2d_shared, viscousFields2dX::viscousFields2d_shared, thermoX::THERMOPHYSICS, 
		UCNodes::SharedArray{Float64,2}, UConsDiffCells::SharedArray{Float64,2}, UConsDiffNodes::SharedArray{Float64,2})



		@sync @distributed for p in workers()	
	
			beginCell::Int64 = cellsThreadsX[p-1,1];
			endCell::Int64 = cellsThreadsX[p-1,2];
			
			nodesGradientReconstructionUconsFastSA(beginCell,endCell, testMeshDistrX, UCNodes, viscousFields2dX);	
			
		end
		
		
		
		@sync @distributed for p in workers()	
	
			beginCell::Int64 = cellsThreadsX[p-1,1];
			endCell::Int64 = cellsThreadsX[p-1,2];
			
			nodesDivergenceReconstructionFastSA(beginCell,endCell, testMeshDistrX,  viscousFields2dX.cdUdxCells,viscousFields2dX.cdUdyCells, viscousFields2dX.laplasUCuCells);
			nodesDivergenceReconstructionFastSA(beginCell,endCell, testMeshDistrX,  viscousFields2dX.cdVdxCells,viscousFields2dX.cdVdyCells, viscousFields2dX.laplasUCvCells);
			nodesDivergenceReconstructionFastSA(beginCell,endCell, testMeshDistrX,  viscousFields2dX.cdEdxCells,viscousFields2dX.cdEdyCells, viscousFields2dX.laplasUCeCells);

			
			
		end
		
		
		##display(viscousFields2dX.laplasUCeCells)

		
		@sync @distributed for p in workers()	
	
			beginCell::Int64 = cellsThreadsX[p-1,1];
			endCell::Int64 = cellsThreadsX[p-1,2];
			
		
			for i = beginCell:endCell
			
				avEpsilon::Float64 = 1e-5;
					
				if (viscousFields2dX.artViscosityCells[i] > avEpsilon)
					
					Pr::Float64 = 3.0/4.0;
					UConsDiffCells[i,2] = viscousFields2dX.artViscosityCells[i]*viscousFields2dX.laplasUCuCells[i];
					UConsDiffCells[i,3] = viscousFields2dX.artViscosityCells[i]*viscousFields2dX.laplasUCvCells[i];
					UConsDiffCells[i,4] = viscousFields2dX.artViscosityCells[i]*thermoX.Cp/Pr*viscousFields2dX.laplasUCeCells[i];

				end

			end
		
		
					
		end
		


end
