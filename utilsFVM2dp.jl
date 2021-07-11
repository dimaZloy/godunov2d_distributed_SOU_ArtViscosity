
# utilities for FVM 


@everywhere @inline function phs2dcns2dcellsSA(
	ACons::SharedArray{Float64,2}, testFields::fields2d_shared, gamma::Float64)

	N::Int64 = size(testFields.densityCells,1);
	
	#ACons = zeros(Float64,N,4);

	for i = 1:N
		ACons[i,1] = testFields.densityCells[i];
		ACons[i,2] = testFields.densityCells[i]*testFields.UxCells[i];
		ACons[i,3] = testFields.densityCells[i]*testFields.UyCells[i];
		ACons[i,4] = testFields.pressureCells[i]/(gamma-1.0) + 0.5*testFields.densityCells[i]*(	testFields.UxCells[i]*testFields.UxCells[i] +  testFields.UyCells[i]*testFields.UyCells[i] );

	end #for
	
	#return ACons;
end




@everywhere  @inline function cells2nodesSolutionReconstructionWithStencilsSA(
		testMesh::mesh2d,cell_solution::SharedArray{Float64,1} ) ::Array{Float64,1}

node_solution = zeros(Float64,testMesh.nNodes); 

for J=1:testMesh.nNodes
	det::Float64 = 0.0;
	for j = 1:testMesh.nNeibCells
		neibCell::Int64 = testMesh.cell_clusters[J,j]; 
		if (neibCell !=0)
			wi::Float64 = testMesh.node_stencils[J,j];
			node_solution[J] += cell_solution[neibCell]*wi;
			det += wi;
		end
	end
	if (det!=0)
		node_solution[J] = node_solution[J]/det; 
	end
end

return node_solution;	

end




@everywhere  @inline function cells2nodesSolutionReconstructionWithStencilsUConsSA(nodesThreadsX:: SharedArray{Int64,2}, testMesh::mesh2d_shared, 
	cell_solution::SharedArray{Float64,2}, node_solution::SharedArray{Float64,2} )

	
	
		@sync @distributed for p in workers()	
		
			beginNode::Int64 = nodesThreadsX[p-1,1];
			endNode::Int64 = nodesThreadsX[p-1,2];
			
			for J=beginNode:endNode
			
				det::Float64 = 0.0;
				
				nNeibCells = size(testMesh.cell_clusters,2);
				
				node_solution[J,1] = 0.0;
				node_solution[J,2] = 0.0;
				node_solution[J,3] = 0.0;
				node_solution[J,4] = 0.0;
				
				for j = 1:nNeibCells
				
					neibCell::Int64 = testMesh.cell_clusters[J,j]; 
					
					if (neibCell !=0)
						wi::Float64 = testMesh.node_stencils[J,j];
						node_solution[J,1] += cell_solution[neibCell,1]*wi;
						node_solution[J,2] += cell_solution[neibCell,2]*wi;
						node_solution[J,3] += cell_solution[neibCell,3]*wi;
						node_solution[J,4] += cell_solution[neibCell,4]*wi;
					
						det += wi;
					end
				end
				
				if (det!=0)
					node_solution[J,1] = node_solution[J,1]/det; 
					node_solution[J,2] = node_solution[J,2]/det; 
					node_solution[J,3] = node_solution[J,3]/det; 
					node_solution[J,4] = node_solution[J,4]/det; 
				end
				
				
			end ## for

		

		end ## p workers 


end


@everywhere function cells2nodesSolutionReconstructionWithStencilsImplicitSA(nodesThreadsX:: SharedArray{Int64,2}, 
	testMesh::mesh2d_shared, testFields::fields2d_shared, dummy::SharedArray{Float64,2})	

		#nNodes = size(testMesh.cell_clusters,1);
		#node_solution = SharedArray{Float64}(nNodes,4); 
	
		@sync @distributed for p in workers()	
		
			beginNode::Int64 = nodesThreadsX[p-1,1];
			endNode::Int64 = nodesThreadsX[p-1,2];
	
			
			for J=beginNode:endNode
	
				dummy[J,1] = 0.0;
				dummy[J,2] = 0.0;
				dummy[J,3] = 0.0;
				dummy[J,4] = 0.0;
	
				det::Float64 = 0.0;
				nNeibCells = size(testMesh.node_stencils,2);
				
				for j = 1:nNeibCells
		
					neibCell::Int64 = testMesh.cell_clusters[J,j]; 
			
					if (neibCell !=0)
						wi::Float64 = testMesh.node_stencils[J,j];
						#node_solution[J,:] += cell_solution[neibCell,:];
						dummy[J,1] += testFields.densityCells[neibCell]*wi;
						dummy[J,2] += testFields.UxCells[neibCell]*wi;
						dummy[J,3] += testFields.UyCells[neibCell]*wi;
						dummy[J,4] += testFields.pressureCells[neibCell]*wi;
				 
						det += wi;
					end
				end
				
				if (det!=0)
					dummy[J,1] = dummy[J,1]/det; 
					dummy[J,2] = dummy[J,2]/det; 
					dummy[J,3] = dummy[J,3]/det; 
					dummy[J,4] = dummy[J,4]/det; 
				end
				
			end

	
			for J = beginNode:endNode
	
				testFields.densityNodes[J] = dummy[J,1]; 
				testFields.UxNodes[J] 	   = dummy[J,2]; 
				testFields.UyNodes[J] 	   = dummy[J,3]; 
				testFields.pressureNodes[J] =  dummy[J,4]; 
		
			end
	
	
	
	
	
		end ## p workers
	

end