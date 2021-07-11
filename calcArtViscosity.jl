

@everywhere @inline function calcArtificialViscositySA( cellsThreadsX::SharedArray{Int64,2}, testMesh::mesh2d_shared, 
	testfields2d::fields2d_shared, viscfields2dX::viscousFields2d_shared)
	
	
	@sync @distributed for p in workers()	
	
		beginCell::Int64 = cellsThreadsX[p-1,1];
		endCell::Int64 = cellsThreadsX[p-1,2];
		
		
		
		nodesGradientReconstructionFastPerThread(beginCell, endCell, testMesh, testfields2dX.UxNodes, viscfields2dX.dUdxCells,viscfields2dX.dUdyCells);
		nodesGradientReconstructionFastPerThread(beginCell, endCell, testMesh, testfields2dX.UyNodes, viscfields2dX.dVdxCells,viscfields2dX.dVdyCells);

		calcArtificialViscosityPerThread( beginCell, endCell, testMesh, testfields2d, viscfields2dX);
					
	end
	
	
end




@everywhere @inline function calcArtificialViscosityPerThread( beginCell::Int64, endCell::Int64, 
	testMesh::mesh2d_shared, testfields2d::fields2d_shared, viscfields2dX::viscousFields2d_shared)



     for i = beginCell:endCell
	 
		
		divU::Float64 = viscfields2dX.dUdxCells[i] + viscfields2dX.dVdyCells[i];
	 
     
        Cth::Float64 = 0.05;
        Cav::Float64 = 0.5;
        h::Float64 = testMesh.HX[i]/sqrt(2.0);
		a::Float64 = Cth*testfields2d.aSoundCells[i]/h; 
		
         if (-divU >  a)
         
             ## tmp = divU[i]*divU[i] - (Cth*aSound[i]/h)*(Cth*aSound[i]/h);
			 tmp::Float64 = divU*divU - a*a;
             PSI::Float64 = 1.0e-5;
			 
			 ##PSI::Float64 = 1.0e-5;
			 
             if (tmp > 0.0)
                viscfields2dX.artViscosityCells[i] = Cav*testfields2d.densityCells[i]*h*h*sqrt(tmp)*PSI;
			 end

         end 


     end ## for
	 
	 

end