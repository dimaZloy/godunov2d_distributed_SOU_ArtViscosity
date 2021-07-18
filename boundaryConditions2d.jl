
@everywhere @inline function ComputeUPhysFromBoundaries(i::Int32,k::Int32,neib_cell::Int32, 
      cur_cell::Array{Float64,1}, nx::Float64,ny::Float64, y::Float64, gamma::Float64, t::Float64 )::Array{Float64,1}

		bnd_cell = zeros(Float64,4);

		if (neib_cell == -3) #inlet

            bnd_cell[1] = 1.4;
            bnd_cell[2] = 300.00;
            bnd_cell[3] = 0.0;
            bnd_cell[4] = 10000.0;

		elseif (neib_cell == -2) #walls


			bnd_cell = updateVelocityFromCurvWall(i,k,cur_cell,nx,ny);

		elseif (neib_cell == -1) # outlet

			bnd_cell = cur_cell;	
					
		end	

	return bnd_cell; 
end


@everywhere @inline function updateVelocityFromCurvWall(i::Int32, k::Int32, U::Array{Float64,1}, nx::Float64, ny::Float64)

# High-Order Accurate Implementation of Solid Wall Boundary Conditions in Curved Geometries, 
# Lilia Krivodonova and Marsha Berger, Courant Institute of Mathematical Sciences, New York, NY 10012

# a = U[1]*(ny*ny - nx*nx) - 2.0*nx*ny*U[2];
# b = U[2]*(nx*nx - ny*ny) - 2.0*nx*ny*U[1];


	Un = deepcopy(U); 

        Un[2] = U[2]*(ny*ny - nx*nx) - 2.0*nx*ny*U[3];
        Un[3] = U[3]*(nx*nx - ny*ny) - 2.0*nx*ny*U[2];


	return Un;	
end

