

using Distributed;
using PyPlot;

const numThreads = 2;


if (numThreads != 1)

	if (nprocs() == 1)
		addprocs(numThreads,lazy=false); 
		display(workers());
	end
		
end


@everywhere using PyPlot;
@everywhere using WriteVTK;
@everywhere using CPUTime;
@everywhere using DelimitedFiles;
@everywhere using Printf
@everywhere using BSON: @load
@everywhere using BSON: @save
@everywhere using SharedArrays;

using HDF5;
using ProfileView;


include("primeObjects.jl");
include("thermo.jl"); #setup thermodynamics
include("utilsIO.jl");
include("RoeFlux2d.jl")
include("AUSMflux2d.jl"); #AUSM+ inviscid flux calculation 
include("utilsFVM2dp.jl"); #FVM utililities
## utilsFVM2dp::cells2nodesSolutionReconstructionWithStencilsImplicitSA
## utilsFVM2dp::cells2nodesSolutionReconstructionWithStencilsSA
## utilsFVM2dp::phs2dcns2dcellsSA

#include("calcCells2NodesMatrix.jl");
include("calcGrad.jl");
include("calcDiv.jl");
include("calcArtViscosity.jl");
include("calcDiffterm.jl");

include("boundaryConditions2d.jl"); 

include("initfields2d.jl");
## initfields2d::distibuteCellsInThreadsSA()
## initfields2d::createFields2d_shared()

include("evaluate2d.jl"); 
## propagate2d::updateResidualSA()
## propagate2d::updateVariablesSA()
## propagate2d::updateOutputSA()

include("limiters.jl");
include("computeslope2d.jl");
include("SOUscheme.jl");
## computeslope2d:: computeInterfaceSlope()
## SOUscheme:: SecondOrderUpwindM2()

include("propagate2d.jl");
## propagate:: calcOneStage() expilict Euler first order
## propagate:: doExplicitRK3TVD() expilict RK3-TVD






function godunov2dthreads(pname::String, numThreads::Int64, outputfile::String, coldrun::Bool)


	flag2loadPreviousResults = true;

	testMesh = readMesh2dHDF5(pname);
		
	cellsThreads = distibuteCellsInThreadsSA(numThreads, testMesh.nCells); ## partition mesh 
	nodesThreads = distibuteNodesInThreadsSA(numThreads, testMesh.nNodes); ## partition mesh 
	

	include("setupSolver2d.jl"); #setup FVM and numerical schemes
	
	
	## init primitive variables 
	println("set initial and boundary conditions ...");
	
	testfields2d = createFields2d_shared(testMesh, thermo);
	solInst = solutionCellsT(
		0.0,
		0.0,
		testMesh.nCells,
		testfields2d.densityCells,
		testfields2d.UxCells,
		testfields2d.UyCells,
		testfields2d.pressureCells,
	);
	
	
	#(testfields2d, solInst) = createFields2dLoadPrevResults_shared(testMesh, thermo, "zzz13700", dynControls);
	
	
	
	viscfields2d = createViscousFields2d_shared(testMesh.nCells, testMesh.nNodes);
	
	println("nCells:\t", testMesh.nCells);
	println("nNodes:\t", testMesh.nNodes);
	
	## init conservative variables 
	UconsCellsOldX = SharedArray{Float64}(testMesh.nCells,4);
	UconsNodesOldX = SharedArray{Float64}(testMesh.nNodes,4);
	UconsCellsNewX = SharedArray{Float64}(testMesh.nCells,4);
	
	# UconsCellsNew1X = SharedArray{Float64}(testMesh.nCells,4);
	# UconsCellsNew2X = SharedArray{Float64}(testMesh.nCells,4);
	# UconsCellsNew3X = SharedArray{Float64}(testMesh.nCells,4);
		
	
	UConsDiffCellsX = SharedArray{Float64}(testMesh.nCells,4);
	UConsDiffNodesX = SharedArray{Float64}(testMesh.nNodes,4);
	
	artViscosityNodes =	zeros(Float64,testMesh.nNodes);
	
	DeltaX = SharedArray{Float64}(testMesh.nCells,4);
	iFLUXX  = SharedArray{Float64}(testMesh.nCells,4);
	dummy  = SharedArray{Float64}(testMesh.nNodes,4);
	
	
	phs2dcns2dcellsSA(UconsCellsOldX,testfields2d, thermo.Gamma);	
	phs2dcns2dcellsSA(UconsCellsNewX,testfields2d, thermo.Gamma);	
	
	cells2nodesSolutionReconstructionWithStencilsSerialUCons(testMesh, UconsCellsOldX,  UconsNodesOldX );	

	
	
	testMeshDistr = createMesh2dShared(testMesh);
	
	
	#@everywhere trianglesX = $triangles;
	@everywhere testMeshDistrX = $testMeshDistr; 
	@everywhere testMeshX = $testMesh; 
	@everywhere thermoX   = $thermo;
	@everywhere cellsThreadsX = $cellsThreads;
	@everywhere nodesThreadsX = $nodesThreads;
	@everywhere testfields2dX  = $testfields2d;
	@everywhere viscfields2dX  = $viscfields2d;
		
	@everywhere dynControlsX = $dynControls;
	@everywhere solControlsX = $solControls;
	@everywhere pControlsX = $pControls;
	@everywhere outputX = $output;

	timeVector = [];
	residualsVector1 = []; 
	residualsVector2 = []; 
	residualsVector3 = []; 
	residualsVector4 = []; 
	residualsVectorMax = ones(Float64,4);
	convergenceCriteria= [1e-5;1e-5;1e-5;1e-5;];
	
	
	debugSaveInit = false;
	if (debugSaveInit)
	
		rhoNodes = zeros(Float64,testMesh.nNodes);
		uxNodes = zeros(Float64,testMesh.nNodes);
		uyNodes = zeros(Float64,testMesh.nNodes);
		pNodes = zeros(Float64,testMesh.nNodes);
	
		cells2nodesSolutionReconstructionWithStencilsImplicitSA(nodesThreadsX, testMeshDistrX, testfields2dX, dummy); 
	
		for i = 1:testMesh.nNodes
			rhoNodes[i] = testfields2dX.densityNodes[i];
			uxNodes[i] = testfields2dX.UxNodes[i];
			uyNodes[i] = testfields2dX.UyNodes[i];
			pNodes[i] = testfields2dX.pressureNodes[i];
		end
		
		outputfileZero = string(outputfile,"_t=0");
		println("Saving  solution to  ", outputfileZero);
			#saveResults2VTK(outputfile, testMesh, densityF);
			saveResults4VTK(outputfileZero, testMesh, rhoNodes, uxNodes, uyNodes, pNodes);
		println("done ...  ");	
		
		
		@save outputfileZero solInst
		
	end
	
	
	
	maxEdge,id = findmax(testMeshX.HX);
	

	dt::Float64 =  solControls.dt;  
	@everywhere dtX = $dt; 
	@everywhere maxEdgeX = $maxEdge; 

	debug = true;	
	useArtViscoistyDapming = true;

	
	println("Start calculations ...");
	println(output.header);
	
	##if (!coldrun)
	
	
		#for l = 1:2
		while (dynControlsX.isRunSimulation == 1)
		
			
			##CPUtic();	
			start = time();
			
			
			# PROPAGATE STAGE: 
			(dynControlsX.velmax,id) = findmax(testfields2dX.VMAXCells);
			# #dynControls.tau = solControls.CFL * testMesh.maxEdgeLength/(max(dynControls.velmax,1.0e-6)); !!!!
			dynControlsX.tau = solControlsX.CFL * maxEdge/(max(dynControlsX.velmax,1.0e-6));
		
			
			if (useArtViscoistyDapming)
			
				calcArtificialViscositySA( cellsThreadsX, testMeshDistrX, testfields2dX, viscfields2dX);
				
				#cells2nodesSolutionReconstructionWithStencilsSerialUCons(testMeshX, UconsCellsOldX,  UconsNodesOldX );			
		
				calcDiffTerm(cellsThreadsX, testMeshDistrX, testfields2dX, viscfields2dX, thermoX, UconsNodesOldX, UConsDiffCellsX, UConsDiffNodesX);
			
			end
	
				
			## Explicit Euler first-order	
			calcOneStage(1.0, dtX, dynControlsX.flowTime, testMeshDistrX , testfields2dX, thermoX, cellsThreadsX,  UconsCellsOldX, iFLUXX, UConsDiffCellsX,  UconsCellsNewX);
			
			#doExplicitRK3TVD(1.0, dtX, testMeshDistrX , testfields2dX, thermoX, cellsThreadsX,  UconsCellsOldX, iFLUXX,  UConsDiffCellsX, 
			#  UconsCellsNew1X,UconsCellsNew2X,UconsCellsNew3X,UconsCellsNewX);
			
						
			
			@sync @distributed for p in workers()	
	
				beginCell::Int32 = cellsThreadsX[p-1,1];
				endCell::Int32 = cellsThreadsX[p-1,2];
				#println("worker: ",p,"\tbegin cell: ",beginCell,"\tend cell: ", endCell);
														
				updateVariablesSA(beginCell, endCell, thermoX.Gamma,  UconsCellsNewX, UconsCellsOldX, DeltaX, testfields2dX);
		
			end
			
			@everywhere finalize(updateVariablesSA);	
			
			
			# @sync @distributed for p in workers()	
	
				# beginNode::Int32 = nodesThreadsX[p-1,1];
				# endNode::Int32 = nodesThreadsX[p-1,2];
				
														
				# cells2nodesSolutionReconstructionWithStencilsDistributed(beginNode, endNode, 
				# testMeshDistrX, testfields2dX, viscfields2dX, UconsCellsOldX,  UconsNodesOldX);
		
			# end
			
			# @everywhere finalize(cells2nodesSolutionReconstructionWithStencilsDistributed);	
			
			cells2nodesSolutionReconstructionWithStencilsSerial(testMeshX,testfields2dX, viscfields2dX, UconsCellsOldX,  UconsNodesOldX);
								
			
			
	
			(dynControlsX.rhoMax,id) = findmax(testfields2dX.densityCells);
			(dynControlsX.rhoMin,id) = findmin(testfields2dX.densityCells);
			

			push!(timeVector, dynControlsX.flowTime); 
			dynControlsX.curIter += 1; 
			dynControlsX.verIter += 1;
				
			
			
			
			updateResidualSA(DeltaX, 
				residualsVector1,residualsVector2,residualsVector3,residualsVector4, residualsVectorMax,  
				convergenceCriteria, dynControlsX);
			
			
			updateOutputSA(timeVector,residualsVector1,residualsVector2,residualsVector3,residualsVector4, residualsVectorMax, 
				testMeshX, testfields2dX, viscfields2dX,  solControlsX, outputX, dynControlsX, solInst);
	
			
			# EVALUATE STAGE:
			
			dynControls.flowTime += dt; 
			##flowTimeX += dt;
			
			# if (solControlsX.timeStepMethod == 1)
				# dynControlsX.flowTime += dynControlsX.tau;  	
			# else
				# dynControlsX.flowTime += solControlsX.dt;  
			# end
			

	

			if (flowTime>= solControlsX.stopTime || dynControlsX.isSolutionConverged == 1)
				dynControlsX.isRunSimulation = 0;
		
				if (dynControlsX.isSolutionConverged == true)
					println("Solution converged! ");
				else
					println("Simultaion flow time reached the set Time!");
				end
			
				if (outputX.saveResiduals == 1)
					#println("Saving Residuals ... ");
					#cd(dynControlsX.localTestPath);
					#saveResiduals(output.fileNameResiduals, timeVector, residualsVector1, residualsVector2, residualsVector3, residualsVector4);
					#cd(dynControlsX.globalPath);
				end
				if (outputX.saveResults == 1)
					#println("Saving Results ... ");
					#cd(dynControlsX.localTestPath);
					#saveSolution(output.fileNameResults, testMeshX.xNodes, testMeshX.yNodes, UphysNodes);
					#cd(dynControlsX.globalPath);
				end
			
				
			
			end

			#dynControlsX.cpuTime  += CPUtoq(); 
			elapsed = time() - start;
			dynControlsX.cpuTime  += elapsed ; 
			
			if (dynControlsX.flowTime >= solControls.stopTime)
				dynControlsX.isRunSimulation = 0;
			end
			
		end ## end while
		 
		 
		
		
		solInst.dt = dtX;
		solInst.flowTime = dynControls.flowTime;
		for i = 1 : solInst.nCells
			solInst.densityCells[i] = testfields2dX.densityCells[i];
			solInst.UxCells[i] = testfields2dX.UxCells[i];
			solInst.UyCells[i] = testfields2dX.UyCells[i];
			solInst.pressureCells[i] = testfields2dX.pressureCells[i];
		end
		

		rhoNodes = zeros(Float64,testMesh.nNodes);
		uxNodes = zeros(Float64,testMesh.nNodes);
		uyNodes = zeros(Float64,testMesh.nNodes);
		pNodes = zeros(Float64,testMesh.nNodes);
	
		cells2nodesSolutionReconstructionWithStencilsImplicitSA(nodesThreadsX, testMeshDistrX, testfields2dX, dummy); 
	
		for i = 1:testMesh.nNodes
			rhoNodes[i] = testfields2dX.densityNodes[i];
			uxNodes[i] = testfields2dX.UxNodes[i];
			uyNodes[i] = testfields2dX.UyNodes[i];
			pNodes[i] = testfields2dX.pressureNodes[i];
		end
				
		println("Saving  solution to  ", outputfile);
			saveResults4VTK(outputfile, testMesh, rhoNodes, uxNodes, uyNodes, pNodes);
			@save outputfile solInst
		println("done ...  ");	
		
		 
		 
		
	#end ## if debug
	
end




##@time godunov2dthreads("2dmixinglayerUp_delta3.bson", numThreads, "2dMixingLayer_delta3", false); 
##@profview godunov2dthreads("2dmixinglayerUp_delta2.bson", numThreads, "2dMixingLayer_delta2", false); 


godunov2dthreads("testStep2dBaseTriSmooth", numThreads, "testStep2dBaseTriSmooth", false); 




