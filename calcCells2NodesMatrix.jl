

function computeCells2Nodes2D(testMesh::mesh2d, cells2nodes::SharedArray{Int64,2} )


  #int n = get_num_cells();
  #cells2nodes.SetDim(n,8);
  #cells2nodes.Fill();

  #cout << "Computing cells2nodes matrix ... ";

  for C=1:testMesh.nCells
  

     
      num_nodes::Int64 = testMesh.mesh_connectivity[C,3];
      nodesC1::Int64 = testMesh.mesh_connectivity[C,4];
      nodesC2::Int64 = testMesh.mesh_connectivity[C,5];
      nodesC3::Int64 = testMesh.mesh_connectivity[C,6];
      nodesC4::Int64 = testMesh.mesh_connectivity[C,7];


      if (num_nodes == 3) ## triangle
      
        for T = 1:num_nodes
          

              neib_cell::Int64 = testMesh.cell_stiffness[C,T];
              node1::Int64 = 0;
              node2::Int64 = 0;

			  
			  nodesP = [];

              if (neib_cell >0) ## internal cell
              
                  nodesT1::Int64 = testMesh.mesh_connectivity[neib_cell,4];
                  nodesT2::Int64 = testMesh.mesh_connectivity[neib_cell,5];
                  nodesT3::Int64 = testMesh.mesh_connectivity[neib_cell,6];

                  if (nodesC1 == nodesT1 ||  nodesC1 == nodesT2 || nodesC1 == nodesT3)
				    push!(nodesP, nodesC1);
                    ##nodesP.push_back(nodesC1);
				  end
                  if (nodesC2 == nodesT1 ||  nodesC2 == nodesT2 || nodesC2 == nodesT3)
					push!(nodesP, nodesC2);
                    ##nodesP.push_back(nodesC2);
				  end
                  if (nodesC3 == nodesT1 ||  nodesC3 == nodesT2 || nodesC3 == nodesT3)
					push!(nodesP, nodesC3);
                    ##nodesP.push_back(nodesC3);
				  end

                  if (size(nodesP,1) == 2)
                  
                    node1 = nodesP[1];
                    node2 = nodesP[2];
                  
                  else
                  
                      display("something wrong in creating cells2node matrix ... ");
                      ##throw(-1);
                  end

              
              else ## boundary cell
              

                  if (T == 1)
                  
                      node1 = nodesC1;
                      node2 = nodesC2;
                  
                  elseif (T == 2)
                  
                       node1 = nodesC2;
                       node2 = nodesC3;
                  
                  elseif (T == 3)
                  
                       node1 = nodesC3;
                       node2 = nodesC1;
                  end


              end ##  end if 


              if (T == 1)
              
                  cells2nodes[C,1] = node1;
                  cells2nodes[C,2] = node2;
              
              elseif (T == 2)
              
                  cells2nodes[C,3] = node1;
                  cells2nodes[C,4] = node2;
              
              elseif (T == 3)
              
                  cells2nodes[C,5] = node1;
                  cells2nodes[C,6] = node2;
              
              else
              
                  display("something wrong in creating cells2node matrix ... ");
                  ##throw(-1);
              end




		end ## end for particular triangular cell

		## end for triangle cells
	  
	  
      elseif (num_nodes == 4) ## quad element
      

          for T = 1:num_nodes
          

              neib_cell::Int64 = testMesh.cell_stiffness[C,T];
              node1::Int64 = 0;
              node2::Int64 = 0;


              nodesP = [];


              if (neib_cell >0) ## internal cell
              
                  nodesT1::Int64 = testMesh.mesh_connectivity[neib_cell,4];
                  nodesT2::Int64 = testMesh.mesh_connectivity[neib_cell,5];
                  nodesT3::Int64 = testMesh.mesh_connectivity[neib_cell,6];
                  nodesT4::Int64 = testMesh.mesh_connectivity[neib_cell,7];

                  if (nodesC1 == nodesT1 ||  nodesC1 == nodesT2 || nodesC1 == nodesT3 || nodesC1 == nodesT4)
				    push!(nodesP,nodesC1);
                    #nodesP.push_back(nodesC1);
				  end
                  if (nodesC2 == nodesT1 ||  nodesC2 == nodesT2 || nodesC2 == nodesT3 || nodesC2 == nodesT4)
				    push!(nodesP,nodesC2);
                    ##nodesP.push_back(nodesC2);
				  end
                  if (nodesC3 == nodesT1 ||  nodesC3 == nodesT2 || nodesC3 == nodesT3 || nodesC3 == nodesT4)
				    push!(nodesP,nodesC3); 
                    ##nodesP.push_back(nodesC3);
				  end
                  if (nodesC4 == nodesT1 ||  nodesC4 == nodesT2 || nodesC4 == nodesT3 || nodesC4 == nodesT4)
					push!(nodesP,nodesC4);
                    ##nodesP.push_back(nodesC4);
				  end


                  if (size(nodesP,1)  == 2)
                  
                    node1 = nodesP[1];
                    node2 = nodesP[2];
                  
                  else
                  
                      display("something wrong in creating cells2node matrix ... ");
                      ##throw(-1);
                  end

              
              else

                  if (T == 1)
                  
                      node1 = nodesC1;
                      node2 = nodesC2;
                  
                  elseif (T == 2)
                  
                       node1 = nodesC2;
                       node2 = nodesC3;
                  
                  elseif (T == 3)
                  
                       node1 = nodesC3;
                       node2 = nodesC4;
                  
                  elseif (T == 4)
                  
                       node1 = nodesC4;
                       node2 = nodesC1;
                  end


              end  ##  end boundary cell




              if (T == 1)
              
                  cells2nodes[C,1] = node1;
                  cells2nodes[C,2] = node2;
              
              elseif (T == 2)
              
                  cells2nodes[C,3] = node1;
                  cells2nodes[C,4] = node2;
            
              elseif (T == 3)
			  
                  cells2nodes[C,5] = node1;
                  cells2nodes[C,6] = node2;

              elseif (T == 4)
              
                  cells2nodes[C,7] = node1;
                  cells2nodes[C,8] = node2;

              else
              
                  display("something wrong in creating cells2node matrix ... ");
                  throw(-1);
              end

            end ## 
			
			## end for quad element

      else
      
          display("something wrong in creating cells2node matrix ... ");
          display("unknown element type: must be triangle or quad ... ");
		  
          ##throw(-1);
      end




  end ## end global loop for cells

  ## cout << "done " << endl;

end ## function
