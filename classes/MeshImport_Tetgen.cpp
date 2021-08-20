#include "MeshImport_Tetgen.hpp"
#include "Output.hpp"

using namespace MainApplication;

namespace GeDiM
{
	// ***************************************************************************
	MeshImport_Tetgen::MeshImport_Tetgen() : GenericMeshImportInterface()
	{
		inputMeshPointer = NULL;
		outputMeshPointer = NULL;
		tetgenOptions = "Qpqfezna";
	}
	MeshImport_Tetgen::~MeshImport_Tetgen()
	{
		if (inputMeshPointer != NULL)
		{
			delete[] inputMeshPointer->pointlist; inputMeshPointer->pointlist = NULL;
			delete[] inputMeshPointer->pointmarkerlist; inputMeshPointer->pointmarkerlist = NULL;

			delete[] inputMeshPointer->edgelist; inputMeshPointer->edgelist = NULL;
			delete[] inputMeshPointer->edgemarkerlist; inputMeshPointer->edgemarkerlist = NULL;

			for (int f = 0; f < inputMeshPointer->numberoffacets; f++)
			{
				tetgenio::facet* tetgenFace = &inputMeshPointer->facetlist[f];

				tetgenio::polygon* tetgenPolygon =  &tetgenFace->polygonlist[0];
				delete[] tetgenPolygon->vertexlist; tetgenPolygon->vertexlist = NULL;

				delete[] tetgenFace->polygonlist; tetgenFace->polygonlist = NULL;
			}

			delete[] inputMeshPointer->facetlist; inputMeshPointer->facetlist = NULL;
		}

		delete inputMeshPointer; inputMeshPointer = NULL;
		delete outputMeshPointer; outputMeshPointer = NULL;
	}
	// ***************************************************************************
	Output::ExitCodes MeshImport_Tetgen::CreateTetgenInput(const GenericDomain& domain)
	{
		delete inputMeshPointer; inputMeshPointer = NULL;
		inputMeshPointer = new tetgenio();

		const GenericDomain& domain3D = domain;

		const unsigned int& numberOfVertices = domain3D.TotalNumberVertices();
		const unsigned int& numberOfConstrainedPoints = constrainedPoints.rows();
		const unsigned int& numberOfEdges = domain3D.TotalNumberEdges();
		const unsigned int& numberOfFaces = domain3D.TotalNumberFaces();
		const unsigned int& numberOfConstrainedFacets = constrainedFacets.size();

		if (numberOfVertices == 0 || numberOfEdges == 0 || numberOfFaces == 0)
		{
			Output::PrintErrorMessage("Wrong initialization of the domain %d, no vertices or faces", false, domain3D.GlobalId());
			return Output::GenericError;
		}

		inputMeshPointer->firstnumber = 0;
		inputMeshPointer->numberofpoints = numberOfVertices + numberOfConstrainedPoints;
		inputMeshPointer->pointlist = new REAL[(numberOfVertices + numberOfConstrainedPoints) * 3];
		inputMeshPointer->pointmarkerlist = new int[numberOfVertices + numberOfConstrainedPoints];

		inputMeshPointer->numberofedges = numberOfEdges;
		inputMeshPointer->edgelist = new int[numberOfEdges * 2];
		inputMeshPointer->edgemarkerlist = new int[numberOfEdges];

		inputMeshPointer->numberoffacets = numberOfFaces + numberOfConstrainedFacets;
		inputMeshPointer->facetlist = new tetgenio::facet[numberOfFaces + numberOfConstrainedFacets];
		inputMeshPointer->facetmarkerlist = new int[numberOfFaces + numberOfConstrainedFacets];

		double* point_list = inputMeshPointer->pointlist;
		int* point_markerlist = inputMeshPointer->pointmarkerlist;

		int* edge_list = inputMeshPointer->edgelist;
		int* edge_markerlist = inputMeshPointer->edgemarkerlist;

		tetgenio::facet* face_list = inputMeshPointer->facetlist;
		int* face_markerlist = inputMeshPointer->facetmarkerlist;


		for (unsigned int v = 0; v < numberOfVertices; v++)
		{
			point_list[3 * v] = domain3D.Vertex(v)(0);
			point_list[3 * v + 1] = domain3D.Vertex(v)(1);
			point_list[3 * v + 2] = domain3D.Vertex(v)(2);

			point_markerlist[v] = 2;
		}

		for (unsigned int e = 0; e < numberOfEdges; e++)
		{
			edge_list[2 * e] = domain3D.EdgeOriginIndex(e);
			edge_list[2 * e + 1] = domain3D.EdgeEndIndex(e);

			edge_markerlist[e]= 2;
		}

		for (unsigned int f = 0; f < numberOfFaces; f++)
		{
			tetgenio::facet* tetgenFace = &face_list[f];

			tetgenFace->numberofpolygons = 1;
			tetgenFace->polygonlist = new tetgenio::polygon[1];
			tetgenFace->numberofholes = 0;
			tetgenFace->holelist = NULL;

			tetgenio::polygon* tetgenPolygon =  &tetgenFace->polygonlist[0];

			const size_t numberFacePoints = domain3D.NumberFacePoints(f);
			tetgenPolygon->numberofvertices = numberFacePoints;
			tetgenPolygon->vertexlist = new int[tetgenPolygon->numberofvertices];

			for (unsigned int v = 0; v < numberFacePoints; v++)
				tetgenPolygon->vertexlist[v] = domain3D.FacePointIndex(f, v);

			face_markerlist[f] = 2;
		}

		if (!unidimensionalVertexMarkers.empty())
			memcpy(point_markerlist, unidimensionalVertexMarkers.data(), numberOfVertices * sizeof(int));
		if (!unidimensionalEdgeMarkers.empty())
			memcpy(edge_markerlist, unidimensionalEdgeMarkers.data(), numberOfEdges * sizeof(int));
		if (!unidimensionalFaceMarkers.empty())
		{
			memcpy(face_markerlist, unidimensionalFaceMarkers.data(), numberOfFaces * sizeof(int));
			CreateMappedMarkers();
		}

		if(numberOfConstrainedPoints > 0)
		{
			MatrixXd pointsCoordinate = constrainedPoints.leftCols(3);

			for (unsigned int j = 0; j < numberOfConstrainedPoints; j++)
			{
				point_list[3 * (numberOfVertices + j)] = pointsCoordinate(j, 0);
				point_list[3 * (numberOfVertices + j) + 1] = pointsCoordinate(j, 1);
				point_list[3 * (numberOfVertices + j) + 2] = pointsCoordinate(j, 2);

				point_markerlist[(numberOfVertices + j)] = 2;
			}
		}

		for(unsigned int numFac = 0; numFac < numberOfConstrainedFacets; numFac++)
		{
			tetgenio::facet* tetgenFace = &face_list[numFac + numberOfFaces];

			tetgenFace->numberofpolygons = 1;
			tetgenFace->polygonlist = new tetgenio::polygon[1];
			tetgenFace->numberofholes = 0;
			tetgenFace->holelist = NULL;

			tetgenio::polygon* tetgenPolygon =  &tetgenFace->polygonlist[0];

			const vector<int>& idPoints = constrainedFacets[numFac];
			const size_t numberFacePoints = idPoints.size();
			tetgenPolygon->numberofvertices = numberFacePoints;
			tetgenPolygon->vertexlist = new int[tetgenPolygon->numberofvertices];

			for (unsigned int v = 0; v < numberFacePoints; v++)
                tetgenPolygon->vertexlist[v] = idPoints[v] + numberOfVertices;

			face_markerlist[numFac + numberOfFaces] = 2;
		}


		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes MeshImport_Tetgen::CreateTetgenOutput(const GenericDomain& domain)
	{
		if (minimumNumberOfCells == 0 && maximumCellSize <= 0)
		{
			Output::PrintErrorMessage("Wrong initialization of the minimumNumberOfCells or minimumCellSize", false);
			return Output::GenericError;
		}

		if (inputMeshPointer == NULL)
		{
			Output::PrintErrorMessage("No Tetgen input in domain %d", false, domain.GlobalId());
			return Output::GenericError;
		}

		const GenericDomain& domain3D = domain;

		if (minimumNumberOfCells > 0 && domain3D.Measure() == 0)
		{
			Output::PrintErrorMessage("Wrong initialization of the domain %d", false, domain3D.GlobalId());
			return Output::GenericError;
		}

		const double& domainVolume = domain3D.Measure();
		double cellVolume = minimumNumberOfCells == 0 ? maximumCellSize : domainVolume / (double)minimumNumberOfCells;

		tetgenbehavior b;

		ostringstream options;
		options.precision(16);
		options<< tetgenOptions;
		options<< cellVolume;
		size_t sizeOptions = options.str().size();
		char* optionPointer = new char[sizeOptions + 1];
		options.str().copy(optionPointer, sizeOptions);
		optionPointer[sizeOptions] = '\0';

		b.parse_commandline(optionPointer);

		delete outputMeshPointer; outputMeshPointer = NULL;
		outputMeshPointer = new tetgenio();

		tetrahedralize(&b, inputMeshPointer, outputMeshPointer);

		delete[] optionPointer;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes MeshImport_Tetgen::CreateMesh(const GenericDomain& domain, GenericMesh& mesh) const
	{
		/// <ul>

		if (outputMeshPointer == NULL)
		{
			Output::PrintErrorMessage("No Tetgen ouput in domain %d", false, domain.GlobalId());
			return Output::GenericError;
		}

		const GenericDomain& domain3D = domain;
		const tetgenio& tetgenMesh = *outputMeshPointer;

		/// <li>	Fill mesh structures
		unsigned int numberOfCellsMesh = tetgenMesh.numberoftetrahedra;
		unsigned int numberOfFacesMesh = tetgenMesh.numberoftrifaces;
		unsigned int numberOfEgdesMesh = tetgenMesh.numberofedges;
		unsigned int numberOfPointsMesh = tetgenMesh.numberofpoints;

		mesh.InitializeCells(numberOfCellsMesh);
		mesh.InitializeFaces(numberOfFacesMesh);
		mesh.InitializePoints(numberOfPointsMesh);
		mesh.InitializeEdges(numberOfEgdesMesh);
		vector<GenericCell*> cells(numberOfCellsMesh);
		vector<GenericFace*> faces(numberOfFacesMesh);
		vector<GenericEdge*> edges(numberOfEgdesMesh);
		vector<GenericPoint*> points(numberOfPointsMesh);

		/// <li> Set Points
		for (unsigned int p = 0; p < numberOfPointsMesh; p++)
		{
			points[p] = mesh.CreatePoint();

			GenericPoint* point = points[p];

			Vector3d point3d(tetgenMesh.pointlist[3 * p], tetgenMesh.pointlist[3 * p + 1], tetgenMesh.pointlist[3 * p + 2]);

			point->SetCoordinates(point3d);
			if(!unidimensionalVertexMarkers.empty())
			{
				int position = tetgenMesh.pointmarkerlist[p];
				(position > 0)? point->SetMarkers(mappedMarkers[position - 1]) : point->SetMarkers(vector<unsigned int>(markerDimension, 0));
			}
			else
			{
				point->SetMarkers(vector<unsigned int>(markerDimension, 2));
			}
			point->InitializeCells(4);
			point->InitializeFaces(4);
			point->InitializeEdges(4);


			mesh.AddPoint(point);
		}

		/// <li> Set Edges
		map< vector<unsigned int>, unsigned int > endpointsToEdge;
		vector<unsigned int> endpoints(2);

		for(unsigned int ed = 0; ed < numberOfEgdesMesh; ed++)
		{
			edges[ed] = mesh.CreateEdge();

			GenericEdge* edge = edges[ed];

			if(!unidimensionalEdgeMarkers.empty())
			{
				int position = tetgenMesh.edgemarkerlist[ed];
				(position > 0)? edge->SetMarkers(mappedMarkers[position - 1]) : edge->SetMarkers(vector<unsigned int>(markerDimension, 0));
			}
			else
			{
				edge->SetMarkers(vector<unsigned int>(markerDimension, 2));
			}
			edge->InitializeCells(10);
			edge->InitializeFaces(10);

			for(int i = 0; i < 2; i++)
			{
				GenericPoint* point = points[tetgenMesh.edgelist[2 * ed + i]];
				edge->AddPoint(point);
				point->AddEdge(edge);
			}

			endpoints[0] = tetgenMesh.edgelist[2 * ed];
			endpoints[1] = tetgenMesh.edgelist[2 * ed + 1];
			sort(endpoints.begin(), endpoints.end());
			endpointsToEdge.insert(std::pair< vector<unsigned int>, unsigned int>(endpoints, edge->Id()));

			mesh.AddEdge(edge);
		}

		/// <li> Set Faces
		map< vector<unsigned int> , unsigned int > verticesToFace;
		vector<unsigned int> vertices(3);
		vector<unsigned int> edgeEndPoints(2);

		for(unsigned int f = 0; f < numberOfFacesMesh; f++)
		{
			faces[f] = mesh.CreateFace();

			GenericFace* face = faces[f];

			face->AllocateCells(2);
			face->InitializePoints(3);
			face->InitializeEdges(3);

			if(!unidimensionalFaceMarkers.empty())
			{
				int position = tetgenMesh.trifacemarkerlist[f];
				(position > 0)? face->SetMarkers(mappedMarkers[position - 1]) : face->SetMarkers(vector<unsigned int>(markerDimension, 0));
			}
			else
			{
				for(unsigned int numMark = 0; numMark < markerDimension; numMark++)
					face->SetMarker(2, numMark);

				face->SetMarkers(vector<unsigned int>(markerDimension, 2));
			}

			for (unsigned int j = 0; j < 3; j++)
			{
				GenericPoint* point = points[tetgenMesh.trifacelist[3 * f + j]];
				face->AddPoint(point);
				point->AddFace(face);
			}

			for (unsigned int j = 0; j < 3; j++)
				vertices[j] = tetgenMesh.trifacelist[3 * f + j];

			sort(vertices.begin(),vertices.end());
			verticesToFace.insert(std::pair< vector<unsigned int>, unsigned int >(vertices, face->Id()));

			for (unsigned int j = 0; j < 3; j++)
			{
				edgeEndPoints[0] = face->Point(j)->Id();
				edgeEndPoints[1] = face->Point((j + 1) % 3)->Id();

				sort(edgeEndPoints.begin(),edgeEndPoints.end());

				map< vector<unsigned int>, unsigned int>::iterator edgeFound = endpointsToEdge.find(edgeEndPoints);
				if (edgeFound == endpointsToEdge.end())
				{
					Output::PrintErrorMessage("Error retrieving an edge when adding edges to face %d in domain %d", true, face->GlobalId(), domain3D.GlobalId());
					return Output::GenericError;
				}

				face->AddEdge(edges[edgeFound->second]);
				edges[edgeFound->second]->AddFace(face);
			}

			mesh.AddFace(face);
		}

		/// <li> Set Cells
		vector<unsigned int> faceVertices(3);
		for (unsigned int c = 0; c < numberOfCellsMesh; c++)
		{
			cells[c] = mesh.CreateCell();

			GenericCell* cell = cells[c];

			cell->InitializePoints(4);
			cell->AllocateCells(4);
			cell->InitializeFaces(4);
			cell->InitializeEdges(6);

			for (int i = 0; i < 4; i++)
			{
				GenericPoint* point = points[tetgenMesh.tetrahedronlist[tetgenMesh.numberofcorners * c + i]];

				cell->AddPoint(point);
				point->AddCell(cell);
			}

			for (unsigned int j = 0; j < 4; j++)
			{
				for(unsigned int k = 0; k < 3; k++)
					faceVertices[k] = cell->Point((j + k) % 4)->Id();

				sort(faceVertices.begin(),faceVertices.end());
				map< vector<unsigned int>, unsigned int>::iterator verticesToFaceFound = verticesToFace.find(faceVertices);

				if (verticesToFaceFound == verticesToFace.end())
				{
					Output::PrintErrorMessage("Error retrieving faces when adding faces to cell %d in domain %d", true, cell->GlobalId(), domain3D.GlobalId());
					return Output::GenericError;
				}

				unsigned int faceId = verticesToFaceFound->second;
				cell->AddFace(faces[faceId]);

				if(faces[faceId]->Cell(0) == NULL)
					faces[faceId]->InsertCell(cell, 0);
				else
					faces[faceId]->InsertCell(cell, 1);

				if (j < 3)
				{
					for(unsigned int k = 0; k < 3; k++)
					{
						bool flag = true;
						unsigned int pos = 0;

						while (pos < cell->NumberOfEdges() && flag)
							flag = cell->Face(j)->Edge(k) != cell->Edge(pos++);

						if(flag)
						{
							cell->AddEdge(cell->Face(j)->Edge(k));
							edges[cell->Face(j)->Edge(k)->Id()]->AddCell(cell);
						}
					}
				}
			}

			mesh.AddCell(cell);
		}
		if ( outputMeshPointer->neighborlist != NULL)
		{
			for (unsigned int c = 0; c < numberOfCellsMesh; c++)
			{
				GenericCell* cell = cells[c];

				for (int i = 0; i < 4; i++)
				{
					if (outputMeshPointer->neighborlist[tetgenMesh.numberofcorners * c + i] > -1)
					{
						GenericCell* cellNeigh = cells[tetgenMesh.neighborlist[4 * c + i]];
						cell->InsertCell(cellNeigh, (i+1)%4);
					}
				}
			}
		}

		return Output::Success;

		/// </ul>
	}
	// ***************************************************************************
	Output::ExitCodes MeshImport_Tetgen::ExportTetgenMesh(const string& nameFolder, const string& nameFile) const
	{
		if (outputMeshPointer == NULL)
			return Output::GenericError;

		ostringstream nameFolderStream, nameFileStream;

		nameFolderStream<< nameFolder<< "/";
		nameFolderStream<< "Tetgen/";

		Output::CreateFolder(nameFolderStream.str());

		nameFileStream<< nameFolderStream.str()<< nameFile;

		Output::CreateFolder(nameFolderStream.str());

		outputMeshPointer->firstnumber = 0;
		outputMeshPointer->save_nodes((char*)nameFileStream.str().c_str());
		outputMeshPointer->save_elements((char*)nameFileStream.str().c_str());
		outputMeshPointer->save_faces((char*)nameFileStream.str().c_str());
		outputMeshPointer->save_edges((char*)nameFileStream.str().c_str());

		return Output::Success;
	}
	// ***************************************************************************
	/*Output::ExitCodes TetgenVemInterface::TetgenToVemMesh(const tetgenio& tetgenMesh, VemMesh& mesh) const
		{
		/// Add Points
		mesh.InitializePoints(tetgenMesh.numberofpoints);
		for(unsigned int i = 0; i < tetgenMesh.numberofpoints; i++)
		{
		VemPoint& point = *(mesh.CreatePoint());
		point.SetCoordinates(tetgenMesh.pointlist[3*i],tetgenMesh.pointlist[3*i+1],tetgenMesh.pointlist[3*i+2]);
		point.SetMarker(tetgenMesh.pointmarkerlist[i]);
		mesh.AddPoint(&point);
		}

		/// Add Edges
		mesh.InitializeEdges(tetgenMesh.numberofedges);
		map< vector<unsigned int> , unsigned int > endpointsToEdge;
		vector<unsigned int> endpoints (2);
		for(unsigned int i = 0; i < tetgenMesh.numberofedges; i++)
		{
		VemEdge& edge = *(mesh.CreateEdge());
		endpoints[0] = tetgenMesh.edgelist[2*i];
		endpoints[1] = tetgenMesh.edgelist[2*i+1];
		edge.AddPoint(mesh.Point(endpoints[0]));
		edge.AddPoint(mesh.Point(endpoints[1]));
		sort(endpoints.begin(),endpoints.end());
		endpointsToEdge.insert(std::pair< vector<unsigned int>, unsigned int>(endpoints, edge.Id()));
		edge.SetMarker(tetgenMesh.edgemarkerlist[i]);
		mesh.AddEdge(&edge);
		}

		/// Add Faces
		map< vector<unsigned int> , unsigned int > verticesToFace;
		vector<unsigned int> vertices (3);
		mesh.InitializeFaces(tetgenMesh.numberoftrifaces);
		for(unsigned int i = 0; i < tetgenMesh.numberoftrifaces; i++)
		{
		VemFace& face = *(mesh.CreateFace());
		face.InitializePoints(3);
		/// Points
		for(unsigned int j = 0; j < 3; j++)
		{
		vertices[j] = tetgenMesh.trifacelist[3*i+j];
		face.AddPoint(mesh.Point(vertices[j]));
		}
		sort(vertices.begin(),vertices.end());
		verticesToFace.insert(std::pair< vector<unsigned int>, unsigned int >(vertices, face.Id()));
		/// Edges
		vector<unsigned int> edgeEndPoints (2);
		face.InitializeEdges(3);
		for(unsigned int j = 0; j < 3; j++)
		{
		edgeEndPoints[0] = face.Point(j)->Id();
		edgeEndPoints[1] = face.Point((j+1)%3)->Id();
		sort(edgeEndPoints.begin(),edgeEndPoints.end());
		if(endpointsToEdge.find(edgeEndPoints)==endpointsToEdge.end())
		{
		Output::PrintErrorMessage("%s: error retrieving an edge when adding edges to faces", true, __func__);
		exit(-1);
		}
		face.AddEdge(mesh.Edge(endpointsToEdge.find(edgeEndPoints)->second));
		}
		face.SetMarker(tetgenMesh.trifacemarkerlist[i]);
		face.InitializeCells(2);
		mesh.AddFace(&face);
		}

		/// Add Cells
		mesh.InitializeCells(tetgenMesh.numberoftetrahedra);
		map<vector<unsigned int> , unsigned int>::iterator verticesToFace_iterator;
		for(unsigned int i = 0; i < tetgenMesh.numberoftetrahedra; i++)
		{
		VemCell& cell = *(mesh.CreateCell3D());
		/// Points
		cell.InitializePoints(4);
		bool flag = false;
		cell.AddPoint(mesh.Point(tetgenMesh.tetrahedronlist[tetgenMesh.numberofcorners*i]));
		cell.AddPoint(mesh.Point(tetgenMesh.tetrahedronlist[tetgenMesh.numberofcorners*i+1]));
		cell.AddPoint(mesh.Point(tetgenMesh.tetrahedronlist[tetgenMesh.numberofcorners*i+2]));
		cell.AddPoint(mesh.Point(tetgenMesh.tetrahedronlist[tetgenMesh.numberofcorners*i+3]));
		/// Faces and Edges
		cell.InitializeFaces(4);
		cell.InitializeEdges(6);
		vector<unsigned int> faceVertices (3);
		for(unsigned int j = 0; j < 4; j++)
		{
		for(unsigned int k = 0; k < 3; k++)
		faceVertices[k] = cell.Point((j+k)%4)->Id();
		sort(faceVertices.begin(),faceVertices.end());
		verticesToFace_iterator = verticesToFace.find(faceVertices);
		if(verticesToFace_iterator != verticesToFace.end())
		{
		unsigned int faceId = verticesToFace_iterator->second;
		cell.AddFace(mesh.Face(faceId));
		mesh.AddFaceNeighbourCell( faceId, &cell );
		}
		else
		{
		Output::PrintErrorMessage("%s: error retrieving faces when adding faces to cells", true, __func__);
		exit(-1);
		}
		if(j < 3)
		{
		for(unsigned int k = 0; k < 3; k++)
		{
		bool flag = true;
		unsigned int pos = 0;
		while(pos < cell.NumberOfEdges() && flag)
		flag = cell.Face(j)->Edge(k) != cell.Edge(pos++);
		if(flag)
		cell.AddEdge(cell.Face(j)->Edge(k));
		}
		}
		}
		if(cell.NumberOfEdges() != 6)
		{
		Output::PrintErrorMessage("%s: error retrieving edges when adding edges to cells", true, __func__);
		exit(-1);
		}
		mesh.AddCell(&cell);
		}
		return Output::Success;
		}*/
	// ***************************************************************************
}
