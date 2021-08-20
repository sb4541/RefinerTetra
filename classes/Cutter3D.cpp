#include "Cutter3D.hpp"

#include "MacroDefinitions.hpp"
#include "Intersector2D1D.hpp"
#include "Intersector2D2D.hpp"
#include "algorithm"
#include <set>

namespace GeDiM
{
    // ***************************************************************************
	// Cutter2D Implementation
	// ***************************************************************************


	//****************************************************************************************************
	bool CutterByDomains::CompareByPlane(const GenericDomain* first, const GenericDomain* second)
	{
		//		Vector3d tangentVectorEdge = first->Vertex(1) - first->Vertex(0);
		Vector3d tangentVectorEdgeSecond = second->Vertex(1) - second->Vertex(0);

		//		Vector3d tangentVectorDiffFirst = second->Vertex(0) - first->Vertex(0);
		//		Vector3d tangentVectorDiffSecond = second->Vertex(1) - first->Vertex(0);
		Vector3d tangentVectorDiffThird = - second->Vertex(0) + first->Vertex(0);
		Vector3d tangentVectorDiffFourth = - second->Vertex(0) + first->Vertex(1);


		double toll = 1.0E-5;

		//true at the right of the tangentVector false at the left
		//		bool positionThirdPoint = (tangentVectorEdge.x() * tangentVectorDiffFirst.y() - tangentVectorDiffFirst.x() * tangentVectorEdge.y()) > toll;
		//		bool positionFourthPoint = (tangentVectorEdge.x() * tangentVectorDiffSecond.y() - tangentVectorDiffSecond.x() * tangentVectorEdge.y()) > toll;

		bool positionFirstPoint = (tangentVectorEdgeSecond.x() * tangentVectorDiffThird.y() - tangentVectorDiffThird.x() * tangentVectorEdgeSecond.y()) > toll;
		bool positionSecondPoint = (tangentVectorEdgeSecond.x() * tangentVectorDiffFourth.y() - tangentVectorDiffFourth.x() * tangentVectorEdgeSecond.y()) > toll;

		if(positionFirstPoint == positionSecondPoint)
			return true;
		else
			return false;
	}
	//****************************************************************************************************
	const Output::ExitCodes CutterByDomains::Initialize(const int& numberOfDomain)
	{
		if(domains.size() > 0)
		{
			Output::PrintErrorMessage("Cutter is already initialized", true);
			return Output::GenericError;
		}
		domains.reserve(numberOfDomain);
		return Output::Success;
	}
	//****************************************************************************************************
	CutterByDomains::CutterByDomains()
	{
		mesh = NULL;
	}
	CutterByDomains::~CutterByDomains()
	{
		domains.clear();
	}

	const Output::ExitCodes CutterByDomains::AddDomain(GenericDomain& domain)
	{
		if(domains.size() + 1 > domains.capacity())
		{
			Output::PrintErrorMessage("All domains already added", false);
			return Output::GenericError;
		}
		positionDomains.insert(pair<unsigned int, unsigned int>(domain.GlobalId(), domains.size()));
		domains.push_back(&domain);

		return Output::Success;
	}

	// ***************************************************************************
	// Cutter3D Implementation
	// ***************************************************************************

	//*************************
	CutterMesh3D::CutterMesh3D()
	{

		intersector2D1D = new Intersector2D1D();
	}

	CutterMesh3D::~CutterMesh3D()
	{
		valuePoints.clear();
		idPositionPoints.clear();

		valueEdges.clear();
		idPositionEdges.clear();

		positivePoints.clear();
		negativePoints.clear();
		planePoints.clear();

		positiveEdges.clear();
		negativeEdges.clear();
		planeEdges.clear();

		positiveFaces.clear();
		negativeFaces.clear();
		planeFaces.clear();

		idFacesToUpdate.clear();
		idCellsToUpdate.clear();

		delete intersector2D1D;
	}


	//***********************************************************************************************************************
	const Output::ExitCodes CutterMesh3D::CutFace(GenericFace& face, const vector<const GenericPoint*>& domainPointsInFace, const double& toll)
	{
		ostringstream key;
		key << "TestDomain";
		face.SetState(false);
		face.InitializeChilds(2);

		GenericFace& newFace_1 = *(mesh->CreateFace());
		mesh->AddFace(&newFace_1);
		face.AddChild(&newFace_1);
		newFace_1.SetFather(&face);
		newFace_1.InheritPropertiesByFather();
		newFace_1.AllocateCells(2);
		newFace_1.InsertCell(face.Cell(0), 0);
		newFace_1.InsertCell(face.Cell(1), 1);
		positiveFaces.push_back(newFace_1.Id());

		//add face to point and to edge

		GenericFace& newFace_2 = *(mesh->CreateFace());
		mesh->AddFace(&newFace_2);
		face.AddChild(&newFace_2);
		newFace_2.SetFather(&face);
		newFace_2.InheritPropertiesByFather();
		newFace_2.AllocateCells(2);
		newFace_2.InsertCell(face.Cell(0), 0);
		newFace_2.InsertCell(face.Cell(1), 1);
		negativeFaces.push_back(newFace_2.Id());

		newFace_1.InitializeEdges(face.NumberOfEdges()+1);
		newFace_2.InitializeEdges(face.NumberOfEdges()+1);

		//add face to point and to edge
		for(unsigned int i = 0; i < face.NumberOfEdges(); i++)
		{
			if(face.Edge(i)->IsActive())
			{
				unsigned int idEdge = face.Edge(i)->Id();
				GenericEdge& edge = *mesh->Edge(idEdge);
				unsigned int& position = idPositionEdges.at(idEdge);

				if(valueEdges[position] == PositionInRelationToSurface::Positive)
				{
					for(unsigned int numFac = 0; numFac < edge.NumberOfFaces(); numFac++)
					{
						if(edge.Face(numFac)->Id() == face.Id())
							edge.InsertFace(&newFace_1, numFac);
					}
					newFace_1.AddEdge(&edge);
					//					edge.AddFace(&newFace_1);
					//					edge.ShrinkFaces();
				}
				else if(valueEdges[position] == PositionInRelationToSurface::Negative)
				{
					for(unsigned int numFac = 0; numFac < edge.NumberOfFaces(); numFac++)
					{
						if(edge.Face(numFac)->Id() == face.Id())
							edge.InsertFace(&newFace_2, numFac);
					}
					newFace_2.AddEdge(&edge);
					//					edge.AddFace(&newFace_2);
					//					edge.ShrinkFaces();
				}
				else
					return Output::GenericError;
			}
			else
			{
				unsigned int idEdge = face.Edge(i)->Child(0)->Id();
				unsigned int idEdgeOne = face.Edge(i)->Child(1)->Id();
				GenericEdge& edge = *mesh->Edge(idEdge);
				GenericEdge& edgeOne = *mesh->Edge(idEdgeOne);

				unsigned int& position = idPositionEdges.at(idEdge);
				if(valueEdges[position] == PositionInRelationToSurface::Positive)
				{
					newFace_1.AddEdge(&edge);
					for(unsigned int numFac = 0; numFac < edge.NumberOfFaces(); numFac++)
					{
						if(edge.Face(numFac)->Id() == face.Id())
							edge.InsertFace(&newFace_1, numFac);
					}

					newFace_2.AddEdge(&edgeOne);
					for(unsigned int numFac = 0; numFac < edgeOne.NumberOfFaces(); numFac++)
					{
						if(edgeOne.Face(numFac)->Id() == face.Id())
							edgeOne.InsertFace(&newFace_2, numFac);
					}
				}
				else if(valueEdges[position] == PositionInRelationToSurface::Negative)
				{
					newFace_1.AddEdge(&edgeOne);
					for(unsigned int numFac = 0; numFac < edgeOne.NumberOfFaces(); numFac++)
					{
						if(edgeOne.Face(numFac)->Id() == face.Id())
							edgeOne.InsertFace(&newFace_1, numFac);
					}

					newFace_2.AddEdge(&edge);
					for(unsigned int numFac = 0; numFac < edge.NumberOfFaces(); numFac++)
					{
						if(edge.Face(numFac)->Id() == face.Id())
							edge.InsertFace(&newFace_2, numFac);
					}
				}
				else
					return Output::GenericError;
			}
		}

		GenericEdge& newEdge = *(mesh->CreateEdge());

		mesh->AddEdge(&newEdge);
		newEdge.InitializeCells(10);
		newEdge.InitializeFaces(10);
		newEdge.InitializeProperty(key.str());
		newEdge.AddProperty(key.str(), new int(0) );

		planeEdges.insert(newEdge.Id());
		newEdge.AddPoint(domainPointsInFace[0]);
		newEdge.AddPoint(domainPointsInFace[1]);
		newEdge.AddFace(&newFace_1);
		newEdge.AddFace(&newFace_2);
		mesh->Point(domainPointsInFace[0]->Id())->AddEdge(&newEdge);
		mesh->Point(domainPointsInFace[1]->Id())->AddEdge(&newEdge);
		mesh->Point(domainPointsInFace[0]->Id())->ShrinkEdges();
		mesh->Point(domainPointsInFace[1]->Id())->ShrinkEdges();

		newFace_1.AddEdge(&newEdge);
		newFace_2.AddEdge(&newEdge);

		//		newFace_1.ShrinkEdges();
		//		newFace_2.ShrinkEdges();

		vector<bool> addedFacesToPoint(2, false);

		for(unsigned int i = 0; i < face.NumberOfPoints(); i++)
		{
			unsigned int idPoint = face.Point(i)->Id();
			map<unsigned int, unsigned int>::iterator it = idPositionPoints.find(idPoint);

			GenericPoint& point = *mesh->Point(idPoint);
			if(it != idPositionPoints.end())
			{
				if (valuePoints[it->second] == PositionInRelationToSurface::Positive)
					point.AddFace(&newFace_1);
				else if(valuePoints[it->second] == PositionInRelationToSurface::Negative)
					point.AddFace(&newFace_2);
				else if(valuePoints[it->second] == PositionInRelationToSurface::Plane)
				{
					point.AddFace(&newFace_1);
					point.AddFace(&newFace_2);
					if(idPoint == domainPointsInFace[0]->Id())
						addedFacesToPoint[0] = true;
					else
						addedFacesToPoint[1] = true;
				}
			}
			else
				return Output::GenericError;
			point.ShrinkFaces();
		}

		if(!addedFacesToPoint[0])
		{
			GenericPoint& point = *mesh->Point(domainPointsInFace[0]->Id());
			point.AddFace(&newFace_1);
			point.AddFace(&newFace_2);
			point.ShrinkFaces();
		}

		if(!addedFacesToPoint[1])
		{
			GenericPoint& point = *mesh->Point(domainPointsInFace[1]->Id());
			point.AddFace(&newFace_1);
			point.AddFace(&newFace_2);
			point.ShrinkFaces();
		}

		newFace_1.InitializePoints(newFace_1.NumberOfEdges());
		newFace_2.InitializePoints(newFace_2.NumberOfEdges());

		//Adding points to newFace_1

		newFace_1.AddPoint(newEdge.Point(0));
		newFace_1.AddPoint(newEdge.Point(1));
		unsigned int previousId = newEdge.Point(1)->Id();
		unsigned int j = 0;
		while(newFace_1.NumberOfPoints() < newFace_1.NumberOfEdges())
		{
			for(unsigned int i = j; i < newFace_1.NumberOfEdges() - 1; i++)
			{
				const GenericEdge* tmpEdge = newFace_1.Edge(i);
				unsigned int idPointZero = newFace_1.Edge(i)->Point(0)->Id();
				unsigned int idPointOne= newFace_1.Edge(i)->Point(1)->Id();

				if(idPointZero ==  previousId)
				{
					newFace_1.AddPoint(mesh->Point(idPointOne));
					previousId = idPointOne;

					newFace_1.InsertEdge(newFace_1.Edge(j), i);
					newFace_1.InsertEdge(tmpEdge, j);
					break;
				}
				else if(idPointOne ==  previousId)
				{
					newFace_1.AddPoint(mesh->Point(idPointZero));
					previousId = idPointZero;

					newFace_1.InsertEdge(newFace_1.Edge(j), i);
					newFace_1.InsertEdge(tmpEdge, j);
					break;
				}
			}
			j++;
		}

		unsigned int numEdgesNewFace = newFace_1.NumberOfEdges();
		const Vector3d& normalFace = face.Normal();
		newFace_1.ComputeNormal();
		Vector3d& normalNewFace = newFace_1.Normal();

		//the orientation is negative, therefore the vector of the points of the face is reversed
		if(normalFace.dot(normalNewFace) < - toll)
		{
			unsigned int forIteration =  ((numEdgesNewFace - 1) * 0.5);
			for(unsigned int i = 0; i < forIteration; i++)
			{
				const GenericPoint* tmpPoint = newFace_1.Point(i + 1);
				newFace_1.InsertPoint(newFace_1.Point(numEdgesNewFace - i - 1), i + 1);
				newFace_1.InsertPoint(tmpPoint, numEdgesNewFace - i - 1);
			}
			normalNewFace *= -1.0;
		}
		for(unsigned int k = 0; k < numEdgesNewFace; k++)
		{
			const GenericPoint* tmpPoint = newFace_1.Point(k);
			const GenericPoint* tmpPointNext = newFace_1.Point((k + 1)%numEdgesNewFace);
			const GenericEdge* tempEdge = newFace_1.Edge(k);
			const GenericPoint* tmpPointEd = tempEdge->Point(0);
			const GenericPoint* tmpPointEdNext = tempEdge->Point(1);
			if((tmpPoint == tmpPointEd && tmpPointNext == tmpPointEdNext) ||
				 (tmpPoint == tmpPointEdNext && tmpPointNext == tmpPointEd))
				continue;
			else
			{
				for(unsigned int ed = k; ed < numEdgesNewFace; ed++)
				{
					const GenericEdge* tempEdge2 = newFace_1.Edge(ed);
					const GenericPoint* tmpPointEd = tempEdge2->Point(0);
					const GenericPoint* tmpPointEdNext = tempEdge2->Point(1);
					if((tmpPoint == tmpPointEd && tmpPointNext == tmpPointEdNext) ||
						 (tmpPoint == tmpPointEdNext && tmpPointNext == tmpPointEd))
					{
						newFace_1.InsertEdge(tempEdge2, k);
						newFace_1.InsertEdge(tempEdge, ed);
					}
					else
						continue;
				}
			}
		}

		newFace_1.ComputeGeometricalProperties();

		//Adding points to newFace_2
		newFace_2.AddPoint(newEdge.Point(0));
		newFace_2.AddPoint(newEdge.Point(1));
		previousId = newEdge.Point(1)->Id();
		j = 0;
		while(newFace_2.NumberOfPoints() < newFace_2.NumberOfEdges())
		{
			for(unsigned int i = j; i < newFace_2.NumberOfEdges() - 1; i++)
			{
				const GenericEdge* tmpEdge = newFace_2.Edge(i);
				unsigned int idPointZero = newFace_2.Edge(i)->Point(0)->Id();
				unsigned int idPointOne= newFace_2.Edge(i)->Point(1)->Id();

				if(idPointZero ==  previousId)
				{
					newFace_2.AddPoint(mesh->Point(idPointOne));
					previousId = idPointOne;

					newFace_2.InsertEdge(newFace_2.Edge(j), i);
					newFace_2.InsertEdge(tmpEdge, j);
					break;

				}
				else if(idPointOne ==  previousId)
				{
					newFace_2.AddPoint(mesh->Point(idPointZero));
					previousId = idPointZero;

					newFace_2.InsertEdge(newFace_2.Edge(j), i);
					newFace_2.InsertEdge(tmpEdge, j);
					break;
				}
			}
			j++;
		}

		unsigned int numEdgesNewFace2 = newFace_2.NumberOfEdges();
		newFace_2.ComputeNormal();
		Vector3d& normalNewFace2 = newFace_2.Normal();
		//checking if the orientation of the points of the face is positive with regard to the normal to the domain

		//the orientation is negative, therfore the vector of the points of the face is reversed
		if(normalFace.dot(normalNewFace2) < - toll)
		{
			unsigned int forIteration = ((numEdgesNewFace2 - 1) * 0.5);
			for(unsigned int i = 0; i < forIteration; i++)
			{
				const GenericPoint* tmpPoint = newFace_2.Point(i + 1);
				newFace_2.InsertPoint(newFace_2.Point(numEdgesNewFace2 - i - 1), i + 1);
				newFace_2.InsertPoint(tmpPoint, numEdgesNewFace2 - i - 1);
			}
			normalNewFace2 *= -1.0;
		}
		for(unsigned int k = 0; k < numEdgesNewFace2; k++)
		{
			const GenericPoint* tmpPoint = newFace_2.Point(k);
			const GenericPoint* tmpPointNext = newFace_2.Point((k + 1)%numEdgesNewFace2);
			const GenericEdge* tempEdge = newFace_2.Edge(k);
			const GenericPoint* tmpPointEd = tempEdge->Point(0);
			const GenericPoint* tmpPointEdNext = tempEdge->Point(1);
			if((tmpPoint == tmpPointEd && tmpPointNext == tmpPointEdNext) ||
				 (tmpPoint == tmpPointEdNext && tmpPointNext == tmpPointEd))
				continue;
			else
			{
				for(unsigned int ed = k; ed < numEdgesNewFace2; ed++)
				{
					const GenericEdge* tempEdge2 = newFace_2.Edge(ed);
					const GenericPoint* tmpPointEd = tempEdge2->Point(0);
					const GenericPoint* tmpPointEdNext = tempEdge2->Point(1);
					if((tmpPoint == tmpPointEd && tmpPointNext == tmpPointEdNext) ||
						 (tmpPoint == tmpPointEdNext && tmpPointNext == tmpPointEd))
					{
						newFace_2.InsertEdge(tempEdge2, k);
						newFace_2.InsertEdge(tempEdge, ed);
					}
					else
						continue;
				}
			}
		}
		newFace_2.ComputeGeometricalProperties();
		return Output::Success;
	}
	//****************************************************************************************************
	const Output::ExitCodes CutterMesh3D::CutEdge(GenericEdge& edge, unsigned int& idNewPoint)
	{
		edge.SetState(false);
		edge.InitializeChilds(2);
		const GenericPoint& firstPoint = *edge.Point(0);
		const GenericPoint& secondPoint = *edge.Point(1);

		Vector3d tangent = secondPoint.Coordinates() -firstPoint.Coordinates();

		Intersector2D1D& intersector = *intersector2D1D;
		intersector.SetLine(firstPoint.Coordinates(), tangent);
		intersector.ComputeIntersection();
		Vector3d intersection = intersector.IntersectionPoint();

		GenericPoint& newPoint = *(mesh->CreatePoint());
		mesh->AddPoint(&newPoint);
		newPoint.SetCoordinates(intersection);
		idNewPoint = newPoint.Id();

		newPoint.InitializeFaces(edge.NumberOfFaces());
		for(unsigned int i = 0; i < edge.NumberOfFaces();i++)
			newPoint.AddFace(edge.Face(i));

		GenericEdge& newEdge_1 = *(mesh->CreateEdge());
		mesh->AddEdge(&newEdge_1);
		newEdge_1.AddPoint(edge.Point(0));
		newEdge_1.AddPoint(&newPoint);

		edge.AddChild(&newEdge_1);
		newEdge_1.SetFather(&edge);
		newEdge_1.InheritPropertiesByFather();

		GenericEdge& newEdge_2 = *(mesh->CreateEdge());
		mesh->AddEdge(& newEdge_2);
		newEdge_2.AddPoint(&newPoint);
		newEdge_2.AddPoint(edge.Point(1));

		edge.AddChild(&newEdge_2);
		newEdge_2.SetFather(&edge);
		newEdge_2.InheritPropertiesByFather();

		newPoint.InitializeEdges(2);

		if((positivePoints.find(firstPoint.Id()) != positivePoints.end() ))
		{
			newPoint.AddEdge(&newEdge_1);
			positiveEdges.insert(newEdge_1.Id());
			newPoint.AddEdge(&newEdge_2);
			negativeEdges.insert(newEdge_2.Id());
		}
		else
		{
			newPoint.AddEdge(&newEdge_2);
			positiveEdges.insert(newEdge_2.Id());
			newPoint.AddEdge(&newEdge_1);
			negativeEdges.insert(newEdge_1.Id());
		}

		mesh->Point(edge.Point(0)->Id())->AddEdge(&newEdge_1);
		mesh->Point(edge.Point(1)->Id())->AddEdge(&newEdge_2);

		newEdge_1.InitializeCells(10);
		newEdge_2.InitializeCells(10);
		for(unsigned int cel = 0; cel < edge.NumberOfCells(); cel++)
		{
			if( edge.Cell(cel) != NULL)
			{
				const GenericCell* cellTemp = edge.Cell(cel);
				newEdge_1.AddCell(cellTemp);
				newEdge_2.AddCell(cellTemp);
			}
		}

		newEdge_1.InitializeFaces(10);
		newEdge_2.InitializeFaces(10);
		for(unsigned int fac = 0; fac < edge.NumberOfFaces(); fac++)
		{
			if(edge.Face(fac) != NULL && edge.Face(fac)->IsActive())
			{
				const GenericFace* faceTemp = edge.Face(fac);
				idFacesToUpdate.insert(faceTemp->Id());
				newEdge_1.AddFace(faceTemp);
				newEdge_2.AddFace(faceTemp);
			}
		}

		return Output::Success;
	}

	//****************************************************************************************************
	const Output::ExitCodes CutterMesh3D::CutCell(GenericCell& cell, const Vector3d& normalPlane, const double& translation, const int& idDomain, const double& toll)
	{
		Reset();

		ostringstream key;
		key << "TestDomain";
		valuePoints.resize(cell.NumberOfPoints());
		for(unsigned int i = 0; i < cell.NumberOfPoints(); i++)
		{
			unsigned int positionOfPoint = PointDomainSplitter(*cell.Point(i), normalPlane, translation);

			idPositionPoints.insert(pair<unsigned int, unsigned int>(cell.Point(i)->Id(), i));
			switch(positionOfPoint)
			{
				case Positive:
				{
					positivePoints.insert(cell.Point(i)->Id());
					valuePoints[i] = PositionInRelationToSurface::Positive;
					break;
				}
				case Negative:
				{
					negativePoints.insert(cell.Point(i)->Id());
					valuePoints[i] = PositionInRelationToSurface::Negative;
					break;
				}
				case Plane:
				{
					planePoints.insert(cell.Point(i)->Id());
					valuePoints[i] = PositionInRelationToSurface::Plane;
					break;
				}
				default:
					return Output::GenericError;
			}
		}//pointsSplitted

		if(positivePoints.size() == 0 || negativePoints.size() == 0)
		{
			if(planePoints.size() > 2)
			{
				for(unsigned int numFac = 0; numFac < cell.NumberOfFaces(); numFac++)
				{
					const GenericFace& face = *cell.Face(numFac);
					if(planePoints.size() != face.NumberOfPoints())
						continue;
					set<unsigned int> pointFace;
					for(unsigned int numPnt = 0; numPnt < face.NumberOfPoints(); numPnt++)
						pointFace.insert(face.Point(numPnt)->Id());

					set<unsigned int>::iterator itFacePoint = pointFace.begin();
					bool foundFace = true;
					for(set<unsigned int>::iterator it = planePoints.begin(); it!= planePoints.end(); it++)
					{
						if(*itFacePoint != *it)
						{
							foundFace = false;
							break;
						}
						itFacePoint++;
					}
					if(foundFace && idDomain != -1)
					{
						GenericFace& faceModified = *mesh->Face(face.Id());
						faceModified.InitializeProperty(key.str());
						faceModified.AddProperty(key.str(), &DomainById(idDomain));
					}
					return Output::Success;
				}
			}
			return Output::Success;
		}

		for(unsigned int i = 0; i < cell.NumberOfEdges(); i++)
		{
			unsigned int idEdge = cell.Edge(i)->Id();
			GenericEdge& edge = *mesh->Edge(idEdge);

			CutterMesh3D::PositionInRelationToSurface positionOfEdge = EdgeDomainSplitter(edge);
			switch(positionOfEdge)
			{
				case Positive:
				{
					positiveEdges.insert(idEdge);
					break;
				}
				case Negative:
				{
					negativeEdges.insert(idEdge);
					break;
				}
				case ToBreak:
				{
					unsigned int idNewPoint;
					if(!edge.HasChilds())
						CutEdge(edge, idNewPoint);

					for(unsigned int cel = 0; cel < edge.NumberOfCells(); cel++)
					{
						if( edge.Cell(cel) != NULL)
						{
							const unsigned int& idCell = edge.Cell(cel)->Id();
							if(idCell != cell.Id())
								idCellsToUpdate.insert(idCell);
						}
					}
					planePoints.insert(idNewPoint);
					break;
				}
				case Plane:
				{
					planeEdges.insert(idEdge);
					break;
				}
				default:
					return Output::GenericError;
			}
		} //EdgesSplitted

		valueEdges.resize(planeEdges.size() + negativeEdges.size() + positiveEdges.size());
		unsigned int positionEdge = 0;
		for(set<unsigned int>::iterator it = positiveEdges.begin(); it != positiveEdges.end(); it++)
		{
			valueEdges[positionEdge] = PositionInRelationToSurface::Positive;
			idPositionEdges.insert(pair<unsigned int, unsigned int>(*it, positionEdge++));
		}
		for(set<unsigned int>::iterator it = negativeEdges.begin(); it != negativeEdges.end(); it++)
		{
			valueEdges[positionEdge] = PositionInRelationToSurface::Negative;
			idPositionEdges.insert(pair<unsigned int, unsigned int>(*it, positionEdge++));
		}
		for(set<unsigned int>::iterator it = planeEdges.begin(); it != planeEdges.end(); it++)
		{
			valueEdges[positionEdge] = PositionInRelationToSurface::Plane;
			idPositionEdges.insert(pair<unsigned int, unsigned int>(*it, positionEdge++));
		}

		for(unsigned int numFace = 0; numFace < cell.NumberOfFaces(); numFace++)
		{
			vector<const GenericPoint*> domainPointsInFace;
			unsigned int idFace = cell.Face(numFace)->Id();
			GenericFace& faceToCut = *mesh->Face(idFace);
			if(FaceDomainSplitter2(faceToCut, domainPointsInFace))
			{
				CutFace(faceToCut, domainPointsInFace);
				idFacesToUpdate.erase(idFace);
			}
			for(unsigned int cel = 0; cel < faceToCut.NumberOfCells(); cel++)
			{
				if( faceToCut.Cell(cel) != NULL)
				{
					const unsigned int& idCell = faceToCut.Cell(cel)->Id();
					if(idCell != cell.Id())
						idCellsToUpdate.insert(idCell);
				}
			}
		}//FacesSplitted

		for(set<unsigned int>::iterator it = idFacesToUpdate.begin(); it != idFacesToUpdate.end(); it++)
		{
			GenericFace& face = *mesh->Face(*it);
			if(face.InitializeProperty("OriginCell") == Output::Success)
			{
				const GenericTreeNode* node = &cell;
				while(node->HasFather())
					node = node->Father();
				face.AddProperty("OriginCell", new unsigned int(node->Id()));
			}
			mesh->UpdateFace(*it);
		}

		if(planePoints.size() != planeEdges.size())
		{
			Output::PrintErrorMessage("Number plane points not equal to number plane edges", false);
			return Output::GenericError;
		}

		if(!planeEdges.empty())
		{
			if(planeFaces.empty())
			{
				GenericFace& newFace = *(mesh->CreateFace());
				mesh->AddFace(&newFace);
				planeFaces.push_back(newFace.Id());

				if(newFace.InitializeProperty("OriginCell") == Output::Success)
				{
					const GenericTreeNode* node = &cell;
					while(node->HasFather())
						node = node->Father();
					newFace.AddProperty("OriginCell", new unsigned int(node->Id()));
				}

				newFace.InitializeEdges(planeEdges.size());
				if(idDomain != -1)
				{
					newFace.InitializeProperty(key.str());
					newFace.AddProperty(key.str(), &DomainById(idDomain));
				}

				for(set<unsigned int>::iterator it = planeEdges.begin(); it != planeEdges.end(); ++it)
				{
					GenericEdge& edge = *mesh->Edge(*it);
					newFace.AddEdge(&edge);
					edge.AddFace(&newFace);
					edge.ShrinkFaces();
				}

				//Adding points to newFace
				unsigned int numberPointNewFace = newFace.NumberOfEdges();
				newFace.InitializePoints(numberPointNewFace);
				newFace.AddPoint(newFace.Edge(0)->Point(0));
				newFace.AddPoint(newFace.Edge(0)->Point(1));
				unsigned int previousId = newFace.Edge(0)->Point(1)->Id();
				unsigned int j = 1;
				unsigned int i;

				while(newFace.NumberOfPoints() < newFace.NumberOfEdges() && j < newFace.NumberOfEdges())
				{
					i = j;
					while(i < newFace.NumberOfEdges())
					{
						if(newFace.Edge(i)->Point(0)->Id() ==  previousId)
						{
							newFace.AddPoint(newFace.Edge(i)->Point(1));
							previousId = newFace.Edge(i)->Point(1)->Id();

							const GenericEdge * tmpEdge = newFace.Edge(i);
							newFace.InsertEdge(newFace.Edge(j), i);
							newFace.InsertEdge(tmpEdge, j);
							break;
						}
						else if(newFace.Edge(i)->Point(1)->Id() ==  previousId)
						{
							newFace.AddPoint(newFace.Edge(i)->Point(0));
							previousId = newFace.Edge(i)->Point(0)->Id();

							const GenericEdge * tmpEdge = newFace.Edge(i);
							newFace.InsertEdge(newFace.Edge(j), i);
							newFace.InsertEdge(tmpEdge, j);
							break;
						}
						i++;
					}
					j++;
				}

				//checking if the orientation of the points of the face is positive with regard to the normal to the domain

				newFace.ComputeNormal();
				if(idDomain == -1)
				{
					Vector3d& normalFaceNewFace = newFace.Normal();
					unsigned int numEdgesNewFace = newFace.NumberOfEdges();
					if(normalFaceNewFace.dot(normalPlane) < - toll)	//the orientation is negative, therfore the vector of the points of the face is reversed
					{
						unsigned int numberPoints = numberPointNewFace;
						unsigned int forIteration =  ((numberPoints - 1) * 0.5);
						for(unsigned int i = 0; i < forIteration ; i++)
						{
							const GenericPoint* tmpPoint = newFace.Point(i + 1);
							newFace.InsertPoint(newFace.Point(numberPoints - i - 1), i + 1);
							newFace.InsertPoint(tmpPoint, numberPoints - i - 1);
						}
						for(unsigned int i = 0; i <= forIteration; i++)
						{
							const GenericEdge* tempEdge = newFace.Edge(i);
							newFace.InsertEdge(newFace.Edge(numberPoints - i - 1), i);
							newFace.InsertEdge(tempEdge, numberPoints - i - 1);
						}
						normalFaceNewFace *= -1.0;
					}

					for(unsigned int k = 0; k < numEdgesNewFace; k++)
					{
						const GenericPoint* tmpPoint = newFace.Point(k);
						const GenericPoint* tmpPointNext = newFace.Point((k + 1)%numEdgesNewFace);
						const GenericEdge* tempEdge = newFace.Edge(k);
						const GenericPoint* tmpPointEd = tempEdge->Point(0);
						const GenericPoint* tmpPointEdNext = tempEdge->Point(1);
						if((tmpPoint == tmpPointEd && tmpPointNext == tmpPointEdNext) ||
							 (tmpPoint == tmpPointEdNext && tmpPointNext == tmpPointEd))
							continue;
						else
						{
							for(unsigned int ed = k; ed < numEdgesNewFace; ed++)
							{
								const GenericEdge* tempEdge2 = newFace.Edge(ed);
								const GenericPoint* tmpPointEd = tempEdge2->Point(0);
								const GenericPoint* tmpPointEdNext = tempEdge2->Point(1);
								if((tmpPoint == tmpPointEd && tmpPointNext == tmpPointEdNext) ||
									 (tmpPoint == tmpPointEdNext && tmpPointNext == tmpPointEd))
								{
									newFace.InsertEdge(tempEdge2, k);
									newFace.InsertEdge(tempEdge, ed);
								}
								else
									continue;
							}
						}
					}
				}

				if(newFace.NumberOfPoints() != newFace.NumberOfEdges())
				{
					Output::PrintErrorMessage("Errror in adding of the points in the face %d", false, newFace.Id());
					return Output::GenericError;
				}

				for(set<unsigned int>::iterator it = planePoints.begin(); it != planePoints.end(); ++it )
				{
					GenericPoint& point = *mesh->Point(*it);
					point.AddFace(&newFace);
					point.ShrinkFaces();
				}

				cell.SetState(false);
				cell.InitializeChilds(2);

				//positive
				unsigned int positionCell = 1;
				GenericCell& newCell_1 = *(mesh->CreateCell());
				mesh->AddCell(&newCell_1);
				cell.AddChild(&newCell_1);

				newCell_1.SetFather(&cell);
				newCell_1.InheritPropertiesByFather();
				newCell_1.InitializeProperty("OriginCell");
				if(cell.HasProperty("OriginCell"))
					newCell_1.AddProperty("OriginCell", cell.GetProperty("OriginCell"));
				else
				{
					const GenericTreeNode* node = &cell;
					while(node->HasFather())
						node = node->Father();
					newCell_1.AddProperty("OriginCell", new unsigned int(node->Id()));
				}
				if(idDomain != -1)
				{
					newCell_1.InitializeProperty(key.str());
					newCell_1.AddProperty(key.str(), &DomainById(idDomain));
				}
				newCell_1.InitializeFaces(positiveFaces.size()+1);
				newCell_1.AllocateCells(positiveFaces.size()+1);

				//negative
				GenericCell& newCell_2 = *(mesh->CreateCell());
				mesh->AddCell(&newCell_2);
				cell.AddChild(&newCell_2);
				newCell_2.SetFather(&cell);
				newCell_2.InheritPropertiesByFather();
				newCell_2.InitializeProperty("OriginCell");

				if(cell.HasProperty("OriginCell"))
					newCell_2.AddProperty("OriginCell", cell.GetProperty("OriginCell"));
				else
					newCell_2.AddProperty("OriginCell", new int(cell.Id()));

				if(idDomain != -1)
				{
					newCell_1.InitializeProperty(key.str());
					newCell_2.AddProperty(key.str(), &DomainById(idDomain));
				}
				newCell_2.InitializeFaces(negativeFaces.size()+1);
				newCell_2.AllocateCells(negativeFaces.size()+1);

				newCell_2.AddFace(&newFace);
				newCell_1.AddFace(&newFace);

				newCell_1.InsertCell(&newCell_2, 0);
				newCell_2.InsertCell(&newCell_1, 0);

				for(list<unsigned int>::iterator it = positiveFaces.begin(); it != positiveFaces.end(); ++it )
				{
					GenericFace* face = mesh->Face(*it);
					newCell_1.AddFace(face);
					for(unsigned int numCell = 0; numCell < 2; numCell++)
					{
						if(face->Cell(numCell) != NULL)
						{
							if(face->Cell(numCell)->Id() == cell.Id() )
								face->InsertCell(&newCell_1, numCell);
							else
								newCell_1.InsertCell(face->Cell(numCell), positionCell);
						}
					}
					positionCell++;
				}

				newCell_1.InitializeEdges(positiveEdges.size() + planeEdges.size());
				for(set<unsigned int>::iterator it = positiveEdges.begin(); it != positiveEdges.end(); ++it)
				{
					GenericEdge& edge = *mesh->Edge(*it);
					newCell_1.AddEdge(&edge);
					for(unsigned int numCell = 0; numCell < edge.NumberOfCells(); numCell++)
					{
						if(edge.Cell(numCell)->Id() == cell.Id())
							edge.InsertCell(&newCell_1, numCell);
					}
				}

				for(set<unsigned int>::iterator it = positivePoints.begin(); it != positivePoints.end(); ++it)
				{
					GenericPoint& point = *mesh->Point(*it);
					newCell_1.AddPoint(&point);
					for(unsigned int numCell = 0; numCell < point.NumberOfCells(); numCell++)
					{
						if(point.Cell(numCell)->Id() == cell.Id())
							point.InsertCell(&newCell_1, numCell);
					}
				}

				positionCell = 1;
				for(list<unsigned int>::iterator it = negativeFaces.begin(); it != negativeFaces.end(); ++it )
				{
					GenericFace* face = mesh->Face(*it);
					newCell_2.AddFace(face);
					for(unsigned int numCell = 0; numCell < 2; numCell++)
					{
						if(face->Cell(numCell) != NULL)
						{
							if(face->Cell(numCell)->Id() == cell.Id() )
								face->InsertCell(&newCell_2, numCell);
							else
								newCell_2.InsertCell(face->Cell(numCell), positionCell);
						}
					}
					positionCell++;
				}

				newCell_2.InitializeEdges(negativeEdges.size() + planeEdges.size());
				for(set<unsigned int>::iterator it = negativeEdges.begin(); it != negativeEdges.end(); ++it)
				{
					GenericEdge& edge = *mesh->Edge(*it);
					newCell_2.AddEdge(&edge);
					for(unsigned int numCell = 0; numCell < edge.NumberOfCells(); numCell++)
					{
						if(edge.Cell(numCell)->Id() == cell.Id())
							edge.InsertCell(&newCell_2,numCell);
					}
				}

				for(set<unsigned int>::iterator it = negativePoints.begin();it!=negativePoints.end(); ++it )
				{
					GenericPoint& point = *mesh->Point(*it);
					newCell_2.AddPoint(&point);
					for(unsigned int numCell = 0; numCell < point.NumberOfCells(); numCell++)
					{
						if(point.Cell(numCell)->Id() == cell.Id())
							point.InsertCell(&newCell_2, numCell);
					}
				}

				for(set<unsigned int>::iterator it = planeEdges.begin(); it != planeEdges.end(); ++it)
				{
					GenericEdge& edge = *mesh->Edge(*it);
					newCell_2.AddEdge(&edge);
					newCell_1.AddEdge(&edge);
					edge.AddCell(&newCell_2);
					edge.AddCell(&newCell_1);
					edge.ShrinkCells();
				}

				for(set<unsigned int>::iterator it = planePoints.begin();it != planePoints.end(); ++it )
				{
					GenericPoint& point = *mesh->Point(*it);
					newCell_2.AddPoint(&point);
					newCell_1.AddPoint(&point);
					point.AddCell(&newCell_2);
					point.AddCell(&newCell_1);
					point.ShrinkCells();
				}

				newFace.ComputeGeometricalProperties();
				double inverseNumberPoints = 1.0/newCell_1.NumberOfPoints();
				double planeTranslation = newFace.Normal().dot(newFace.Centroid());
				Vector3d barycenter(0.0,0.0,0.0);
				for(unsigned int numCellPnt = 0; numCellPnt < newCell_1.NumberOfPoints(); numCellPnt++)
					barycenter += inverseNumberPoints * newCell_1.Point(numCellPnt)->Coordinates();

				newFace.AllocateCells(2);
				if((barycenter.dot(newFace.Normal()) - planeTranslation) > toll)
				{
					newFace.InsertCell(&newCell_2, 0);
					newFace.InsertCell(&newCell_1, 1);
				}
				else
				{
					newFace.InsertCell(&newCell_1, 0);
					newFace.InsertCell(&newCell_2, 1);
				}
			}
		}
		return Output::Success;
	}

	const Output::ExitCodes CutterMesh3D::Reset()
	{
		valuePoints.clear();
		idPositionPoints.clear();
		valueEdges.clear();
		idPositionEdges.clear();

		positivePoints.clear();
		negativePoints.clear();
		planePoints.clear();

		positiveEdges.clear();
		negativeEdges.clear();
		planeEdges.clear();

		positiveFaces.clear();
		negativeFaces.clear();
		planeFaces.clear();
		return Output::Success;
	}

	//****************************************************************************************************
	const CutterMesh3D::PositionInRelationToSurface CutterMesh3D::PointDomainSplitter(const GenericPoint& point, const Vector3d& normalPlane, const double& translation, const double& toll)
	{
		double r = normalPlane.transpose()*point.Coordinates() - translation;
		if(r  > toll)
			return PositionInRelationToSurface::Positive;
		else if(r  < -toll)
			return PositionInRelationToSurface::Negative;
		else
			return PositionInRelationToSurface::Plane;
	}

	//*************************
	const CutterMesh3D::PositionInRelationToSurface CutterMesh3D::EdgeDomainSplitter(const GenericEdge& edge)
	{
		const unsigned int& idZero = edge.Point(0)->Id();
		const unsigned int& idOne = edge.Point(1)->Id();

		unsigned int positionZero = idPositionPoints.at(idZero);
		unsigned int positionOne = idPositionPoints.at(idOne);

		PositionInRelationToSurface valueZero = valuePoints[positionZero];
		PositionInRelationToSurface valueOne = valuePoints[positionOne];

		bool pointP0 = valueZero == PositionInRelationToSurface::Positive;
		bool pointP1 = valueOne == PositionInRelationToSurface::Positive;
		bool pointN0 =  valueZero == PositionInRelationToSurface::Negative;
		bool pointN1 = valueOne == PositionInRelationToSurface::Negative;

		if( (pointP0	&& pointN1) ||	(pointN0 && pointP1))
			return ToBreak;
		else if( pointP0 || pointP1)
			return Positive;
		else if( pointN0 || pointN1)
			return Negative;
		else
			return Plane;

	}

	//*************************
	const bool CutterMesh3D::FaceDomainSplitter(GenericFace& face, vector<const GenericPoint*>& domainPointsInFace)
	{
		bool split = true;
		for(set<unsigned int>::iterator it = planePoints.begin(); it != planePoints.end(); it++)
		{
			const GenericPoint& point = *mesh->Point(*it);
			for(unsigned int numFace = 0; numFace < point.NumberOfFaces(); numFace++)
			{
				if(point.Face(numFace)->Id() == face.Id())
					domainPointsInFace.push_back(&point);
			}
		}

		if((domainPointsInFace.size() == 2))
		{
			bool domainThroughEdge = false;
			unsigned int numberOfEdges = face.NumberOfEdges();
			if(planeEdges.size() > 0)
				for(unsigned int i = 0; i < numberOfEdges; i++)
				{
					if(planeEdges.find(face.Edge(i)->Id()) != planeEdges.end())
					{
						domainThroughEdge = true;
						split = false;
						break;
					}
				}

			if(domainThroughEdge)
			{
				unsigned int j = 0;
				while(j < numberOfEdges && planeEdges.find(face.Edge(j)->Id()) != planeEdges.end())
					j++;
				unsigned int positionId = idPositionEdges.at(face.Edge(j)->Id());

				if(j == numberOfEdges)
					planeFaces.push_back(face.Id());
				else if(valueEdges[positionId] == PositionInRelationToSurface::Positive)
					positiveFaces.push_back(face.Id());
				else if(valueEdges[positionId] == PositionInRelationToSurface::Negative)
					negativeFaces.push_back(face.Id());

				idFacesToUpdate.erase(face.Id());

			}
			return split;
		}

		for(unsigned int numEdg = 0; numEdg < face.NumberOfEdges(); numEdg++)
		{
			if(face.Edge(numEdg)->IsActive())
			{
				unsigned int positionId = idPositionEdges.at(face.Edge(numEdg)->Id());
				if(valueEdges[positionId] == PositionInRelationToSurface::Positive)
					positiveFaces.push_back(face.Id());
				else if(valueEdges[positionId] == PositionInRelationToSurface::Negative)
					negativeFaces.push_back(face.Id());
				break;
			}
		}
		domainPointsInFace.clear();
		return false;
	}
	//*************************
	const bool CutterMesh3D::FaceDomainSplitter2(GenericFace& face, vector<const GenericPoint*>& domainPointsInFace)
	{
		bool split = true;
		bool domainThroughEdge = false;
		for(set<unsigned int>::iterator it = planePoints.begin(); it != planePoints.end(); it++)
		{
			const GenericPoint& point = *mesh->Point(*it);
			for(unsigned int numFace = 0; numFace < point.NumberOfFaces(); numFace++)
			{
				if(point.Face(numFace)->Id() == face.Id())
					domainPointsInFace.push_back(&point);
			}
		}

		if((domainPointsInFace.size() == 2))
		{
			unsigned int numberOfEdges = face.NumberOfEdges();
			if(planeEdges.size() > 0)
				for(unsigned int i = 0; i < numberOfEdges; i++)
				{
					if(planeEdges.find(face.Edge(i)->Id()) != planeEdges.end())
					{
						domainThroughEdge = true;
						split = false;
						break;
					}
				}

			if(domainThroughEdge)
			{
				unsigned int j = 0;
				while(j < numberOfEdges && planeEdges.find(face.Edge(j)->Id()) != planeEdges.end())
					j++;
				unsigned int positionId = idPositionEdges.at(face.Edge(j)->Id());

				if(j == numberOfEdges)
					planeFaces.push_back(face.Id());
				else if(valueEdges[positionId] == PositionInRelationToSurface::Positive)
					positiveFaces.push_back(face.Id());
				else if(valueEdges[positionId] == PositionInRelationToSurface::Negative)
					negativeFaces.push_back(face.Id());

				idFacesToUpdate.erase(face.Id());
			}
			return split;
		}

		for(unsigned int numEdg = 0; numEdg < face.NumberOfEdges(); numEdg++)
		{
			if(face.Edge(numEdg)->IsActive())
			{
				unsigned int positionId = idPositionEdges.at(face.Edge(numEdg)->Id());
				if(valueEdges[positionId] == PositionInRelationToSurface::Positive)
				{
					positiveFaces.push_back(face.Id());
					break;
				}
				else if(valueEdges[positionId] == PositionInRelationToSurface::Negative)
				{
					negativeFaces.push_back(face.Id());
					break;
				}
			}
		}
		domainPointsInFace.clear();
		return false;
	}
	// ***************************************************************************
	// CutterMesh3D Implementation
	// ***************************************************************************

	const Output::ExitCodes CutterMesh3D::CutMesh()
	{
		double inverseSize = 1.0/domains.size();
		for(unsigned int i = 0; i < domains.size(); i++)
		{
			GenericDomain2D& domain = static_cast<GenericDomain2D&>(*domains[i]);
			Output::PrintGenericMessageOnLine(" >> Test %f %% of domains", (i * 100) * inverseSize);
			Intersector2D1D& intersector = *intersector2D1D;
			intersector.SetPlane(domain.PlaneNormal(), domain.PlaneTranslation());

			for(unsigned int numCell = 0; numCell < mesh->NumberOfCells(); numCell++)
			{
				GenericCell& cell = *mesh->Cell(numCell);
				if(cell.IsActive())
					CutCell(cell, domain.PlaneNormal(), domain.PlaneTranslation(), domain.GlobalId());
				else
				{
					const GenericTreeNode* cellChild = cell.Child(0);
					while(cellChild->HasChilds())
						cellChild = cellChild->Child(0);

					unsigned int idChild = cellChild->Id();
					GenericCell& cellChildToCut = *mesh->Cell(idChild);
					CutCell(cellChildToCut, domain.PlaneNormal(), domain.PlaneTranslation(), domain.GlobalId());
				}

				for(set<unsigned int>::iterator itCellUpdate = idCellsToUpdate.begin(); itCellUpdate != idCellsToUpdate.end(); itCellUpdate++)
					mesh->UpdateCell(*itCellUpdate);

				idCellsToUpdate.clear();
				idFacesToUpdate.clear();
			}

		}
		mesh->CleanInactiveTreeNode();
		for(unsigned int numFac = 0; numFac < mesh->NumberOfFaces(); numFac++)
		{
			GenericFace& face = *mesh->Face(numFac);
			if(face.IsActive())
			{
				if(face.NumberOfCells() != 2)
				{
					Output::PrintErrorMessage("Error in Cutter", false);
					return Output::GenericError;
				}
			}
		}
		return Output::Success;
	}
}

