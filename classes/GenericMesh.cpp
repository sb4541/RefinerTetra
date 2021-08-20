#include "GenericMesh.hpp"
#include "Output.hpp"
#include <set>
#include "DefineNumbers.hpp"

using namespace MainApplication;

namespace GeDiM
{
	// ***************************************************************************
	GenericTreeNode::GenericTreeNode(const unsigned int& _id)
	{
		father = NULL;
		childs = vector<const GenericTreeNode*>{};

		id = _id;
		globalId = _id;
		isActive = true;

		markers = vector<unsigned int>{0};
	}

	GenericTreeNode::GenericTreeNode(const GenericTreeNode& treeNode)
	{
		father = NULL;
		childs.resize(treeNode.childs.size(), NULL); ///< Child nodes in the tree.
		for(unsigned int numChild = 0; numChild < treeNode.childs.size(); numChild++)
			childs[numChild] = treeNode.childs[numChild];

		id = treeNode.id; ///< Internal id of the node
		globalId = treeNode.globalId; ///< Global id of the node, generally equal to id
		isActive = treeNode.isActive; ///< Tells if the node is active for system matrix computation

		markers = treeNode.markers; ///< Marker for Dirichlet and Neumann conditions
		properties.insert(treeNode.properties.begin(), treeNode.properties.end());
	}
	GenericTreeNode::~GenericTreeNode()
	{
		properties.clear();
		childs.clear();
		father = NULL;
	}
	// ***************************************************************************
	Output::ExitCodes GenericTreeNode::InheritPropertiesByFather()
	{
		if(father == NULL)
			return Output::GenericError;

		const map< string, void* >& propertiesFather = father->GetAllProperties();
		if(propertiesFather.size() == 0)
			return Output::Success;

		for(map< string, void* >::const_iterator iteratorProperties = propertiesFather.begin(); iteratorProperties != propertiesFather.end(); iteratorProperties++)
			properties[iteratorProperties->first] = iteratorProperties->second;

		if(father->HasProperty("Cells"))
		{
			const vector<const GenericCell*>& cells = *static_cast< const vector<const GenericCell*>*>(father->GetProperty("Cells"));
			properties["Cells"]	= new vector<const GenericCell*>({cells[0], cells[1]});
		}

		return Output::Success;
	}

	// ***************************************************************************
	const Output::ExitCodes GenericTreeNode::InitializeProperty(const string& key)
	{
		map<string, void*>::iterator finder;
		finder = properties.find(key);
		if(finder == properties.end())
		{
			pair <string, void*> property(key, NULL);
			properties.insert(property);
			return Output::Success;
		}
		else
			return Output::GenericError;
	}

	// ***************************************************************************
	GenericCell::GenericCell(const unsigned int& _id) : GenericTreeNode(_id)
	{
	}

	GenericCell::GenericCell(const GenericCell& cell) : GenericTreeNode(cell)
	{
		edges.resize(cell.edges.size(), NULL);
		faces.resize(cell.faces.size(), NULL);
		cells.resize(cell.cells.size(), NULL);
		points.resize(cell.points.size(), NULL);
		childs.resize(cell.childs.size(), NULL);
		father = NULL;
	}
	GenericCell::~GenericCell()
	{
		//		DestructorProperties();
		points.clear();
		edges.clear();
		faces.clear();
		cells.clear();
	}
	// ***************************************************************************
	Output::ExitCodes GenericCell::AddCell(const GenericCell* cell)
	{
		if (cell == NULL)
			return Output::GenericError;

		cells.push_back(cell);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericCell::AddFace(const GenericFace* face)
	{
		if (face == NULL)
			return Output::GenericError;

		faces.push_back(face);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericCell::AddEdge(const GenericEdge* edge)
	{
		if (edge == NULL)
			return Output::GenericError;

		edges.push_back(edge);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericCell::AddPoint(const GenericPoint* point)
	{
		if (point == NULL)
			return Output::GenericError;

		points.push_back(point);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericCell::InsertCell(const GenericCell* cell, const unsigned int& position)
	{
		if (cell == NULL || position >= cells.size())
			return Output::GenericError;

		cells[position] = cell;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericCell::InsertFace(const GenericFace* face, const unsigned int& position)
	{
		if (face == NULL || position >= faces.size())
			return Output::GenericError;

		faces[position] = face;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericCell::InsertEdge(const GenericEdge* edge, const unsigned int& position)
	{
		if (edge == NULL || position >= edges.size())
			return Output::GenericError;

		edges[position] = edge;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericCell::InsertPoint(const GenericPoint* point, const unsigned int& position)
	{
		if (point == NULL || position >= points.size())
			return Output::GenericError;

		points[position] = point;

		return Output::Success;
	}
	// ***************************************************************************
	bool GenericCell::PointInCellAndIdEdgeBoundary (const Vector3d& point, int& numEdges, const double& toll) const
	{
		for(unsigned int pntCell = 0; pntCell < NumberOfPoints(); pntCell++)
		{
			const GenericPoint* pointFirst = Point(pntCell);
			const GenericPoint* pointSecond = Point((pntCell+1)%NumberOfPoints());

			Vector3d tangentVectorEdge = pointSecond->Coordinates() - pointFirst->Coordinates();
			Vector3d tangentVectorDifference = point - pointFirst->Coordinates();

			double crossProduct = tangentVectorEdge.x() * tangentVectorDifference.y() - tangentVectorDifference.x() * tangentVectorEdge.y();
			if( crossProduct > -toll && crossProduct < toll)
			{
				if((tangentVectorEdge.x() * tangentVectorDifference.x() < -toll) || (tangentVectorEdge.y() * tangentVectorDifference.y() < -toll))
					continue;
				if(tangentVectorEdge.squaredNorm() < tangentVectorDifference.squaredNorm())
					continue;
				numEdges = pntCell;
			}

			if( crossProduct < -toll)
				return false;
		}
		return true;
	}

	// ***************************************************************************
	bool GenericCell::PointInCell (const Vector3d& point, const double& toll) const
	{
		for(unsigned int pntCell = 0; pntCell < NumberOfPoints(); pntCell++)
		{
			const GenericPoint* pointFirst = Point(pntCell);
			const GenericPoint* pointSecond = Point((pntCell+1)%NumberOfPoints());

			Vector3d tangentVectorEdge = pointSecond->Coordinates() - pointFirst->Coordinates();
			Vector3d tangentVectorDifference = point - pointFirst->Coordinates();

			double crossProduct = tangentVectorEdge.x() * tangentVectorDifference.y() - tangentVectorDifference.x() * tangentVectorEdge.y();

			if( crossProduct < -toll)
				return false;
		}
		return true;
	}
	// ***************************************************************************
	const unsigned int GenericCell::RatioMaxMinEdge() const
	{
		double edgeMinCell = numeric_limits<double>::max();
		double edgeMaxCell = 0.0;
		for(unsigned int numEdge = 0; numEdge < NumberOfEdges(); numEdge++)
		{
			const GenericEdge& edge = *Edge(numEdge);
			const Vector3d& firstPoint = edge.Point(0)->Coordinates();
			const Vector3d& secondPoint = edge.Point(1)->Coordinates();
			double length = (secondPoint - firstPoint).norm();
			if(length < edgeMinCell)
				edgeMinCell = length;
			if(length > edgeMaxCell)
				edgeMaxCell = length;
		}
		return edgeMaxCell / edgeMinCell;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericCell::Centroid(Vector3d& centroid) const
	{
		centroid.setZero(3);
		unsigned int numberPoints =  NumberOfPoints();

		if(faces.size() == 0)
		{
			vector<double> valuePositive;
			vector<double> valueNegative;
			double valueP = 0.0;
			double valueN = 0.0;
			vector<double> valuePositive2;
			vector<double> valueNegative2;
			double value2P = 0.0;
			double value2N = 0.0;
			vector<double> measurePositive;
			vector<double> measureNegative;
			double measureP = 0.0;
			double measureN = 0.0;

			valuePositive.reserve(numberPoints);
			valueNegative.reserve(numberPoints);
			valuePositive2.reserve(numberPoints);
			valueNegative2.reserve(numberPoints);
			measurePositive.reserve(numberPoints);
			measureNegative.reserve(numberPoints);
			double resultTemp = 0.0;

			for (unsigned int i = 0; i < numberPoints; i++)
			{
				const Vector3d& vertex = Point(i)->Coordinates();
				const Vector3d& nextVertex = Point((i + 1)%numberPoints)->Coordinates();
				Vector3d sumPoint = (vertex + nextVertex);
				Vector3d normal(0.5*(nextVertex[1] - vertex[1]), 0.5*(vertex[0] - nextVertex[0]), 0.0);

				resultTemp = sumPoint(0) * normal(0);
				if(resultTemp > 1.0E-16)
					measurePositive.push_back(resultTemp);
				else
					measureNegative.push_back(resultTemp);

				resultTemp = (vertex(0)*vertex(0) +
											sumPoint(0)*sumPoint(0)	+
											nextVertex(0)*nextVertex(0)) * normal(0);

				if(resultTemp > 1.0E-16)
					valuePositive.push_back(resultTemp);
				else
					valueNegative.push_back(resultTemp);

				resultTemp = (vertex(1)*vertex(1) +
											sumPoint(1)*sumPoint(1)	+
											nextVertex(1)*nextVertex(1)) * normal(1);

				if(resultTemp > 1.0E-16)
					valuePositive2.push_back(resultTemp);
				else
					valueNegative2.push_back(resultTemp);
			}
			sort(valuePositive.begin(), valuePositive.end());
			sort(valueNegative.begin(), valueNegative.end());
			sort(valuePositive2.begin(), valuePositive2.end());
			sort(valueNegative2.begin(), valueNegative2.end());
			sort(measurePositive.begin(), measurePositive.end());
			sort(measureNegative.begin(), measureNegative.end());

			for_each(valuePositive.begin(), valuePositive.end(), [&] (double n) {valueP += n;});
			for_each(valueNegative.begin(), valueNegative.end(), [&] (double n) {valueN += n;});
			for_each(valuePositive2.begin(), valuePositive2.end(), [&] (double n) {value2P += n;});
			for_each(valueNegative2.begin(), valueNegative2.end(), [&] (double n) {value2N += n;});
			for_each(measurePositive.begin(), measurePositive.end(), [&] (double n) {measureP += n;});
			for_each(measureNegative.begin(), measureNegative.end(), [&] (double n) {measureN += n;});

			centroid[0] = valueP + valueN ;
			centroid[1] = value2P + value2N;

			double invertMeasure = 1.0/(measureP + measureN);
			centroid[0] = centroid[0] * invertMeasure * RAT1_6;
			centroid[1] = centroid[1] * invertMeasure * RAT1_6;
		}
		else
		{

			double measure = 0.0;
			for(unsigned int numFac = 0; numFac < NumberOfFaces(); numFac++)
			{
				const GenericFace& face = *faces[numFac];
				if(!face.CheckGeometricalProperties())
				{
					Output::PrintErrorMessage("Compute Faces Geometrical Properties", false);
					return Output::GenericError;
				}

				unsigned int numberVerticesFace = face.NumberOfPoints();
				const Vector3d& centroidFace = face.Centroid();
				Vector3d normal = face.DirectionNormal(id) * face.Normal();

				MatrixXd matrixB(3,2);
				for(unsigned int numPnt = 0; numPnt < numberVerticesFace; numPnt++)
				{
					const Vector3d& PointI = face.Point(numPnt)->Coordinates();
					const Vector3d& PointI_1 = face.Point((numPnt+1)%numberVerticesFace)->Coordinates();
					Vector3d tangentFirst = PointI_1 - PointI;
					Vector3d tangentSecond = centroidFace - PointI;
					matrixB.col(0) = tangentFirst;
					matrixB.col(1) = tangentSecond;

					double doubleAreaSubTriangles = (tangentFirst.cross(tangentSecond)).norm();
					Vector2d first(RAT1_6, RAT1_6);
					Vector2d second(RAT2_3, RAT1_6);
					Vector2d third(RAT1_6, RAT2_3);
					Vector3d firstPoint = PointI + matrixB * first;
					Vector3d secondPoint = PointI + matrixB * second;
					Vector3d thirdPoint = PointI + matrixB * third;

					Vector3d point =  PointI + matrixB * Vector2d(RAT1_3, RAT1_3);

					measure += point(0) * normal(0) * doubleAreaSubTriangles * 0.5;
					centroid[0] += normal(0) * doubleAreaSubTriangles * (firstPoint(0) * firstPoint(0) + secondPoint(0)*secondPoint(0) + thirdPoint(0)*thirdPoint(0));
					centroid[1] += normal(1) * doubleAreaSubTriangles * (firstPoint(1) * firstPoint(1) + secondPoint(1)*secondPoint(1) + thirdPoint(1)*thirdPoint(1));
					centroid[2] += normal(2) * doubleAreaSubTriangles * (firstPoint(2) * firstPoint(2) + secondPoint(2)*secondPoint(2) + thirdPoint(2)*thirdPoint(2));
				}
			}
			centroid *= 0.5 * RAT1_6 / measure;
		}
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericCell::Radius(double& radius, const Vector3d& centroid, const double squaredTolerance) const
	{

		double squaredRadius = 0.0;
		radius = 0.0;
		const size_t numberVertices = NumberOfPoints();
		for (unsigned int v = 0; v < numberVertices; v++)
		{
			double squaredDistance = (centroid - Point(v)->Coordinates()).squaredNorm();
			if (squaredRadius < squaredDistance)
				squaredRadius = squaredDistance;
		}
		radius = sqrt(squaredRadius);

		if (fabs(squaredRadius) < squaredTolerance)
		{
			squaredRadius = 0.0;
			Output::PrintErrorMessage("Cell %d. Squared Radius too small: %.2e - tolerance %.2e", false, globalId, squaredRadius, squaredTolerance);
			return Output::GenericError;
		}
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericCell::Measure(double& measure) const
	{
		unsigned int numberPoints =  NumberOfPoints();
		vector<double> measurePositive;
		vector<double> measureNegative;
		double measureP = 0.0;
		double measureN = 0.0;

		measurePositive.reserve(numberPoints);
		measureNegative.reserve(numberPoints);
		double resultTemp = 0.0;

		for (unsigned int i = 0; i < numberPoints; i++)
		{
			const Vector3d& vertex = Point(i)->Coordinates();
			const Vector3d& nextVertex = Point((i + 1)%numberPoints)->Coordinates();
			Vector3d sumPoint = (vertex + nextVertex);
			Vector3d normal(0.5*(nextVertex[1] - vertex[1]), 0.5*(vertex[0] - nextVertex[0]), 0.0);

			resultTemp = sumPoint(0) * normal(0);
			if(resultTemp > 1.0E-16)
				measurePositive.push_back(resultTemp);
			else
				measureNegative.push_back(resultTemp);
		}
		sort(measurePositive.begin(), measurePositive.end());
		sort(measureNegative.begin(), measureNegative.end());

		for_each(measurePositive.begin(), measurePositive.end(), [&] (double n) {measureP += n;});
		for_each(measureNegative.begin(), measureNegative.end(), [&] (double n) {measureN += n;});

		measure = measureP + measureN;
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericCell::MeasureCentroid(double& measure, Vector3d& centroid) const
	{
		centroid.setZero(3);
		unsigned int numberPoints =  NumberOfPoints();
		vector<double> valuePositive;
		vector<double> valueNegative;
		double valueP = 0.0;
		double valueN = 0.0;
		vector<double> valuePositive2;
		vector<double> valueNegative2;
		double value2P = 0.0;
		double value2N = 0.0;
		vector<double> measurePositive;
		vector<double> measureNegative;
		double measureP = 0.0;
		double measureN = 0.0;

		valuePositive.reserve(numberPoints);
		valueNegative.reserve(numberPoints);
		valuePositive2.reserve(numberPoints);
		valueNegative2.reserve(numberPoints);
		measurePositive.reserve(numberPoints);
		measureNegative.reserve(numberPoints);
		double resultTemp = 0.0;

		for (unsigned int i = 0; i < numberPoints; i++)
		{
			const Vector3d& vertex = Point(i)->Coordinates();
			const Vector3d& nextVertex = Point((i + 1)%numberPoints)->Coordinates();
			Vector3d sumPoint = (vertex + nextVertex);
			Vector3d normal(0.5*(nextVertex[1] - vertex[1]), 0.5*(vertex[0] - nextVertex[0]), 0.0);

			resultTemp = sumPoint(0) * normal(0);
			if(resultTemp > 1.0E-16)
				measurePositive.push_back(resultTemp);
			else
				measureNegative.push_back(resultTemp);

			resultTemp = (vertex(0)*vertex(0) +
										sumPoint(0)*sumPoint(0)	+
										nextVertex(0)*nextVertex(0)) * normal(0);

			if(resultTemp > 1.0E-16)
				valuePositive.push_back(resultTemp);
			else
				valueNegative.push_back(resultTemp);

			resultTemp = (vertex(1)*vertex(1) +
										sumPoint(1)*sumPoint(1)	+
										nextVertex(1)*nextVertex(1)) * normal(1);

			if(resultTemp > 1.0E-16)
				valuePositive2.push_back(resultTemp);
			else
				valueNegative2.push_back(resultTemp);
		}
		sort(valuePositive.begin(), valuePositive.end());
		sort(valueNegative.begin(), valueNegative.end());
		sort(valuePositive2.begin(), valuePositive2.end());
		sort(valueNegative2.begin(), valueNegative2.end());
		sort(measurePositive.begin(), measurePositive.end());
		sort(measureNegative.begin(), measureNegative.end());

		for_each(valuePositive.begin(), valuePositive.end(), [&] (double n) {valueP += n;});
		for_each(valueNegative.begin(), valueNegative.end(), [&] (double n) {valueN += n;});
		for_each(valuePositive2.begin(), valuePositive2.end(), [&] (double n) {value2P += n;});
		for_each(valueNegative2.begin(), valueNegative2.end(), [&] (double n) {value2N += n;});
		for_each(measurePositive.begin(), measurePositive.end(), [&] (double n) {measureP += n;});
		for_each(measureNegative.begin(), measureNegative.end(), [&] (double n) {measureN += n;});

		centroid[0] = valueP + valueN ;
		centroid[1] = value2P + value2N;
		measure = measureP + measureN;
		double invertMeasure = 1.0/(measureP + measureN);
		centroid[0] = centroid[0] * invertMeasure * RAT1_6;
		centroid[1] = centroid[1] * invertMeasure * RAT1_6;
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericCell::ComputeGeometricalProperties(Vector3d& normal, Matrix3d& rotationMatrix, const double& rotationTolerance) const
	{
		if(faces.size() != 0)
		{
			Output::PrintErrorMessage("Is not possible to compute geometrical properties for 3D cells", false);
			return Output::GenericError;
		}

		unsigned int totalNumberVertices = points.size();
		normal.setZero();
		for(unsigned int i = 0; i < totalNumberVertices; i++)
		{
			int previousEdgeIndex = (i == 0) ? totalNumberVertices - 1 : i - 1;

			Vector3d edge = points[(i+1) % totalNumberVertices]->Coordinates() - points[i]->Coordinates();
			Vector3d edgePrevious = points[previousEdgeIndex]->Coordinates() - points[i]->Coordinates();

			normal.noalias() += edge.cross(edgePrevious);
		}
		normal.normalize();

		MatrixXd Z(3, totalNumberVertices);
		MatrixXd W(3, totalNumberVertices);
		Matrix3d H;
		Vector3d V1mV0 = points[1]->Coordinates() - points[0]->Coordinates();
		double normVectorOne = V1mV0.norm();
		Z.col(0) = V1mV0;
		W.col(0) << normVectorOne * cos(0.0), normVectorOne * sin(0.0), 0;
		for (unsigned int i = 2; i < totalNumberVertices; i++)
		{
			Vector3d VimV0 = points[i]->Coordinates() - points[0]->Coordinates();
			Z.col(i - 1) = VimV0;

			double normVectorI = VimV0.norm();
			double prodScalar = VimV0.dot(V1mV0) / (normVectorOne * normVectorI);
			double angleBetweenVectors = 0.0;
			if(prodScalar > 1.0 - rotationTolerance)
			{
				angleBetweenVectors = 0.0;
				W.col(i - 1) << normVectorI, 0.0, 0.0;
			}
			else if(prodScalar < -1.0 + rotationTolerance)
			{
				angleBetweenVectors = M_PI;
				W.col(i - 1) << -normVectorI, 0.0, 0.0;
			}
			else if((prodScalar > -rotationTolerance) && (prodScalar < rotationTolerance))
			{
				angleBetweenVectors = M_PI * 0.5;
				W.col(i - 1) << 0.0, normVectorI, 0.0;
			}
			else
			{
				angleBetweenVectors = acos(prodScalar);
				W.col(i - 1) << normVectorI * cos(angleBetweenVectors), normVectorI * sin(angleBetweenVectors), 0;
			}
		}
		Z.col(totalNumberVertices - 1) = normal;
		W.col(totalNumberVertices - 1)<< 0, 0, 1;
		H = W * Z.transpose();
		JacobiSVD<MatrixXd> svd(H, ComputeThinU | ComputeThinV);
		rotationMatrix = svd.matrixV() * (svd.matrixU()).transpose();
		return Output::Success;
	}
	//****************************************************************************************************
	void GenericCell::AllignedEdgesPoints()
	{
		unsigned int numberOfEdges = NumberOfEdges(); // equal to number of points
		vector<const GenericEdge*> positionEdges(numberOfEdges);
		unsigned int firstIdPointCell = Point(0)->Id();
		unsigned int secondIdPointCell = Point(1)->Id();
		bool check = true;

		for(unsigned int edg = 0; edg < numberOfEdges; edg++)
		{
			const GenericEdge& edgeTemp = *Edge(edg);
			unsigned int idFirstPoint = edgeTemp.Point(0)->Id();
			unsigned int idSecondPoint = edgeTemp.Point(1)->Id();
			if((firstIdPointCell == idFirstPoint && secondIdPointCell == idSecondPoint) || (firstIdPointCell == idSecondPoint && secondIdPointCell == idFirstPoint))
			{
				if(edg == 0)
				{
					check = false;
					break;
				}

				positionEdges[0] = &edgeTemp;
				for(unsigned int i = 1; i < numberOfEdges; i++)
					positionEdges[i] = Edge((edg+i)%numberOfEdges);
				break;
			}
		}

		if(check)
			for(unsigned int edg = 0; edg < numberOfEdges; edg++)
				InsertEdge(positionEdges[edg], edg);
		return;
	}
	// ***************************************************************************
	GenericFace::GenericFace(const unsigned int& _id) : GenericTreeNode(_id)
	{
		normal = NULL;
		centroid = NULL;
		rotationMatrix = NULL;
	}

	GenericFace::GenericFace(const GenericFace& face): GenericTreeNode(face)
	{
		edges.resize(face.edges.size(), NULL);
		faces.resize(face.faces.size(), NULL);
		cells.resize(face.cells.size(), NULL);
		points.resize(face.points.size(), NULL);
		childs.resize(face.childs.size(), NULL);
		father = NULL;

		if(face.CheckGeometricalProperties())
		{
			normal = face.normal;
			centroid = face.centroid;
			rotationMatrix = face.rotationMatrix;
		}
	}
	GenericFace::~GenericFace()
	{
		//		DestructorProperties();
		points.clear();
		edges.clear();
		faces.clear();
		cells.clear();

		for(unsigned int numPnt = 0; numPnt < points.size(); numPnt++)
		{
			if(rotatedPoints[numPnt] != NULL)
				delete rotatedPoints[numPnt];
		}
		rotatedPoints.clear();

		if(normal != NULL)
			delete normal;

		if(centroid != NULL)
			delete centroid;

		if(rotationMatrix != NULL)
			delete rotationMatrix;
	}
	// ***************************************************************************
	const int GenericFace::DirectionNormal(unsigned int idCell) const
	{
		const GenericCell* cellTemp;
		if(cells[1] != NULL && cells[1]->Id() == idCell)
			cellTemp = cells[1];
		else if(cells[0] != NULL && cells[0]->Id() == idCell)
			cellTemp = cells[0];
		else
		{
			Output::PrintErrorMessage("Face %d has not cells neighs ", false, this->Id());
		}

		double inverseNumberPoints = 1.0/cellTemp->NumberOfPoints();
		Vector3d barycenter;
		for(unsigned int numCellPnt = 0; numCellPnt < cellTemp->NumberOfPoints(); numCellPnt++)
			barycenter += inverseNumberPoints * cellTemp->Point(numCellPnt)->Coordinates();

		double planeTranslation = (*normal).dot(*centroid);
		if((barycenter.dot(*normal) - planeTranslation) > 1.0E-16)
			return -1;
		else
			return 1;
	}
	// ***************************************************************************
	Output::ExitCodes GenericFace::AddCell(const GenericCell* cell)
	{
		if (cell == NULL)
			return Output::GenericError;

		cells.push_back(cell);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericFace::AddFace(const GenericFace* face)
	{
		if (face == NULL)
			return Output::GenericError;

		faces.push_back(face);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericFace::AddEdge(const GenericEdge* edge)
	{
		if (edge == NULL)
			return Output::GenericError;

		edges.push_back(edge);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericFace::AddPoint(const GenericPoint* point)
	{
		if (point == NULL)
			return Output::GenericError;

		points.push_back(point);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericFace::InsertCell(const GenericCell* cell, const unsigned int& position)
	{
		if (cell == NULL || position >= cells.size())
			return Output::GenericError;

		cells[position] = cell;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericFace::InsertFace(const GenericFace* face, const unsigned int& position)
	{
		if (face == NULL || position >= faces.size())
			return Output::GenericError;

		faces[position] = face;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericFace::InsertEdge(const GenericEdge* edge, const unsigned int& position)
	{
		if (edge == NULL || position >= edges.size())
			return Output::GenericError;

		edges[position] = edge;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericFace::InsertPoint(const GenericPoint* point, const unsigned int& position)
	{
		if (point == NULL || position >= points.size())
			return Output::GenericError;

		points[position] = point;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericFace::InitializeGeometricalProperties()
	{
		measure = 0.0;
		rotatedPoints.resize(points.size());
		for(unsigned int numPnt = 0; numPnt < points.size(); numPnt++)
			rotatedPoints[numPnt] = new Vector3d();
		centroid = new Vector3d();
		rotationMatrix = new Matrix3d();

		return Output::Success;
	}

	// ***************************************************************************
	Output::ExitCodes GenericFace::ComputeNormal()
	{
		unsigned int totalNumberVertices = points.size();
		normal = new Vector3d();
		Vector3d& normalRef = *normal;
		normalRef.setZero();
		for(unsigned int i = 0; i < totalNumberVertices; i++)
		{
			int previousEdgeIndex = (i == 0) ? totalNumberVertices - 1 : i - 1;

			Vector3d edge = points[(i+1) % totalNumberVertices]->Coordinates() - points[i]->Coordinates();
			Vector3d edgePrevious = points[previousEdgeIndex]->Coordinates() - points[i]->Coordinates();

			normalRef.noalias() += edge.cross(edgePrevious);
		}
		normalRef.normalize();
		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericFace::ComputeGeometricalProperties(const double rotationTolerance)
	{
		if(centroid == NULL)
		{
			InitializeGeometricalProperties();
		}
		unsigned int totalNumberVertices = points.size();

		if(normal == NULL)
			ComputeNormal();

		Vector3d& normalRef = *normal;
		MatrixXd Z(3, totalNumberVertices);
		MatrixXd W(3, totalNumberVertices);
		Matrix3d H;
		Z.setZero();
		W.setZero();
		H.setZero();
		Vector3d V1mV0 = points[1]->Coordinates() - points[0]->Coordinates();
		double normVectorOne = V1mV0.norm();
		Z.col(0) = V1mV0;
		W.col(0) << normVectorOne * cos(0.0), normVectorOne * sin(0.0), 0;
		for (unsigned int i = 2; i < totalNumberVertices; i++)
		{
			Vector3d VimV0 = points[i]->Coordinates() - points[0]->Coordinates();
			Z.col(i - 1) = VimV0;

			double normVectorI = VimV0.norm();
			double prodScalar = VimV0.dot(V1mV0) / (normVectorOne * normVectorI);
			double angleBetweenVectors = 0.0;

			if((prodScalar > -1.0 + rotationTolerance) && (prodScalar < 1.0 - rotationTolerance))
			{
				angleBetweenVectors = acos(prodScalar);
			}
			else if(prodScalar > 1.0 - rotationTolerance)
			{
				angleBetweenVectors = 0.0;
			}
			else
			{
				angleBetweenVectors = M_PI;
			}

			W.col(i - 1) << normVectorI * cos(angleBetweenVectors), normVectorI * sin(angleBetweenVectors), 0;
		}
		Z.col(totalNumberVertices - 1) = normalRef;
		W.col(totalNumberVertices - 1)<< 0, 0, 1;
		H = W * Z.transpose();
		JacobiSVD<MatrixXd> svd(H, ComputeThinU | ComputeThinV);
		Matrix3d& rotMatrix = *rotationMatrix;
		rotMatrix = svd.matrixV() * (svd.matrixU()).transpose();
		Vector3d translation = points[0]->Coordinates();


		vector<double> valuePositive;
		vector<double> valueNegative;
		double valueP = 0.0;
		double valueN = 0.0;
		vector<double> valuePositive2;
		vector<double> valueNegative2;
		double value2P = 0.0;
		double value2N = 0.0;
		vector<double> measurePositive;
		vector<double> measureNegative;

		measure = 0.0;
		double measureP = 0.0;
		double measureN = 0.0;

		unsigned int numberVerticesFace = NumberOfPoints();

		valuePositive.reserve(numberVerticesFace);
		valueNegative.reserve(numberVerticesFace);
		valuePositive2.reserve(numberVerticesFace);
		valueNegative2.reserve(numberVerticesFace);
		measurePositive.reserve(numberVerticesFace);
		measureNegative.reserve(numberVerticesFace);
		double resultTemp = 0.0;

		for (unsigned int i = 0; i < numberVerticesFace; i++)
		{
			Vector3d& rotatedPoint = *rotatedPoints[i];
			rotatedPoint.setZero();
			rotatedPoint = rotMatrix.transpose() * (points[i]->Coordinates() - translation);
		}

		Vector3d& centroidTemp = *centroid;
		centroidTemp.setZero();
		for (unsigned int i = 0; i < numberVerticesFace; i++)
		{
			const Vector3d& vertex = *rotatedPoints[i];
			const Vector3d& nextVertex = *rotatedPoints[(i+1)%numberVerticesFace];
			Vector3d sumPoint = (vertex + nextVertex);
			Vector3d normalTemp(0.5*(nextVertex[1] - vertex[1]), 0.5*(vertex[0] - nextVertex[0]), 0.0);

			resultTemp = sumPoint(0) * normalTemp(0);
			if(resultTemp > 1.0E-16)
				measurePositive.push_back(resultTemp);
			else
				measureNegative.push_back(resultTemp);

			resultTemp = (vertex(0)*vertex(0) +
										sumPoint(0)*sumPoint(0)	+
										nextVertex(0)*nextVertex(0)) * normalTemp(0);

			if(resultTemp > 1.0E-16)
				valuePositive.push_back(resultTemp);
			else
				valueNegative.push_back(resultTemp);

			resultTemp = (vertex(1)*vertex(1) +
										sumPoint(1)*sumPoint(1)	+
										nextVertex(1)*nextVertex(1)) * normalTemp(1);

			if(resultTemp > 1.0E-16)
				valuePositive2.push_back(resultTemp);
			else
				valueNegative2.push_back(resultTemp);
		}
		sort(valuePositive.begin(), valuePositive.end());
		sort(valueNegative.begin(), valueNegative.end());
		sort(valuePositive2.begin(), valuePositive2.end());
		sort(valueNegative2.begin(), valueNegative2.end());
		sort(measurePositive.begin(), measurePositive.end());
		sort(measureNegative.begin(), measureNegative.end());

		for_each(valuePositive.begin(), valuePositive.end(), [&] (double n) {valueP += n;});
		for_each(valueNegative.begin(), valueNegative.end(), [&] (double n) {valueN += n;});
		for_each(valuePositive2.begin(), valuePositive2.end(), [&] (double n) {value2P += n;});
		for_each(valueNegative2.begin(), valueNegative2.end(), [&] (double n) {value2N += n;});
		for_each(measurePositive.begin(), measurePositive.end(), [&] (double n) {measureP += n;});
		for_each(measureNegative.begin(), measureNegative.end(), [&] (double n) {measureN += n;});

		centroidTemp[0] = valueP + valueN ;
		centroidTemp[1] = value2P + value2N;

		measure = (measureP + measureN);

		double invertMeasure = 1.0/measure;
		centroidTemp[0] = centroidTemp[0] * invertMeasure * RAT1_6;
		centroidTemp[1] = centroidTemp[1] * invertMeasure * RAT1_6;
		if(measure < 1.0E-16)
		{
			measure = -measure;
		}

		centroidTemp = rotMatrix * centroidTemp + translation;

		//		double planeTranslation = normalRef.dot(points[0]->Coordinates());
		//		const GenericCell& cellOne = *cells[1];
		//		double inverseNumberPoints = 1.0/cellOne.NumberOfPoints();
		//		Vector3d barycenter;
		//		for(unsigned int numCellPnt = 0; numCellPnt < cellOne.NumberOfPoints(); numCellPnt++)
		//			barycenter += inverseNumberPoints * cellOne.Point(numCellPnt)->Coordinates();

		//		if((barycenter.dot(normalRef) - planeTranslation) > 1.0E-16)
		//			normalRef *= -1.0;

		return Output::Success;
	}
	// ***************************************************************************
	GenericEdge::GenericEdge(const unsigned int& _id) : GenericTreeNode(_id)
	{
		points.reserve(2);
	}
	GenericEdge::GenericEdge(const unsigned int& _id, const GenericPoint* origin, const GenericPoint* end) : GenericEdge(_id)
	{
		AddPoint(origin);
		AddPoint(end);
	}

	GenericEdge::GenericEdge(const GenericEdge& edge) : GenericTreeNode(edge)
	{
		edges.resize(edge.edges.size(), NULL);
		faces.resize(edge.faces.size(), NULL);
		cells.resize(edge.cells.size(), NULL);
		points.resize(2, NULL);
		childs.resize(edge.childs.size(), NULL);
		father = NULL;
	}

	GenericEdge::~GenericEdge()
	{
		//		DestructorProperties();
		points.clear();
		edges.clear();
		faces.clear();
		cells.clear();
	}
	// ***************************************************************************
	Output::ExitCodes GenericEdge::AddCell(const GenericCell* cell)
	{
		if (cell == NULL)
			return Output::GenericError;

		cells.push_back(cell);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericEdge::AddFace(const GenericFace* face)
	{
		if (face == NULL)
			return Output::GenericError;

		faces.push_back(face);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericEdge::AddEdge(const GenericEdge* edge)
	{
		if (edge == NULL)
			return Output::GenericError;

		edges.push_back(edge);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericEdge::AddPoint(const GenericPoint* point)
	{
		if (point == NULL)
			return Output::GenericError;

		points.push_back(point);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericEdge::InsertCell(const GenericCell* cell, const unsigned int& position)
	{
		if (cell == NULL || position >= cells.size())
			return Output::GenericError;

		cells[position] = cell;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericEdge::InsertFace(const GenericFace* face, const unsigned int& position)
	{
		if (face == NULL || position >= faces.size())
			return Output::GenericError;

		faces[position] = face;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericEdge::InsertEdge(const GenericEdge* edge, const unsigned int& position)
	{
		if (edge == NULL || position >= edges.size())
			return Output::GenericError;

		edges[position] = edge;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericEdge::InsertPoint(const GenericPoint* point, const unsigned int& position)
	{
		if (point == NULL || position >= points.size())
			return Output::GenericError;

		points[position] = point;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericEdge::ChangeOrientation()
	{
		const GenericPoint* tempPoint = points[0];
		if(tempPoint == NULL)
			return Output::GenericError;
		points[0] = points[1];
		points[1] = tempPoint;
		return Output::Success;
	}
	// ***************************************************************************
	GenericEdge::PositionPoint GenericEdge::PointOnEdge(const Vector3d& point, const double& toll) const
	{
		Vector3d tangentVectorEdge = Point(1)->Coordinates() - Point(0)->Coordinates();
		Vector3d tangentVectorDiff = point - Point(0)->Coordinates();

		double crossProd = tangentVectorEdge.x() * tangentVectorDiff.y() - tangentVectorDiff.x() * tangentVectorEdge.y();

		double coordinateCurvilinear = (tangentVectorEdge.dot(tangentVectorDiff)) / tangentVectorEdge.squaredNorm();

		if(crossProd > toll)
			return PositionPoint::AtTheLeft;
		if(crossProd < -toll)
			return PositionPoint::AtTheRight;
		if(coordinateCurvilinear < -toll)
			return PositionPoint::Behind;
		if(coordinateCurvilinear > 1.0 + toll)
			return PositionPoint::Beyond;
		if(abs(coordinateCurvilinear) < toll)
			return PositionPoint::AtTheOrigin;
		if(abs(coordinateCurvilinear) > 1.0 - toll)
			return PositionPoint::AtTheEnd;
		return PositionPoint::Between;
	}
	// ***************************************************************************
	GenericPoint::GenericPoint(const unsigned int& _id) : GenericTreeNode(_id)
	{
		coordinates.setZero();
	}

	GenericPoint::GenericPoint(const GenericPoint& point) : GenericTreeNode(point)
	{
		coordinates = point.coordinates;
		edges.resize(point.edges.size(), NULL);
		faces.resize(point.faces.size(), NULL);
		cells.resize(point.cells.size(), NULL);
	}

	GenericPoint::~GenericPoint()
	{
		//		DestructorProperties();
		edges.clear();
		faces.clear();
		cells.clear();
	}
	// ***************************************************************************
	Output::ExitCodes GenericPoint::SetCoordinates(const Vector3d& _coordinates)
	{
		coordinates = _coordinates;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericPoint::SetCoordinates(const double& x, const double& y, const double& z)
	{
		coordinates<< x, y, z;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericPoint::AddCell(const GenericCell* cell)
	{
		if (cell == NULL)
			return Output::GenericError;

		cells.push_back(cell);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericPoint::AddFace(const GenericFace* face)
	{
		if (face == NULL)
			return Output::GenericError;

		faces.push_back(face);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericPoint::AddEdge(const GenericEdge* edge)
	{
		if (edge == NULL)
			return Output::GenericError;

		edges.push_back(edge);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericPoint::InsertCell(const GenericCell* cell, const unsigned int& position)
	{
		if (cell == NULL || position >= cells.size())
			return Output::GenericError;

		cells[position] = cell;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericPoint::InsertFace(const GenericFace* face, const unsigned int& position)
	{
		if (face == NULL || position >= faces.size())
			return Output::GenericError;

		faces[position] = face;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericPoint::InsertEdge(const GenericEdge* edge, const unsigned int& position)
	{
		if (edge == NULL || position >= edges.size())
			return Output::GenericError;

		edges[position] = edge;

		return Output::Success;
	}
	// ***************************************************************************
	GenericMesh::GenericMesh()
	{
	}

	GenericMesh::GenericMesh(const GenericMesh& mesh) : IRotation(mesh)
	{
		unsigned int numberOfCells = mesh.cells.size();
		unsigned int numberOfFaces = mesh.faces.size();
		unsigned int numberOfEdges = mesh.edges.size();
		unsigned int numberOfPoints = mesh.points.size();

		cells.resize(numberOfCells);
		faces.resize(numberOfFaces);
		edges.resize(numberOfEdges);
		points.resize(numberOfPoints);

		for(unsigned int cel = 0; cel < numberOfCells; cel++)
			cells[cel] = new GenericCell(*mesh.cells[cel]);
		for(unsigned int fac = 0; fac < numberOfFaces; fac++)
			faces[fac] = new GenericFace(*mesh.faces[fac]);
		for(unsigned int edg = 0; edg < numberOfEdges; edg++)
			edges[edg] = new GenericEdge(*mesh.edges[edg]);
		for(unsigned int pnt = 0; pnt < numberOfPoints; pnt++)
			points[pnt] = new GenericPoint(*mesh.points[pnt]);

		for(unsigned int cel = 0; cel < numberOfCells; cel++)
		{
			GenericCell& cell = *cells[cel];
			GenericCell& cellMesh = *mesh.cells[cel];
			for(unsigned int pnt = 0; pnt < cellMesh.NumberOfPoints(); pnt++)
				cell.InsertPoint(points[cellMesh.Point(pnt)->Id()], pnt);
			for(unsigned int edg = 0; edg < cellMesh.NumberOfEdges(); edg++)
				cell.InsertEdge(edges[cellMesh.Edge(edg)->Id()], edg);
			for(unsigned int fac = 0; fac < cellMesh.NumberOfFaces(); fac++)
				cell.InsertFace(faces[cellMesh.Face(fac)->Id()], fac);
			for(unsigned int cel = 0; cel < cellMesh.NumberOfCells(); cel++)
			{
				if(cellMesh.Cell(cel) != NULL)
					cell.InsertCell(cells[cellMesh.Cell(cel)->Id()], cel);
			}

			if(cellMesh.HasFather())
				cell.SetFather(cells[cellMesh.Father()->Id()]);

			if(cellMesh.HasChilds())
			{
				for(unsigned int chd = 0; chd < cellMesh.NumberOfChilds(); chd++)
					cell.InsertChild(cells[cellMesh.Child(chd)->Id()], chd);
			}
		}

		for(unsigned int fac = 0; fac < numberOfFaces; fac++)
		{
			GenericFace& face = *faces[fac];
			GenericFace& faceMesh = *mesh.faces[fac];
			for(unsigned int pnt = 0; pnt < faceMesh.NumberOfPoints(); pnt++)
				face.InsertPoint(points[faceMesh.Point(pnt)->Id()], pnt);
			for(unsigned int edg = 0; edg < faceMesh.NumberOfEdges(); edg++)
				face.InsertEdge(edges[faceMesh.Edge(edg)->Id()], edg);
			for(unsigned int fac = 0; fac < faceMesh.NumberOfFaces(); fac++)
				face.InsertFace(faces[faceMesh.Face(fac)->Id()], fac);
			for(unsigned int cel = 0; cel < faceMesh.NumberOfCells(); cel++)
			{
				if(faceMesh.Cell(cel) != NULL)
					face.InsertCell(cells[faceMesh.Cell(cel)->Id()], cel);
			}

			if(faceMesh.HasFather())
				face.SetFather(cells[faceMesh.Father()->Id()]);

			if(faceMesh.HasChilds())
			{
				for(unsigned int chd = 0; chd < faceMesh.NumberOfChilds(); chd++)
					face.InsertChild(cells[faceMesh.Child(chd)->Id()], chd);
			}
		}

		for(unsigned int edg = 0; edg < numberOfEdges; edg++)
		{
			GenericEdge& edge = *edges[edg];
			GenericEdge& edgeMesh = *mesh.edges[edg];
			for(unsigned int pnt = 0; pnt < 2; pnt++)
				edge.InsertPoint(points[edgeMesh.Point(pnt)->Id()], pnt);
			for(unsigned int edg = 0; edg < edgeMesh.NumberOfEdges(); edg++)
				edge.InsertEdge(edges[edgeMesh.Edge(edg)->Id()], edg);
			for(unsigned int fac = 0; fac < edgeMesh.NumberOfFaces(); fac++)
				edge.InsertFace(faces[edgeMesh.Face(fac)->Id()], fac);
			for(unsigned int cel = 0; cel < edgeMesh.NumberOfCells(); cel++)
			{
				if(edgeMesh.Cell(cel) != NULL)
					edge.InsertCell(cells[edgeMesh.Cell(cel)->Id()], cel);
			}

			if(edgeMesh.HasFather())
				edge.SetFather(cells[edgeMesh.Father()->Id()]);

			if(edgeMesh.HasChilds())
			{
				for(unsigned int chd = 0; chd < edgeMesh.NumberOfChilds(); chd++)
					edge.InsertChild(cells[edgeMesh.Child(chd)->Id()], chd);
			}
		}

		for(unsigned int pnt = 0; pnt < numberOfPoints; pnt++)
		{
			GenericPoint& point = *points[pnt];
			GenericPoint& pointMesh = *mesh.points[pnt];
			for(unsigned int edg = 0; edg < pointMesh.NumberOfEdges(); edg++)
				point.InsertEdge(edges[pointMesh.Edge(edg)->Id()], edg);
			for(unsigned int fac = 0; fac < pointMesh.NumberOfFaces(); fac++)
				point.InsertFace(faces[pointMesh.Face(fac)->Id()], fac);
			for(unsigned int cel = 0; cel < pointMesh.NumberOfCells(); cel++)
				point.InsertCell(cells[pointMesh.Cell(cel)->Id()], cel);
		}
	}

	GenericMesh::~GenericMesh()
	{
		for (vector<GenericPoint*>::iterator pointPtr = points.begin(); pointPtr != points.end(); pointPtr++)
			delete *pointPtr;

		for (vector<GenericEdge*>::iterator edgePtr = edges.begin(); edgePtr != edges.end(); edgePtr++)
			delete *edgePtr;

		for (vector<GenericFace*>::iterator facePtr = faces.begin(); facePtr != faces.end(); facePtr++)
			delete *facePtr;

		for (vector<GenericCell*>::iterator cellPtr = cells.begin(); cellPtr != cells.end(); cellPtr++)
			delete *cellPtr;

		points.clear();
		edges.clear();
		faces.clear();
		cells.clear();
	}
	// ***************************************************************************
	Output::ExitCodes GenericMesh::Rotate(const bool& inverse)
	{
		for (vector<GenericPoint*>::iterator pointPtr = points.begin(); pointPtr != points.end(); pointPtr++)
		{
			GenericPoint& point = **pointPtr;
			Vector3d& coordinates = point.Coordinates();

			coordinates = RotatePoint(coordinates, true, inverse);
			if(!inverse)
				coordinates[2] = 0.0;
		}

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericMesh::RotateWithInput(const Matrix3d& rotationMatrix, const Vector3d& translationVector, const bool& inverse)
	{
		for (vector<GenericPoint*>::iterator pointPtr = points.begin(); pointPtr != points.end(); pointPtr++)
		{
			GenericPoint& point = **pointPtr;
			Vector3d& coordinates = point.Coordinates();

			SetRotationMatrix(rotationMatrix);
			SetOriginTranslation(translationVector);
			coordinates = RotatePoint(coordinates, true, inverse);
			if(!inverse)
				coordinates[2] = 0.0;
		}
		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericMesh::RotateCellWithInput(const unsigned int& idCell, const Matrix3d& rotationMatrix, const Vector3d& translationVector, const bool& inverse)
	{
		GenericCell& cell = *Cell(idCell);
		for(unsigned int numPnt = 0; numPnt < cell.NumberOfPoints(); numPnt++)
		{
			GenericPoint& point = *Point(cell.Point(numPnt)->Id());
			Vector3d& coordinates = point.Coordinates();

			coordinates = RotatePointWithInput(coordinates, rotationMatrix, translationVector, inverse);
			if(!inverse)
				coordinates[2] = 0.0;
		}
		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericMesh::RotateFaceWithInput(const unsigned int& idFace, const Matrix3d& rotationMatrix, const Vector3d& translationVector, const bool& inverse)
	{
		GenericFace& face = *Face(idFace);
		for(unsigned int numPnt = 0; numPnt < face.NumberOfPoints(); numPnt++)
		{
			GenericPoint& point = *Point(face.Point(numPnt)->Id());
			Vector3d& coordinates = point.Coordinates();

			coordinates = RotatePointWithInput(coordinates, rotationMatrix, translationVector, inverse);
			if(!inverse)
				coordinates[2] = 0.0;
		}
		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericMesh::AddCell(GenericCell* cell)
	{
		if (cell == NULL)
			return Output::GenericError;

		cells.push_back(cell);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericMesh::AddFace(GenericFace* face)
	{
		if (face == NULL)
			return Output::GenericError;

		faces.push_back(face);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericMesh::AddEdge(GenericEdge* edge)
	{
		if (edge == NULL)
			return Output::GenericError;

		edges.push_back(edge);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericMesh::AddPoint(GenericPoint* point)
	{
		if (point == NULL)
			return Output::GenericError;

		points.push_back(point);

		return Output::Success;
	}
	// ***************************************************************************
	const bool GenericMesh::FindCoordinates(const Vector3d& coordinates, unsigned int& idPoint, const double& toll)
	{
		double squaredToll = toll* toll;
		for(unsigned int numPnt = 0; numPnt < points.size(); numPnt++)
		{
			const GenericPoint& point = *points[numPnt];
			Vector3d diff = coordinates - point.Coordinates();
			if(diff.squaredNorm() < squaredToll)
			{
				idPoint = point.Id();
				return true;
			}
		}
		return false;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CheckDoublePoints(const double& toll)
	{
		double squaredToll = toll*toll;
		for(unsigned int pnt = 0; pnt < points.size() - 1; pnt++)
		{
			GenericPoint& point = *points[pnt];
			for(unsigned int pnt2 = pnt+1; pnt2 < points.size(); pnt2++)
			{
				GenericPoint& point2 = *points[pnt2];
				Vector3d diff = point2.Coordinates() - point.Coordinates();
				if(diff.squaredNorm() < squaredToll)
				{
					Output::PrintWarningMessage("%s :Point %d and Point %d have the same coordinates", true, __func__, point.Id(), point2.Id());
				}
			}
		}
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CheckPointsInCells()
	{
		set<unsigned int> pointsInCells;
		vector<bool> pointsChecked(points.size(), false);
		for(unsigned int cel = 0; cel < cells.size(); cel++)
		{
			GenericCell& cell = *cells[cel];
			for(unsigned int pnt = 0; pnt < cell.NumberOfPoints(); pnt++)
			{
				if(pointsChecked[cell.Point(pnt)->Id()])
					continue;

				pointsInCells.insert(cell.Point(pnt)->Id());
				pointsChecked[cell.Point(pnt)->Id()] = true;
			}
		}
		if(pointsInCells.size() != points.size())
		{
			for(unsigned int numPnt = 0; numPnt < points.size(); numPnt++)
			{
				if(pointsChecked[numPnt])
					continue;
				Output::PrintErrorMessage("%s : There is point with id %d not used in the cells", true, __func__, numPnt);
			}
			return Output::GenericError;
		}

		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CheckPointsInFaces()
	{
		set<unsigned int> pointsInFaces;
		vector<bool> pointsChecked(points.size(), false);
		for(unsigned int fac = 0; fac < faces.size(); fac++)
		{
			GenericFace& face = *faces[fac];
			for(unsigned int pnt = 0; pnt < face.NumberOfPoints(); pnt++)
			{
				if(pointsChecked[face.Point(pnt)->Id()])
					continue;

				pointsInFaces.insert(face.Point(pnt)->Id());
				pointsChecked[face.Point(pnt)->Id()] = true;
			}
		}
		if(pointsInFaces.size() != points.size())
		{
			for(unsigned int numPnt = 0; numPnt < points.size(); numPnt++)
			{
				if(pointsChecked[numPnt])
					continue;
				Output::PrintErrorMessage("%s : There is point with id %d not used in the cells", true, __func__, numPnt);
			}
			return Output::GenericError;
		}
		return Output::Success;
	}

	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CheckPointsEqualsEdgesInCells2D()
	{
		for(unsigned int cel = 0; cel < cells.size(); cel++)
		{
			GenericCell& cell = *cells[cel];
			if(cell.NumberOfPoints() - cell.NumberOfEdges() > 0)
			{
				Output::PrintErrorMessage("%s :There are %d point/s and %d edges in the face %d", true, __func__, cell.NumberOfPoints(), cell.NumberOfEdges(), cell.Id());
				return Output::GenericError;
			}
		}
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CheckPointsEqualsEdgesInFaces()
	{
		for(unsigned int cel = 0; cel < faces.size(); cel++)
		{
			GenericFace& face = *faces[cel];
			if(face.NumberOfPoints() - face.NumberOfEdges() > 0)
			{
				Output::PrintErrorMessage("%s :There are %d point/s and %d edges in the face %d", true, __func__, face.NumberOfPoints(), face.NumberOfEdges(), face.Id());
				return Output::GenericError;
			}
		}
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CheckEdgesInFaces()
	{
		for(unsigned int fac = 0; fac < faces.size(); fac++)
		{
			GenericFace& face = *faces[fac];
			for(unsigned int edg = 0; edg < face.NumberOfEdges(); edg++)
			{
				GenericEdge& edge = *edges[face.Edge(edg)->Id()];
				bool foundFace = false;
				for(unsigned int facEdg = 0; facEdg < edge.NumberOfFaces(); facEdg++)
				{
					if(edge.Face(facEdg) == &face)
					{
						foundFace = true;
						break;
					}
				}
				if(!foundFace)
				{
					Output::PrintErrorMessage("In edge %d there is not the face %d ", edge.Id(), face.Id());
					return Output::GenericError;
				}
			}
		}
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CheckDoubleCells()
	{
		set<unsigned int> firstCell;
		for(unsigned int cel = 0; cel < cells.size()-1; cel++)
		{
			GenericCell& cell = *cells[cel];
			unsigned int numFacesOrEdges = cell.NumberOfFaces();
			bool controlEdges = false;
			if(numFacesOrEdges == 0)
			{
				controlEdges = true;
				numFacesOrEdges = cell.NumberOfEdges();
			}
			if(!controlEdges)
			{
				if(numFacesOrEdges < 4)
				{
					Output::PrintErrorMessage("%s :In cell %d there are %d faces ", true, __func__, cell.Id(), numFacesOrEdges);
					return Output::GenericError;
				}
			}
			else
			{
				if(numFacesOrEdges < 3)
				{
					Output::PrintErrorMessage("%s :In cell %d there are %d edges ", true, __func__, cell.Id(), numFacesOrEdges);
					return Output::GenericError;
				}
			}

			if(!controlEdges)
				for(unsigned int fac = 0; fac < numFacesOrEdges; fac++)
					firstCell.insert(cell.Face(fac)->Id());
			else
				for(unsigned int edg = 0; edg < numFacesOrEdges; edg++)
					firstCell.insert(cell.Edge(edg)->Id());

			if(firstCell.size() != numFacesOrEdges)
			{
				if(!controlEdges)
					Output::PrintErrorMessage("%s :In cell %d there are doubled faces", true, __func__, cell.Id());
				else
					Output::PrintErrorMessage("%s :In cell %d there are doubled edges", true, __func__, cell.Id());
				return Output::GenericError;
			}

			set<unsigned int> secondCell;
			for(unsigned int cel2 = cel+1; cel2 < cells.size(); cel2++)
			{
				GenericCell& cell2 = *cells[cel2];
				if(cell2.NumberOfFaces() != numFacesOrEdges)
					continue;

				if(!controlEdges)
					for(unsigned int fac2 = 0; fac2 < numFacesOrEdges; fac2++)
						secondCell.insert(cell2.Face(fac2)->Id());
				else
					for(unsigned int edg2 = 0; edg2 < numFacesOrEdges; edg2++)
						secondCell.insert(cell2.Edge(edg2)->Id());

				if(secondCell.size() != numFacesOrEdges)
				{
					if(!controlEdges)
						Output::PrintErrorMessage("%s :In cell %d there are doubled faces ", true, __func__, cell2.Id());
					else
						Output::PrintErrorMessage("%s :In cell %d there are doubled edges ", true, __func__, cell2.Id());
					return Output::GenericError;
				}

				set<unsigned int>::iterator it = firstCell.begin();
				set<unsigned int>::iterator it2 = secondCell.begin();
				bool cont = true;
				for(unsigned int i = 0; i < numFacesOrEdges; i++)
				{
					if(*it != *it2)
					{
						cont = false;
						break;
					}
					it++;
					it2++;
				}
				if(cont)
				{
					Output::PrintErrorMessage("%s :There is a doubled cell with id %d and %d ", true, __func__, cell.Id(), cell2.Id());
					return Output::GenericError;
				}
				secondCell.clear();
			}
			firstCell.clear();
		}
		return Output::Success;
	}

	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CheckDoubleFaces()
	{
		set<unsigned int> firstFace;
		for(unsigned int fac = 0; fac < faces.size()-1; fac++)
		{
			GenericFace& face = *faces[fac];
			unsigned int numEdges = face.NumberOfEdges();
			if(numEdges < 3)
			{
				Output::PrintErrorMessage("%s :In face %d there are %d edges ", true, __func__, face.Id(), numEdges);
				return Output::GenericError;
			}

			for(unsigned int edg = 0; edg < numEdges; edg++)
				firstFace.insert(face.Edge(edg)->Id());

			if(firstFace.size() != numEdges)
			{
				Output::PrintErrorMessage("%s :In face %d there are doubled edges", true, __func__, face.Id());
				return Output::GenericError;
			}

			set<unsigned int> secondFace;
			for(unsigned int fac2 = fac+1; fac2 < faces.size(); fac2++)
			{
				GenericFace& face2 = *faces[fac2];
				if(face2.NumberOfFaces() != numEdges)
					continue;

				for(unsigned int edg2 = 0; edg2 < numEdges; edg2++)
					secondFace.insert(face2.Edge(edg2)->Id());

				if(secondFace.size() != numEdges)
				{
					Output::PrintErrorMessage("%s :In face %d there are doubled edges ", true, __func__, face2.Id());
					return Output::GenericError;
				}

				set<unsigned int>::iterator it = firstFace.begin();
				set<unsigned int>::iterator it2 = secondFace.begin();
				bool cont = true;
				for(unsigned int i = 0; i < numEdges; i++)
				{
					if(*it != *it2)
					{
						cont = false;
						break;
					}
					it++;
					it2++;
				}
				if(cont)
				{
					Output::PrintErrorMessage("%s :There is a doubled face with id %d and %d ", true, __func__, face.Id(), face2.Id());
					return Output::GenericError;
				}
				secondFace.clear();
			}
			firstFace.clear();
		}
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CheckDoubleEdges()
	{
		set<unsigned int> firstEdge;
		for(unsigned int edg = 0; edg < edges.size()-1; edg++)
		{
			GenericEdge& edge = *edges[edg];
			unsigned int numPoints = edge.NumberOfPoints();
			if(numPoints != 2)
			{
				Output::PrintErrorMessage("%s :In edge %d there are %d points ", true, edge.Id(), __func__, numPoints);
				return Output::GenericError;
			}

			for(unsigned int edg = 0; edg < 2; edg++)
				firstEdge.insert(edge.Point(edg)->Id());

			if(firstEdge.size() != 2)
			{
				Output::PrintErrorMessage("%s :In edge %d there are doubled points", true, __func__, edge.Id());
				return Output::GenericError;
			}

			set<unsigned int> secondEdge;
			for(unsigned int edg2 = edg+1; edg2 < edges.size(); edg2++)
			{
				GenericEdge& edge2 = *edges[edg2];

				for(unsigned int edg2 = 0; edg2 < 2; edg2++)
					secondEdge.insert(edge2.Point(edg2)->Id());

				if(secondEdge.size() != numPoints)
				{
					Output::PrintErrorMessage("%s :In edge %d there are doubled point ", true, __func__, edge2.Id());
					return Output::GenericError;
				}

				set<unsigned int>::iterator it = firstEdge.begin();
				set<unsigned int>::iterator it2 = secondEdge.begin();
				bool cont = true;
				for(unsigned int i = 0; i < numPoints; i++)
				{
					if(*it != *it2)
					{
						cont = false;
						break;
					}
					it++;
					it2++;
				}
				if(cont)
				{
					Output::PrintErrorMessage("%s :There is a doubled edge with id %d and %d ", true, __func__, edge.Id(), edge2.Id());
					return Output::GenericError;
				}
				secondEdge.clear();
			}
			firstEdge.clear();
		}
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CheckNeigs()
	{
		for(unsigned int cel = 0; cel < cells.size()-1; cel++)
		{
			GenericCell& cell = *cells[cel];
			unsigned int idCell = cell.Id();
			unsigned int idCellToControl = 0;
			unsigned int dimension = Dimension();
			switch(dimension)
			{
				case 3:
				{
					for(unsigned int fac = 0; fac < cell.NumberOfFaces(); fac++)
					{
						const GenericFace& face = *cell.Face(fac);
						if(face.Cell(0) == NULL )
						{
							if(cell.Cell(fac) == NULL)
								continue;
							else
							{
								Output::PrintErrorMessage("%s :Position neigs not right", true, __func__);
								return Output::GenericError;
							}
						}
						else if(face.Cell(0) != NULL && face.Cell(0)->Id() != idCell)
						{
							idCellToControl = face.Cell(0)->Id();
							if(cell.Cell(fac)->Id() == idCellToControl)
								continue;
							else
							{
								Output::PrintErrorMessage("%s :Position neigs not right", true, __func__);
								return Output::GenericError;
							}
						}
						if(face.Cell(1) == NULL )
						{
							if(cell.Cell(fac) == NULL)
								continue;
							else
							{
								Output::PrintErrorMessage("%s :Position neigs not right", true, __func__);
								return Output::GenericError;
							}
						}
						else if(face.Cell(1) != NULL && face.Cell(1)->Id() != idCell)
						{
							idCellToControl = face.Cell(1)->Id();
							if(cell.Cell(fac)->Id() == idCellToControl)
								continue;
							else
							{
								Output::PrintErrorMessage("%s :Position neigs not right", true, __func__);
								return Output::GenericError;
							}
						}
					}

				}
					break;
				case 2:
				{
					for(unsigned int edg = 0; edg < cell.NumberOfEdges(); edg++)
					{
						const GenericEdge& edge = *cell.Edge(edg);
						if(edge.Cell(0) == NULL )
						{
							if(cell.Cell(edg) == NULL)
								continue;
							else
							{
								Output::PrintErrorMessage("%s :Position neigs not right", true, __func__);
								return Output::GenericError;
							}
						}
						else if(edge.Cell(0) != NULL && edge.Cell(0)->Id() != idCell)
						{
							idCellToControl = edge.Cell(0)->Id();
							if(cell.Cell(edg)->Id() == idCellToControl)
								continue;
							else
							{
								Output::PrintErrorMessage("%s :Position neigs not right", true, __func__);
								return Output::GenericError;
							}
						}
						if(edge.Cell(1) == NULL )
						{
							if(cell.Cell(edg) == NULL)
								continue;
							else
							{
								Output::PrintErrorMessage("%s :Position neigs not right", true, __func__);
								return Output::GenericError;
							}
						}
						else if(edge.Cell(1) != NULL && edge.Cell(1)->Id() != idCell)
						{
							idCellToControl = edge.Cell(1)->Id();
							if(cell.Cell(edg)->Id() == idCellToControl)
								continue;
							else
							{
								Output::PrintErrorMessage("%s :Position neigs not right", true, __func__);
								return Output::GenericError;
							}
						}
					}
				}
					break;
				case 1:
				{
					for(unsigned int pnt = 0; pnt < cell.NumberOfPoints(); pnt++)
					{
						const GenericPoint& point = *cell.Point(pnt);
						if(point.Cell(0) == NULL )
						{
							if(cell.Cell(pnt) == NULL)
								continue;
							else
							{
								Output::PrintErrorMessage("%s :Position neigs not right", true, __func__);
								return Output::GenericError;
							}
						}
						else if(point.Cell(0) != NULL && point.Cell(0)->Id() != idCell)
						{
							idCellToControl = point.Cell(0)->Id();
							if(cell.Cell(pnt)->Id() == idCellToControl)
								continue;
							else
							{
								Output::PrintErrorMessage("%s :Position neigs not right", true, __func__);
								return Output::GenericError;
							}
						}
						if(point.Cell(1) == NULL )
						{
							if(cell.Cell(pnt) == NULL)
								continue;
							else
							{
								Output::PrintErrorMessage("%s :Position neigs not right", true, __func__);
								return Output::GenericError;
							}
						}
						else if(point.Cell(1) != NULL && point.Cell(1)->Id() != idCell)
						{
							idCellToControl = point.Cell(1)->Id();
							if(cell.Cell(pnt)->Id() == idCellToControl)
								continue;
							else
							{
								Output::PrintErrorMessage("%s :Position neigs not right", true, __func__);
								return Output::GenericError;
							}
						}
					}
				}
					break;
				default:
					break;
			}
		}
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CutEdgeWithPoints(const unsigned int& idEdge, const vector<Vector3d>& coordinatesPoints, const bool& inheritProperty)
	{
		unsigned int numPoints = coordinatesPoints.size();
		unsigned int numEdges = numPoints +1;

		vector<GenericPoint*> pointTraces(numPoints, NULL);
		vector<GenericEdge*> newEdges(numEdges, NULL);

		GenericEdge& edge = *edges[idEdge];
		edge.SetState(false);

		GenericPoint* firstPoint = points[edge.Point(0)->Id()];
		GenericPoint* secondPoint = points[edge.Point(1)->Id()];

		unsigned int marker = edge.Marker();
		for(unsigned int pnt = 0; pnt < numPoints; pnt++)
		{
			pointTraces[pnt] = CreatePoint();
			pointTraces[pnt]->SetMarker(marker);
			AddPoint(pointTraces[pnt]);
			pointTraces[pnt]->InitializeEdges(2);
			pointTraces[pnt]->SetCoordinates(coordinatesPoints[pnt]);
		}

		unsigned int numberCell = edge.NumberOfCells();

		edge.InitializeChilds(numEdges);
		for(unsigned int edg = 0; edg < numEdges; edg++)
		{
			newEdges[edg] = CreateEdge();
			AddEdge(newEdges[edg]);
			newEdges[edg]->SetFather(&edge);
			newEdges[edg]->AllocateCells(numberCell);
			newEdges[edg]->SetMarker(marker);
			if(inheritProperty)
				newEdges[edg]->InheritPropertiesByFather();
			edge.AddChild(newEdges[edg]);
		}

		newEdges[0]->AddPoint(firstPoint);
		newEdges[0]->AddPoint(pointTraces[0]);
		newEdges[numPoints]->AddPoint(pointTraces[numPoints-1]);
		newEdges[numPoints]->AddPoint(secondPoint);

		for(unsigned int edg = 0; edg < numPoints - 1; edg++)
		{
			newEdges[edg+1]->AddPoint(pointTraces[edg]);
			newEdges[edg+1]->AddPoint(pointTraces[edg+1]);
		}

		for(unsigned int pnt = 0; pnt < numPoints; pnt++)
		{
			pointTraces[pnt]->AddEdge(newEdges[pnt]);
			pointTraces[pnt]->AddEdge(newEdges[pnt+1]);
		}

		for(unsigned int edg = 0; edg < firstPoint->NumberOfEdges(); edg++)
			if(firstPoint->Edge(edg)->Id() == edge.Id())
				firstPoint->InsertEdge(newEdges[0],edg);
		for(unsigned int edg = 0; edg < secondPoint->NumberOfEdges(); edg++)
			if(secondPoint->Edge(edg)->Id() == edge.Id())
				secondPoint->InsertEdge(newEdges[numEdges-1],edg);

		if(edge.Cell(0) != NULL)
		{
			for(unsigned int edg = 0; edg < numEdges; edg++)
				newEdges[edg]->InsertCell(edge.Cell(0),0);
			for(unsigned int pnt = 0; pnt < numPoints; pnt++)
				pointTraces[pnt]->AddCell(edge.Cell(0));
		}
		if(edge.Cell(1) != NULL)
		{
			for(unsigned int edg = 0; edg < numEdges; edg++)
				newEdges[edg]->InsertCell(edge.Cell(1),1);
			for(unsigned int pnt = 0; pnt < numPoints; pnt++)
				pointTraces[pnt]->AddCell(edge.Cell(1));
		}

		if(edge.Cell(0) != NULL)
			UpdateCell(edge.Cell(0)->Id());
		if(edge.Cell(1) != NULL)
			UpdateCell(edge.Cell(1)->Id());

		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CutEdgeWithCoordinateCurvilinears(const unsigned int& idEdge, const vector<double>& coordinatesCurvilinear, const bool& inheritProperty)
	{
		unsigned int numPoints = coordinatesCurvilinear.size();
		GenericEdge& edge = *edges[idEdge];
		const Vector3d& firstPoint = edge.Point(0)->Coordinates();
		const Vector3d& secondPoint = edge.Point(1)->Coordinates();
		vector<Vector3d> coordinatePoints(numPoints);
		Vector3d tangent = secondPoint - firstPoint;
		for(unsigned int numPnt = 0; numPnt < numPoints; numPnt++)
		{
			coordinatePoints[numPnt] = firstPoint + coordinatesCurvilinear[numPnt] * tangent;
		}

		CutEdgeWithPoints(idEdge, coordinatePoints, inheritProperty);
		return Output::Success;
	}

	// ***************************************************************************
	const Output::ExitCodes GenericMesh::UpdateCell(const unsigned int& idCell, const int& idEdge)
	{
		GenericCell& cell = *cells[idCell];
		if(!cell.IsActive())
			return Output::Success;
		unsigned int dimension = Dimension();
		switch(dimension)
		{
			case 2:
			{
				GenericCell& cellChild = *CreateCell();

				cell.InitializeChilds(1);
				cell.AddChild(&cellChild);
				cell.SetState(false);
				cellChild.SetFather(&cell);
				cellChild.InheritPropertiesByFather();

				unsigned int numberEdges = cell.NumberOfEdges(); //equal number points
				unsigned int numberEdgesChild = 0;

				vector<unsigned int> pointsIdCutEdge(2,0); //id dei punti del lato tagliato
				vector<unsigned int> newIdPoint; // nuovo punto da aggiungere

				if(idEdge == -1)
				{
					for(unsigned int numEd = 0; numEd < numberEdges; numEd++)
					{
						unsigned idEdgeCut = cell.Edge(numEd)->Id();
						GenericEdge* edge = edges[idEdgeCut];
						if(edge->HasChilds())
						{
							pointsIdCutEdge[0] = edge->Point(0)->Id();
							pointsIdCutEdge[1] = edge->Point(1)->Id();
							numberEdgesChild = edge->NumberOfChilds();
							break;
						}
					}
				}
				else
				{
					GenericEdge* edge = edges[idEdge];
					pointsIdCutEdge[0] = edge->Point(0)->Id();
					pointsIdCutEdge[1] = edge->Point(1)->Id();
					numberEdgesChild = edge->NumberOfChilds();
				}

				for(unsigned int numPoint = 0; numPoint < numberEdges; numPoint++)
				{
					unsigned int idPoint = cell.Point(numPoint)->Id();
					unsigned int idPointNext = cell.Point((numPoint+1)%numberEdges)->Id();
					if(idPoint == pointsIdCutEdge[1] && idPointNext == pointsIdCutEdge[0])
					{
						pointsIdCutEdge[0] = idPoint;
						pointsIdCutEdge[1] = idPointNext;
						break;
					}
				}

				unsigned int numEdgesCellChild = (numberEdgesChild == 0)? numberEdges : numberEdges - 1 + numberEdgesChild;
				cellChild.InitializeEdges(numEdgesCellChild);
				cellChild.InitializePoints(numEdgesCellChild);
				cellChild.AllocateCells(numEdgesCellChild);
				AddCell(&cellChild);

				newIdPoint.reserve(numEdgesCellChild-1);

				//CICLO SUI LATI DELLA CELLA
				for(unsigned int numEd = 0; numEd < numberEdges; numEd++)
				{
					unsigned int idEdge = cell.Edge(numEd)->Id();
					GenericEdge* edge = edges[idEdge];
					//SE NON HA FIGLI IL LATO:
					//1) AGGIUNGO IL LATO ALLA CELLA
					//2) AGGIORNO LA CELLA DEL LATO PADRE CON LA CELLA FIGLIO
					if(!edge->HasChilds())
					{
						cellChild.AddEdge(edge);
						for(unsigned int neigCell = 0; neigCell < 2; neigCell++)
							if(edge->Cell(neigCell) != NULL)
								if(edge->Cell(neigCell)->Id() == cell.Id())
									edge->InsertCell(&cellChild, neigCell);
					}
					//SE HA FIGLI IL LATO:
					//1) SALVO GLI ID DEI PUNTI DEL LATO PADRE CHE E' STATO TAGLIATO
					//2) AGGIUNGO ALLA CELLA I FIGLI DEL LATO
					//3) AI FIGLI DEL LATO AGGIORNO LA CELLA PADRE CON LA CELLA FIGLIO
					//4) TROVO L'ID DEL NUOVO PUNTO DA AGGIUNGERE
					//5) AGGIORNO GLI ID DEI PUNTI DEL LATO PADRE CON QUELLI DEL FIGLIO TEMPORANEAMENTE
					//   PER TROVARE L'EVENTUALE NUOVO ID DEL PUNTO SE NE SONO DUE
					//6) SALVO GLI ID DEI PUNTI DEL LATO PADRE CHE E' STATO TAGLIATO
					else
					{
						bool inversePushBack = false;
						unsigned int idEdge = edge->Child(0)->Id();
						GenericEdge* edgeChild = edges[idEdge];
						if(edgeChild->Point(0)->Id() == pointsIdCutEdge[1])
							inversePushBack = true;

						unsigned int counter = 0;
						if(inversePushBack)
						{
							for(int numChild = edge->NumberOfChilds() - 1; numChild >= 0 ; numChild--)
							{
								unsigned int idEdge = edge->Child(numChild)->Id();
								GenericEdge* edgeChild = edges[idEdge];
								cellChild.AddEdge(edgeChild);

								for(unsigned int neigCell = 0; neigCell < 2; neigCell++)
									if(edge->Cell(neigCell) != NULL)
										if(edge->Cell(neigCell)->Id() == cell.Id())
											edgeChild->InsertCell(&cellChild, neigCell);

								if(counter < edge->NumberOfChilds() - 1)
								{
									if(edgeChild->Point(0)->Id() != pointsIdCutEdge[0] && edgeChild->Point(0)->Id() != pointsIdCutEdge[1])
										newIdPoint.push_back(edgeChild->Point(0)->Id());
									else
										newIdPoint.push_back(edgeChild->Point(1)->Id());
								}
								pointsIdCutEdge[0] = edgeChild->Point(0)->Id();
								pointsIdCutEdge[1] = edgeChild->Point(1)->Id();
								counter++;
							}
							pointsIdCutEdge[0] = edge->Point(0)->Id();
							pointsIdCutEdge[1] = edge->Point(1)->Id();
						}
						else
						{
							for(unsigned int numChild =  0; numChild < edge->NumberOfChilds() ; numChild++)
							{
								unsigned int idEdge = edge->Child(numChild)->Id();
								GenericEdge* edgeChild = edges[idEdge];
								cellChild.AddEdge(edgeChild);

								for(unsigned int neigCell = 0; neigCell < 2; neigCell++)
									if(edge->Cell(neigCell) != NULL)
										if(edge->Cell(neigCell)->Id() == cell.Id())
											edgeChild->InsertCell(&cellChild, neigCell);

								if(counter < edge->NumberOfChilds() - 1)
								{
									if(edgeChild->Point(0)->Id() != pointsIdCutEdge[0] && edgeChild->Point(0)->Id() != pointsIdCutEdge[1])
										newIdPoint.push_back(edgeChild->Point(0)->Id());
									else
										newIdPoint.push_back(edgeChild->Point(1)->Id());
								}
								pointsIdCutEdge[0] = edgeChild->Point(0)->Id();
								pointsIdCutEdge[1] = edgeChild->Point(1)->Id();
								counter++;
							}
						}
						pointsIdCutEdge[0] = edge->Point(1)->Id();
						pointsIdCutEdge[1] = edge->Point(0)->Id();
					}
				}


				for(unsigned int numPoint = 0; numPoint < numberEdges; numPoint++)
				{
					unsigned int idPoint = cell.Point(numPoint)->Id();
					unsigned int idPointNext = cell.Point((numPoint+1)%numberEdges)->Id();
					GenericPoint* point = points[idPoint];
					cellChild.AddPoint(point);
					for(unsigned int cellPosition = 0; cellPosition < point->NumberOfCells(); cellPosition++)
					{
						//Tolgo il padre dai punti ed inserisco il figlio
						if(cell.Id() == point->Cell(cellPosition)->Id())
							point->InsertCell(&cellChild,cellPosition);
					}
					if((idPoint == pointsIdCutEdge[0] && idPointNext == pointsIdCutEdge[1]) || (idPoint == pointsIdCutEdge[1] && idPointNext == pointsIdCutEdge[0]))
					{
						for(unsigned int newId = 0; newId < newIdPoint.size() ; newId++)
						{
							GenericPoint* newPoint = points[newIdPoint[newId]];
							cellChild.AddPoint(newPoint);
							for(unsigned int cellPosition = 0; cellPosition < newPoint->NumberOfCells(); cellPosition++)
							{
								//Tolgo il padre dai punti ed inserisco il figlio
								if(cell.Id() == newPoint->Cell(cellPosition)->Id())
									newPoint->InsertCell(&cellChild,cellPosition);
							}
						}
					}
				}

				for(unsigned int numEd = 0; numEd < numEdgesCellChild; numEd++)
				{
					const GenericEdge& edge = *cellChild.Edge(numEd);
					for(unsigned int numNeig = 0; numNeig < edge.NumberOfCells(); numNeig++)
					{
						if(edge.Cell(numNeig) != NULL)
						{
							GenericCell& cellNeigh = *cells[edge.Cell(numNeig)->Id()];
							if(cellNeigh.Id() != cellChild.Id())
								cellChild.InsertCell(&cellNeigh, numEd);
							for(unsigned int neig = 0; neig < cellNeigh.NumberOfCells(); neig++)
								if(cellNeigh.Cell(neig) != NULL && cell.Id() == cellNeigh.Cell(neig)->Id())
									cellNeigh.InsertCell(&cellChild, neig);
						}
					}
				}

				if(cellChild.NumberOfEdges() != cellChild.NumberOfPoints())
				{
					Output::PrintErrorMessage("Failed Update Cell %d", false, cell.Id());
					return Output::GenericError;
				}
			}
				break;
			case 3:
			{
				unsigned int counterFaces = 0;
				unsigned int counterEdges = 0;
				for(unsigned int numFace = 0; numFace < cell.NumberOfFaces(); numFace++)
					if(cell.Face(numFace)->IsActive())
						counterFaces++;
				for(unsigned int numEdg = 0; numEdg < cell.NumberOfEdges(); numEdg++)
					if(cell.Edge(numEdg)->IsActive())
						counterEdges++;

				GenericCell& cellChild = *CreateCell();

				cell.InitializeChilds(1);
				cell.AddChild(&cellChild);
				cell.SetState(false);
				cellChild.SetFather(&cell);
				cellChild.InheritPropertiesByFather();

				set<unsigned int> pointsSetChild;
				set<unsigned int> edgesSetChild;
				set<unsigned int> facesSetChild;
				for(unsigned int fac = 0; fac < cell.NumberOfFaces(); fac++)
				{
					//in caso di cvx devono essere le ultime foglie dell'albero
					const GenericFace& face = *cell.Face(fac);
					if(!face.HasChilds())
						facesSetChild.insert(face.Id());
					else
					{
						list<unsigned int> idChildrenToVisit;
						idChildrenToVisit.push_back(face.Id());
						while(!idChildrenToVisit.empty())
						{
							unsigned int id = idChildrenToVisit.front();
							idChildrenToVisit.pop_front();
							GenericTreeNode* faceVisit = faces[id];
							for(unsigned int numChild = 0; numChild < faceVisit->NumberOfChilds(); numChild++)
							{
								const GenericTreeNode* child = faceVisit->Child(numChild);
								if(child->IsActive())
									facesSetChild.insert(child->Id());
								else
									idChildrenToVisit.push_back(child->Id());
							}
						}
					}
				}
				for(set<unsigned int>::iterator fc = facesSetChild.begin(); fc != facesSetChild.end(); ++fc)
				{
					GenericFace& faceReference = *faces[*fc];
					for(unsigned int pt = 0; pt < faceReference.NumberOfPoints(); pt++)
						pointsSetChild.insert(faceReference.Point(pt)->Id());

					for(unsigned int pt = 0; pt < faceReference.NumberOfEdges(); pt++)
						edgesSetChild.insert(faceReference.Edge(pt)->Id());
				}

				cellChild.AllocateCells(facesSetChild.size());
				cellChild.InitializeFaces(facesSetChild.size());
				cellChild.InitializeEdges(edgesSetChild.size());
				cellChild.InitializePoints(pointsSetChild.size());
				AddCell(&cellChild);

				for(set<unsigned int>::iterator it=pointsSetChild.begin(); it!=pointsSetChild.end(); ++it)
				{
					cellChild.AddPoint(points[*it]);
					for(unsigned int pointCell = 0; pointCell < points[*it]->NumberOfCells(); pointCell++)
					{
						if(points[*it]->Cell(pointCell)->Id() == cell.Id())
							points[*it]->InsertCell(&cellChild, pointCell);
					}
				}

				for(set<unsigned int>::iterator it = edgesSetChild.begin(); it != edgesSetChild.end(); ++it)
				{
					cellChild.AddEdge(edges[*it]);
					if(!edges[*it]->HasChilds())
					{
						for(unsigned int edgeCell = 0; edgeCell < edges[*it]->NumberOfCells(); edgeCell++)
						{
							if(edges[*it]->Cell(edgeCell)->Id() == cell.Id())
								edges[*it]->InsertCell(&cellChild,edgeCell);
						}
					}
					else
					{
						edges[*it]->AddCell(&cellChild);
					}
				}


				unsigned int positionNeigCell = 0;
				for(set<unsigned int>::iterator fc = facesSetChild.begin(); fc != facesSetChild.end(); ++fc)
				{
					GenericFace& faceReference = *faces[*fc];

					if(!faceReference.HasChilds())
					{
						cellChild.AddFace(&faceReference);
						for(unsigned int neigCell = 0; neigCell < 2; neigCell++)
						{
							if(faceReference.Cell(neigCell) != NULL)
							{
								GenericCell& cellNeigh = *cells[faceReference.Cell(neigCell)->Id()];
								if(cellNeigh.Id() == cell.Id())
									faceReference.InsertCell(&cellChild, neigCell);
								else
									cellChild.InsertCell(&cellNeigh, positionNeigCell);

								for(unsigned int neig = 0; neig < cellNeigh.NumberOfCells(); neig++)
									if(cellNeigh.Cell(neig) != NULL && cell.Id() == cellNeigh.Cell(neig)->Id())
										cellNeigh.InsertCell(&cellChild, neig);
							}
						}
						positionNeigCell++;
					}
					else
					{
						for(unsigned int fcChild = 0; fcChild < faceReference.NumberOfChilds(); fcChild++)
						{
							GenericFace& faceChildReference = *faces[faceReference.Child(fcChild)->Id()];
							cellChild.AddFace(&faceChildReference);
							//updating cells of the child faces of the face in input
							for(unsigned int neigCell = 0; neigCell < 2; neigCell++)
							{
								if(faceReference.Cell(neigCell) != NULL)
								{
									GenericCell& cellNeigh = *cells[faceReference.Cell(neigCell)->Id()];
									if(cellNeigh.Id() == cell.Id())
										faceReference.InsertCell(&cellChild, neigCell);
									else
										cellChild.InsertCell(faceReference.Cell(neigCell), positionNeigCell);

									for(unsigned int neig = 0; neig < cellNeigh.NumberOfCells(); neig++)
										if(cellNeigh.Cell(neig) != NULL && cell.Id() == cellNeigh.Cell(neig)->Id())
											cellNeigh.InsertCell(&cellChild, neig);
								}
							}
							positionNeigCell++;
						}
					}
				}
			}
				break;
			default:
				break;
		}
		return Output::Success;
	}

	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CreateCellChild2D(GenericCell& cell, GenericCell& cellFather, const list<unsigned int>& idEdgesCell, const list<unsigned int>& idPointCell, const bool& property)
	{
		//CONTROLLARE CHE AI LATI CI SIA TUTTO
		cellFather.AddChild(&cell);
		cellFather.SetState(false);

		cell.SetFather(&cellFather);
		unsigned int numEdgesCell = idEdgesCell.size();
		cell.InitializePoints(numEdgesCell);
		cell.AllocateEdges(numEdgesCell);
		cell.AllocateCells(numEdgesCell);
		if(property)
			cell.InheritPropertiesByFather();

		//CICLO SUGLI ID DEI PUNTI DA AGGIUNGERE ALLA PAGINA
		//1) Aggiungo il punto alla cella
		//2) Aggiorno la cella padre con la cella figlio ed inserisco update = true
		//3) Se Update rimane false significa che il punto è nuovo e non ha la cella padre e bisogna
		//   aggiungere direttamente il figlio
		for(list<unsigned int>::const_iterator iteratorId = idPointCell.begin(); iteratorId != idPointCell.end(); iteratorId++)
		{
			GenericPoint* point = points[*iteratorId];
			cell.AddPoint(point);
			bool update = false;
			for(unsigned int cellPosition = 0; cellPosition < point->NumberOfCells(); cellPosition++)
			{
				//Tolgo il padre dai punti ed inserisco il figlio
				if(cellFather.Id() == point->Cell(cellPosition)->Id())
				{
					update = true;
					point->InsertCell(&cell, cellPosition);
					break;
				}
			}
			if(!update)
				point->AddCell(&cell);
		}

		unsigned int numEdge = 0;
		for(list<unsigned int>::const_iterator iteratorId = idEdgesCell.begin(); iteratorId != idEdgesCell.end(); iteratorId++)
		{
			GenericEdge* edge = edges[*iteratorId];
			cell.InsertEdge(edge, numEdge);

			//CICLO SULLE CELLE VICINE AL LATO:
			//1) AGGIORNARE LA CELLA DEL LATO CON LA CELLA FIGLIO
			//2) AGGIUNGO LA CELLA VICINA ALLA CELLA CHE STO CREANDO
			for(unsigned int neigCellEdge = 0;  neigCellEdge < 2; neigCellEdge++)
			{
				if(edge->Cell(neigCellEdge) != NULL)
				{
					const unsigned int& idNeigCell = edge->Cell(neigCellEdge)->Id();
					if(idNeigCell == cellFather.Id())
						edge->InsertCell(&cell,neigCellEdge);
					else if((idNeigCell != cellFather.Id()) && (idNeigCell != cell.Id()))
					{
						cell.InsertCell(cells[idNeigCell], numEdge);
						GenericCell& cellNeigh = *cells[idNeigCell];
						for(unsigned int edg = 0; edg < cellNeigh.NumberOfEdges(); edg++)
							if(cellNeigh.Edge(edg) == edge)
								cellNeigh.InsertCell(&cell, edg);
					}
				}
			}
			numEdge++;
		}
		return Output::Success;
	}

	// ***************************************************************************
	const Output::ExitCodes GenericMesh::ComputeGeometricalProperties()
	{
		switch(Dimension())
		{
			case 1:
			{

				break;
			}
			case 2:
			{
				break;
			}
			case 3:
			{
				for(unsigned int numFac = 0; numFac < faces.size(); numFac++)
				{
					faces[numFac]->ComputeGeometricalProperties();
				}
			}
				break;
			default:
			{
				Output::PrintErrorMessage("Does not exist a mesh of dimension greater than 3", false);
				return Output::GenericError;
			}
				break;
		}
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::UpdateFace(const unsigned int& idFace, const int& idEdge)
	{
		GenericFace& face = *faces[idFace];
		unsigned int counterEdges = 0;
		for(unsigned int numEdg = 0; numEdg < face.NumberOfEdges(); numEdg++)
			if(face.Edge(numEdg)->IsActive())
				counterEdges++;

		if(counterEdges < face.NumberOfEdges() - 1)
		{
			Output::PrintWarningMessage("Update in face %d with more than two inactive edges ", true, idFace);
		}

		if(counterEdges == face.NumberOfEdges())
			return Output::Success;

		GenericFace& faceChild = *CreateFace();

		face.InitializeChilds(1);
		face.AddChild(&faceChild);
		face.SetState(false);
		faceChild.SetFather(&face);
		faceChild.InheritPropertiesByFather();

		unsigned int numberEdges = face.NumberOfEdges(); //equal number points
		unsigned int numberEdgesChild = 0;

		vector<unsigned int> pointsIdCutEdge(2,0); //id dei punti del lato tagliato
		vector<unsigned int> newIdPoint; // nuovo punto da aggiungere

		//Riempi pointsIdCutEdge e dai valore a numberEdgesChild
		if(idEdge == -1)
		{
			for(unsigned int numEd = 0; numEd < numberEdges; numEd++)
			{
				unsigned idEdgeCut = face.Edge(numEd)->Id();
				GenericEdge* edge = edges[idEdgeCut];
				if(edge->HasChilds())
				{
					pointsIdCutEdge[0] = edge->Point(0)->Id();
					pointsIdCutEdge[1] = edge->Point(1)->Id();
					numberEdgesChild = edge->NumberOfChilds();
					break;
				}
			}
		}
		else
		{
			GenericEdge* edge = edges[idEdge];
			pointsIdCutEdge[0] = edge->Point(0)->Id();
			pointsIdCutEdge[1] = edge->Point(1)->Id();
			numberEdgesChild = edge->NumberOfChilds();
		}

		for(unsigned int numPoint = 0; numPoint < numberEdges; numPoint++)
		{
			unsigned int idPoint = face.Point(numPoint)->Id();
			unsigned int idPointNext = face.Point((numPoint+1)%numberEdges)->Id();
			if(idPoint == pointsIdCutEdge[1] && idPointNext == pointsIdCutEdge[0])
			{
				pointsIdCutEdge[0] = idPoint;
				pointsIdCutEdge[1] = idPointNext;
				break;
			}
		}

		unsigned int numEdgesFaceChild = numberEdges - 1 + numberEdgesChild;
		faceChild.InitializeEdges(numEdgesFaceChild);
		faceChild.InitializePoints(numEdgesFaceChild);
		//		faceChild.AllocateFaces(numEdgesFaceChild);
		faceChild.AllocateCells(2);
		AddFace(&faceChild);

		//Vedi come inserisci figli e quali id prende

		newIdPoint.reserve(numEdgesFaceChild-1);
		//CICLO SUI LATI DELLA FACCIA
		for(unsigned int numEd = 0; numEd < numberEdges; numEd++)
		{
			unsigned int idEdge = face.Edge(numEd)->Id();
			GenericEdge* edge = edges[idEdge];
			//SE NON HA FIGLI IL LATO:
			//1) AGGIUNGO IL LATO ALLA FACCIA
			//2) AGGIORNO LA FACCIA DEL LATO PADRE CON LA FACCIA FIGLIO
			if(!edge->HasChilds())
			{
				faceChild.AddEdge(edge);
				for(unsigned int neigFace = 0; neigFace < edge->NumberOfFaces(); neigFace++)
					if(edge->Face(neigFace) != NULL)
						if(edge->Face(neigFace)->Id() == face.Id())
							edge->InsertFace(&faceChild, neigFace);
			}
			//SE HA FIGLI IL LATO:
			//1) SALVO GLI ID DEI PUNTI DEL LATO PADRE CHE E' STATO TAGLIATO
			//2) AGGIUNGO ALLA FACCIA I FIGLI DEL LATO
			//3) AI FIGLI DEL LATO AGGIORNO LA FACCIA PADRE CON LA FACCIA FIGLIO
			//4) TROVO L'ID DEL NUOVO PUNTO DA AGGIUNGERE
			//5) AGGIORNO GLI ID DEI PUNTI DEL LATO PADRE CON QUELLI DEL FIGLIO TEMPORANEAMENTE
			//   PER TROVARE L'EVENTUALE NUOVO ID DEL PUNTO SE NE SONO DUE
			//6) SALVO GLI ID DEI PUNTI DEL LATO PADRE CHE E' STATO TAGLIATO
			else
			{
				bool inversePushBack = false;
				unsigned int idEdgeChild = edge->Child(0)->Id();
				GenericEdge* edgeChild = edges[idEdgeChild];
				if(edgeChild->Point(0)->Id() == pointsIdCutEdge[1])
					inversePushBack = true;

				unsigned int counter = 0;
				if(inversePushBack)
				{
					for(int numChild = edge->NumberOfChilds() - 1; numChild >= 0 ; numChild--)
					{
						unsigned int id = edge->Child(numChild)->Id();
						GenericEdge* edgeChild = edges[id];
						faceChild.AddEdge(edgeChild);

						for(unsigned int neigFace = 0; neigFace < edge->NumberOfFaces(); neigFace++)
							if(edge->Face(neigFace) != NULL)
								if(edge->Face(neigFace)->Id() == face.Id())
									edgeChild->InsertFace(&faceChild, neigFace);

						if(counter < edge->NumberOfChilds() - 1)
						{
							if(edgeChild->Point(0)->Id() != pointsIdCutEdge[0] && edgeChild->Point(0)->Id() != pointsIdCutEdge[1])
								newIdPoint.push_back(edgeChild->Point(0)->Id());
							else
								newIdPoint.push_back(edgeChild->Point(1)->Id());
						}
						pointsIdCutEdge[0] = edgeChild->Point(0)->Id();
						pointsIdCutEdge[1] = edgeChild->Point(1)->Id();
						counter++;
					}
				}
				else
				{
					for(unsigned int numChild =  0; numChild < edge->NumberOfChilds() ; numChild++)
					{
						unsigned int id = edge->Child(numChild)->Id();
						GenericEdge* edgeChild = edges[id];
						faceChild.AddEdge(edgeChild);

						for(unsigned int neigFace = 0; neigFace < edge->NumberOfFaces(); neigFace++)
							if(edge->Face(neigFace) != NULL)
								if(edge->Face(neigFace)->Id() == face.Id())
									edgeChild->InsertFace(&faceChild, neigFace);

						if(counter < edge->NumberOfChilds() - 1)
						{
							if(edgeChild->Point(0)->Id() != pointsIdCutEdge[0] && edgeChild->Point(0)->Id() != pointsIdCutEdge[1])
								newIdPoint.push_back(edgeChild->Point(0)->Id());
							else
								newIdPoint.push_back(edgeChild->Point(1)->Id());
						}
						pointsIdCutEdge[0] = edgeChild->Point(0)->Id();
						pointsIdCutEdge[1] = edgeChild->Point(1)->Id();
						counter++;
					}
				}
				pointsIdCutEdge[0] = edge->Point(1)->Id();
				pointsIdCutEdge[1] = edge->Point(0)->Id();
			}
		}

//		pointsIdCutEdge[0] = idPoint;
//		pointsIdCutEdge[1] = idPointNext;
		for(unsigned int numPoint = 0; numPoint < numberEdges; numPoint++)
		{
			unsigned int idPoint = face.Point(numPoint)->Id();
			unsigned int idPointNext = face.Point((numPoint+1)%numberEdges)->Id();
			GenericPoint* point = points[idPoint];
			faceChild.AddPoint(point);
			for(unsigned int facePosition = 0; facePosition < point->NumberOfFaces(); facePosition++)
			{
				//Tolgo il padre dai punti ed inserisco il figlio
				if(face.Id() == point->Face(facePosition)->Id())
				{
					point->InsertFace(&faceChild, facePosition);
					break;
				}
			}
			if((idPoint == pointsIdCutEdge[0] && idPointNext == pointsIdCutEdge[1]) || (idPoint == pointsIdCutEdge[1] && idPointNext == pointsIdCutEdge[0]))
			{
				for(unsigned int newId = 0; newId < newIdPoint.size() ; newId++)
				{
					GenericPoint* newPoint = points[newIdPoint[newId]];
					faceChild.AddPoint(newPoint);
					for(unsigned int facePosition = 0; facePosition < newPoint->NumberOfFaces(); facePosition++)
					{
						//Tolgo il padre dai punti ed inserisco il figlio
						if(face.Id() == newPoint->Face(facePosition)->Id())
							newPoint->InsertFace(&faceChild,facePosition);
					}
				}
			}
		}

		//		for(unsigned int numEd = 0; numEd < numEdgesFaceChild; numEd++)
		//		{
		//			const GenericEdge& edge = *faceChild.Edge(numEd);
		//			for(unsigned int numNeig = 0; numNeig < edge.NumberOfFaces(); numNeig++)
		//			{
		//				if(edge.Face(numNeig) != NULL)
		//				{
		//					GenericFace& faceNeigh = *faces[edge.Face(numNeig)->Id()];
		//					if(faceNeigh.Id() != faceChild.Id())
		//						faceChild.InsertFace(&faceNeigh, numEd);
		//					for(unsigned int neig = 0; neig < faceNeigh.NumberOfFaces(); neig++)
		//						if(faceNeigh.Face(neig) != NULL && face.Id() == faceNeigh.Face(neig)->Id())
		//							faceNeigh.InsertFace(&faceChild, neig);
		//				}
		//			}
		//		}

		for(unsigned int numCell = 0; numCell < 2; numCell++)
		{
			if(face.Cell(numCell) != NULL)
			{
				faceChild.InsertCell(face.Cell(numCell), numCell);
			}
		}

		faceChild.SetCentroid(face.Centroid());
		faceChild.SetMeasure(face.Measure());
		faceChild.SetRotationMatrix(face.RotationMatrix());
		faceChild.SetNormal(face.Normal());

//		unsigned int numEdgesNewFace = faceChild.NumberOfEdges();
//		const Vector3d& normalFace = face.Normal();
//		faceChild.ComputeNormal();
//		Vector3d& normalNewFace = faceChild.Normal();
//		//the orientation is negative, therfore the vector of the points of the face is reversed
//		if(normalFace.dot(normalNewFace) < 1.0E-7)
//		{
//			unsigned int forIteration =  ((numEdgesNewFace - 1) * 0.5);
//			for(unsigned int i = 0; i < forIteration; i++)
//			{
//				const GenericPoint* tmpPoint = faceChild.Point(i + 1);
//				faceChild.InsertPoint(faceChild.Point(numEdgesNewFace - i - 1), i + 1);
//				faceChild.InsertPoint(tmpPoint, numEdgesNewFace - i - 1);
//			}
//			normalNewFace *= -1.0;
//		}

		//			for(unsigned int k = 0; k < numEdgesNewFace; k++)
		//			{
		//				const GenericPoint* tmpPoint = faceChild.Point(k);
		//				const GenericPoint* tmpPointNext = faceChild.Point((k + 1)%numEdgesNewFace);
		//				const GenericEdge* tempEdge = faceChild.Edge(k);
		//				const GenericPoint* tmpPointEd = tempEdge->Point(0);
		//				const GenericPoint* tmpPointEdNext = tempEdge->Point(1);
		//				if((tmpPoint == tmpPointEd && tmpPointNext == tmpPointEdNext) ||
		//					 (tmpPoint == tmpPointEdNext && tmpPointNext == tmpPointEd))
		//					continue;
		//				else
		//				{
		//					for(unsigned int ed = k; ed < numEdgesNewFace; ed++)
		//					{
		//						const GenericEdge* tempEdge2 = faceChild.Edge(ed);
		//						const GenericPoint* tmpPointEd = tempEdge2->Point(0);
		//						const GenericPoint* tmpPointEdNext = tempEdge2->Point(1);
		//						if((tmpPoint == tmpPointEd && tmpPointNext == tmpPointEdNext) ||
		//							 (tmpPoint == tmpPointEdNext && tmpPointNext == tmpPointEd))
		//						{
		//							faceChild.InsertEdge(tempEdge2, k);
		//							faceChild.InsertEdge(tempEdge, ed);
		//						}
		//						else
		//							continue;

//		faceChild.ComputeGeometricalProperties();
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CreateFaceChild(GenericFace& face, GenericFace& faceFather, const list<unsigned int>& idEdgesFace, const list<unsigned int>& idPointFace, const bool& property)
	{
		//CONTROLLARE CHE AI LATI CI SIA TUTTO
		faceFather.AddChild(&face);
		faceFather.SetState(false);

		face.SetFather(&faceFather);
		unsigned int numEdgesFace = idEdgesFace.size();
		face.InitializePoints(numEdgesFace);
		face.AllocateCells(2);
		face.AllocateEdges(numEdgesFace);
		if(property)
			face.InheritPropertiesByFather();

		//CICLO SUGLI ID DEI PUNTI DA AGGIUNGERE ALLA PAGINA
		//1) Aggiungo il punto alla faccia
		//2) Aggiorno la faccia padre con la faccia figlio ed inserisco update = true
		//3) Se Update rimane false significa che il punto è nuovo e non ha la faccia padre e bisogna
		//   aggiungere direttamente il figlio
		for(list<unsigned int>::const_iterator iteratorId = idPointFace.begin(); iteratorId != idPointFace.end(); iteratorId++)
		{
			GenericPoint* point = points[*iteratorId];
			face.AddPoint(point);
			bool update = false;
			for(unsigned int facePosition = 0; facePosition < point->NumberOfFaces(); facePosition++)
			{
				//Tolgo il padre dai punti ed inserisco il figlio
				if(faceFather.Id() == point->Face(facePosition)->Id())
				{
					update = true;
					point->InsertFace(&face, facePosition);
					break;
				}
			}
			if(!update)
				point->AddFace(&face);
		}

		unsigned int numEdge = 0;
		for(list<unsigned int>::const_iterator iteratorId = idEdgesFace.begin(); iteratorId != idEdgesFace.end(); iteratorId++)
		{
			GenericEdge* edge = edges[*iteratorId];
			face.InsertEdge(edge, numEdge++);

			//CICLO SULLE FACCE VICINE AL LATO:
			//1) AGGIORNARE LA FACCIA DEL LATO CON LA FACCUA FIGLIO
			for(unsigned int neigFaceEdge = 0;  neigFaceEdge < edge->NumberOfFaces(); neigFaceEdge++)
			{
				if(edge->Face(neigFaceEdge) != NULL)
					if(edge->Face(neigFaceEdge)->Id() == faceFather.Id())
						edge->InsertFace(&face, neigFaceEdge);
			}
		}

		if(faceFather.Cell(0) != NULL)
			face.InsertCell(faceFather.Cell(0), 0);

		if(faceFather.Cell(1) != NULL)
			face.InsertCell(faceFather.Cell(1), 1);

		face.ComputeGeometricalProperties();
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::ActivateFatherNodes()
	{
		for(unsigned int numEdge = 0; numEdge < edges.size(); numEdge++)
		{
			GenericEdge& edge = *edges[numEdge];
			if(edge.HasFather())
			{
				edge.SetState(false);
				const unsigned int idFather = edge.Father()->Id();
				GenericEdge& edgeFather = *edges[idFather];
				edgeFather.SetState(true);
			}
		}
		for(unsigned int numFaces = 0; numFaces < faces.size(); numFaces++)
		{
			GenericFace& face = *faces[numFaces];
			if(face.HasFather())
			{
				face.SetState(false);
				const unsigned int idFather = face.Father()->Id();
				GenericFace& faceFather = *faces[idFather];
				faceFather.SetState(true);
			}
		}
		for(unsigned int numCell = 0; numCell < cells.size(); numCell++)
		{
			GenericCell& cell = *cells[numCell];
			if(cell.HasFather())
			{
				cell.SetState(false);
				const unsigned int idFather = cell.Father()->Id();
				GenericCell& cellFather = *cells[idFather];
				cellFather.SetState(true);
			}
		}
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::ActivateChildrenNodes()
	{
		for(unsigned int numEdge = 0; numEdge < edges.size(); numEdge++)
		{
			GenericEdge& edge = *edges[numEdge];
			if(edge.HasChilds())
			{
				edge.SetState(false);
				for(unsigned int numChild = 0; numChild < edge.NumberOfChilds(); numChild++)
				{
					const unsigned int idFather = edge.Child(numChild)->Id();
					GenericEdge& edgeChild = *edges[idFather];
					edgeChild.SetState(true);
				}
			}
		}
		for(unsigned int numFaces = 0; numFaces < faces.size(); numFaces++)
		{
			GenericFace& face = *faces[numFaces];
			if(face.HasChilds())
			{
				face.SetState(false);
				for(unsigned int numChild = 0; numChild < face.NumberOfChilds(); numChild++)
				{
					const unsigned int idFather = face.Child(numChild)->Id();
					GenericFace& faceChild = *faces[idFather];
					faceChild.SetState(true);
				}
			}
		}
		for(unsigned int numCell = 0; numCell < cells.size(); numCell++)
		{
			GenericCell& cell = *cells[numCell];
			if(cell.HasChilds())
			{
				cell.SetState(false);
				for(unsigned int numChild = 0; numChild < cell.NumberOfChilds(); numChild++)
				{
					const unsigned int idFather = cell.Child(numChild)->Id();
					GenericCell& cellChild = *cells[idFather];
					cellChild.SetState(true);
				}
			}
		}
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CleanInactiveTreeNode()
	{
		vector<GenericPoint*> pointsTemp;
		vector<GenericEdge*> edgesTemp;
		vector<GenericFace*> facesTemp;
		vector<GenericCell*> cellsTemp;

		pointsTemp.reserve(NumberOfPoints());
		edgesTemp.reserve(NumberOfEdges());
		facesTemp.reserve(NumberOfFaces());
		cellsTemp.reserve(NumberOfCells());

		unsigned int pointId = 0;
		unsigned int edgeId = 0;
		unsigned int faceId = 0;
		unsigned int cellId = 0;

		for (vector<GenericPoint*>::iterator pointPtr = points.begin(); pointPtr != points.end(); pointPtr++)
		{
			GenericPoint& point = **pointPtr;
			if(point.IsActive())
			{
				point.SetId(pointId++);
				pointsTemp.push_back(*pointPtr);
				for(int numEdge = (int)point.NumberOfEdges() - 1; numEdge >= 0; numEdge--)
					if(!point.Edge(numEdge)->IsActive())
						point.EraseEdge(numEdge);
				for(int numFace = (int)point.NumberOfFaces() - 1; numFace >= 0; numFace--)
					if(!point.Face(numFace)->IsActive())
						point.EraseFace(numFace);
				for(int numCell = (int)point.NumberOfCells() - 1; numCell >= 0 ; numCell--)
					if(!point.Cell(numCell)->IsActive())
						point.EraseCell(numCell);
			}
			else
				delete &point;
		}

		for (vector<GenericEdge*>::iterator edgePtr = edges.begin(); edgePtr != edges.end(); edgePtr++)
		{
			GenericEdge& edge = **edgePtr;
			if(edge.IsActive())
			{
				edge.SetId(edgeId++);
				edgesTemp.push_back(&edge);

				for(int numPnt = (int)edge.NumberOfPoints() - 1; numPnt >= 0; numPnt--)
					if(!edge.Point(numPnt)->IsActive())
						edge.ErasePoint(numPnt);
				for(int numEdge = (int)edge.NumberOfEdges() - 1; numEdge >= 0; numEdge--)
					if(!edge.Edge(numEdge)->IsActive())
						edge.EraseEdge(numEdge);
				for(int numFace = (int)edge.NumberOfFaces() - 1; numFace >= 0; numFace--)
					if(!edge.Face(numFace)->IsActive())
						edge.EraseFace(numFace);
				for(int numCell = (int)edge.NumberOfCells() - 1; numCell >= 0 ; numCell--)
					if((edge.Cell(numCell) != NULL) && !edge.Cell(numCell)->IsActive())
						edge.EraseCell(numCell);

				edge.SetFather(NULL);
			}
			else
				delete &edge;
		}

		for (vector<GenericFace*>::iterator facePtr = faces.begin(); facePtr != faces.end(); facePtr++)
		{
			GenericFace& face = **facePtr;
			if(face.IsActive())
			{
				face.SetId(faceId++);
				facesTemp.push_back(*facePtr);

				for(int numPnt = (int)face.NumberOfPoints() - 1; numPnt >= 0; numPnt--)
					if(!face.Point(numPnt)->IsActive())
						face.ErasePoint(numPnt);
				for(int numEdge = (int)face.NumberOfEdges() - 1; numEdge >= 0; numEdge--)
					if(!face.Edge(numEdge)->IsActive())
						face.EraseEdge(numEdge);
				for(int numFace = (int)face.NumberOfFaces() - 1; numFace >= 0; numFace--)
					if(!face.Face(numFace)->IsActive())
						face.EraseFace(numFace);
				for(int numCell = (int)face.NumberOfCells() - 1; numCell >= 0 ; numCell--)
					if((face.Cell(numCell) != NULL) && !face.Cell(numCell)->IsActive())
						face.EraseCell(numCell);

				face.SetFather(NULL);
			}
			else
				delete &face;
		}

		for (vector<GenericCell*>::iterator cellPtr = cells.begin(); cellPtr != cells.end(); cellPtr++)
		{
			GenericCell& cell = **cellPtr;
			if(cell.IsActive())
			{
				cell.SetId(cellId++);
				cellsTemp.push_back(*cellPtr);

				for(int numPnt = (int)cell.NumberOfPoints() - 1; numPnt >= 0; numPnt--)
					if(!cell.Point(numPnt)->IsActive())
						cell.ErasePoint(numPnt);
				for(int numEdge = (int)cell.NumberOfEdges() - 1; numEdge >= 0; numEdge--)
					if(!cell.Edge(numEdge)->IsActive())
						cell.EraseEdge(numEdge);
				for(int numFace = (int)cell.NumberOfFaces() - 1; numFace >= 0; numFace--)
					if(!cell.Face(numFace)->IsActive())
						cell.EraseFace(numFace);
				for(int numCell = (int)cell.NumberOfCells() - 1; numCell >= 0 ; numCell--)
					if((cell.Cell(numCell) != NULL) && !cell.Cell(numCell)->IsActive())
						cell.EraseCell(numCell);

				cell.SetFather(NULL);
			}
			else
				delete &cell;
		}

		unsigned int numPoints = pointsTemp.size();
		unsigned int numEdges = edgesTemp.size();
		unsigned int numFaces = facesTemp.size();
		unsigned int numCells = cellsTemp.size();

		points.resize(numPoints, NULL);
		edges.resize(numEdges, NULL);
		faces.resize(numFaces, NULL);
		cells.resize(numCells, NULL);

		for(unsigned int pnt = 0 ; pnt < numPoints; pnt++)
			points[pnt] = pointsTemp[pnt];

		for(unsigned int edg = 0 ; edg < numEdges; edg++)
			edges[edg] = edgesTemp[edg];

		for(unsigned int fac = 0 ; fac < numFaces; fac++)
			faces[fac] = facesTemp[fac];

		for(unsigned int cel = 0 ; cel < numCells; cel++)
			cells[cel] = cellsTemp[cel];

		return Output::Success;
	}
	// ***************************************************************************
	GenericMeshImportInterface::GenericMeshImportInterface()
	{
		maximumCellSize = 0.0;
		minimumNumberOfCells = 0;
		markerDimension = 1;
		vertexMarkers.resize(1);
		edgeMarkers.resize(1);
		faceMarkers.resize(1);
	}
	GenericMeshImportInterface::~GenericMeshImportInterface()
	{
	}
	// ***************************************************************************
	const vector<unsigned int> GenericMeshImportInterface::MarkersPerVertex(const unsigned int& position) const
	{
		vector<unsigned int> temp(markerDimension, 0);
		for(unsigned int numMark = 0; numMark < markerDimension; numMark++)
			temp[numMark] = vertexMarkers[numMark][position];

		return temp;
	}
	// ***************************************************************************
	const vector<unsigned int> GenericMeshImportInterface::MarkersPerEdge(const unsigned int& position) const
	{
		vector<unsigned int> temp(markerDimension, 0);
		for(unsigned int numMark = 0; numMark < markerDimension; numMark++)
			temp[numMark] = edgeMarkers[numMark][position];

		return temp;
	}
	// ***************************************************************************
	const vector<unsigned int> GenericMeshImportInterface::MarkersPerFace(const unsigned int& position) const
	{
		vector<unsigned int> temp(markerDimension, 0);
		for(unsigned int numMark = 0; numMark < markerDimension; numMark++)
			temp[numMark] = faceMarkers[numMark][position];

		return temp;
	}
	// ***************************************************************************
	void GenericMeshImportInterface::SetBoundaryConditions(const vector<unsigned int>& _vertexMarkers, const unsigned int& position)
	{
		if (_vertexMarkers.size() > 0)
		{
			vertexMarkers[position].resize(_vertexMarkers.size());
			memcpy(&vertexMarkers[position][0], &_vertexMarkers[0], _vertexMarkers.size() * sizeof(unsigned int));
		}

		edgeMarkers.clear();
		faceMarkers.clear();

		if(position == 0)
			CreateUniDimMarkers();
	}
	// ***************************************************************************
	void GenericMeshImportInterface::SetBoundaryConditions(const vector<unsigned int>& _vertexMarkers, const vector<unsigned int>& _edgeMarkers, const unsigned int& position)
	{
		if (_vertexMarkers.size() > 0)
		{
			vertexMarkers[position].resize(_vertexMarkers.size());
			memcpy(&vertexMarkers[position][0], &_vertexMarkers[0], _vertexMarkers.size() * sizeof(unsigned int));
		}

		if (_edgeMarkers.size() > 0)
		{
			edgeMarkers[position].resize(_edgeMarkers.size());
			memcpy(&edgeMarkers[position][0], &_edgeMarkers[0], _edgeMarkers.size() * sizeof(unsigned int));
		}
		faceMarkers.clear();

		if(position == 0)
			CreateUniDimMarkers();
	}
	// ***************************************************************************
	void GenericMeshImportInterface::SetBoundaryConditions(const vector<unsigned int>& _vertexMarkers, const vector<unsigned int>& _edgeMarkers, const vector<unsigned int>& _faceMarkers, const unsigned int& position)
	{
		if (_vertexMarkers.size() > 0)
		{
			vertexMarkers[position].resize(_vertexMarkers.size());
			memcpy(&vertexMarkers[position][0], &_vertexMarkers[0], _vertexMarkers.size() * sizeof(unsigned int));
		}

		if (_edgeMarkers.size() > 0)
		{
			edgeMarkers[position].resize(_edgeMarkers.size());
			memcpy(&edgeMarkers[position][0], &_edgeMarkers[0], _edgeMarkers.size() * sizeof(unsigned int));
		}

		if (_faceMarkers.size() > 0)
		{
			faceMarkers[position].resize(_faceMarkers.size());
			memcpy(&faceMarkers[position][0], &_faceMarkers[0], _faceMarkers.size() * sizeof(unsigned int));
		}

		if(position == 0)
			CreateUniDimMarkers();
	}

	// ***************************************************************************
	const Output::ExitCodes GenericMeshImportInterface::CreateUniDimMarkers()
	{
		unsigned int sizeVertexMarker = vertexMarkers[0].size();
		unsigned int sizeEdgeMarker = edgeMarkers[0].size();
		unsigned int sizeFaceMarker = faceMarkers[0].size();

		unidimensionalVertexMarkers.resize(sizeVertexMarker);
		for(unsigned int val = 0; val < sizeVertexMarker; val++)
			unidimensionalVertexMarkers[val] = val + 1;

		if (sizeEdgeMarker > 0)
		{
			unidimensionalEdgeMarkers.resize(sizeEdgeMarker);
			for(unsigned int val = 0; val < sizeEdgeMarker; val++)
				unidimensionalEdgeMarkers[val] = val + 1 + sizeVertexMarker;

		}

		if (sizeFaceMarker > 0)
		{
			unidimensionalFaceMarkers.resize(sizeFaceMarker);
			for(unsigned int val = 0; val < sizeFaceMarker; val++)
				unidimensionalFaceMarkers[val] = val + 1 + sizeVertexMarker + sizeFaceMarker;
		}
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMeshImportInterface::CreateMappedMarkers()
	{
		unsigned int sizeVertexMarker = vertexMarkers[0].size();
		unsigned int sizeEdgeMarker = edgeMarkers[0].size();
		unsigned int sizeFaceMarker = faceMarkers[0].size();

		unsigned int totalSizeVertex = sizeVertexMarker + sizeEdgeMarker + sizeFaceMarker;
		mappedMarkers.resize(totalSizeVertex);
		for(unsigned int val = 0; val < sizeVertexMarker; val++)
			mappedMarkers[val] = MarkersPerVertex(val);

		if (sizeEdgeMarker > 0)
		{
			for(unsigned int val = 0; val < sizeEdgeMarker; val++)
				mappedMarkers[val + sizeVertexMarker] = MarkersPerEdge(val);
		}

		if (sizeFaceMarker > 0)
		{
			for(unsigned int val = 0; val < sizeFaceMarker; val++)
				mappedMarkers[val + sizeVertexMarker + sizeFaceMarker] = MarkersPerFace(val);
		}

		return Output::Success;
	}

	// ***************************************************************************
	InterfaceGenericMesh::InterfaceGenericMesh()
	{

	}

	InterfaceGenericMesh::~InterfaceGenericMesh()
	{
		for(unsigned int numLinkCell = 0; numLinkCell < linkedCells.size(); numLinkCell++)
			linkedCells[numLinkCell].clear();

		for(unsigned int numLinkFace = 0; numLinkFace < linkedFaces.size(); numLinkFace++)
			linkedFaces[numLinkFace].clear();

		for(unsigned int numLinkEdge = 0; numLinkEdge < linkedEdges.size(); numLinkEdge++)
			linkedEdges[numLinkEdge].clear();

		for(unsigned int numLinkPoint = 0; numLinkPoint < linkedPoints.size(); numLinkPoint++)
			linkedPoints[numLinkPoint].clear();

		linkedCells.clear();
		linkedFaces.clear();
		linkedEdges.clear();
		linkedPoints.clear();
	}

	Output::ExitCodes InterfaceGenericMesh::InsertLinkedCell(const vector<const GenericCell*>& linkCells, const unsigned int& position)
	{
		if (linkCells.size() == 0 || position >= linkedCells.size())
			return Output::GenericError;

		linkedCells[position] = linkCells;

		return Output::Success;
	}

	// ***************************************************************************
	Output::ExitCodes InterfaceGenericMesh::InsertLinkedFace(const vector<const GenericFace*>& linkFaces, const unsigned int& position)
	{
		if (linkFaces.size() == 0 || position >= linkedFaces.size())
			return Output::GenericError;

		linkedFaces[position] = linkFaces;

		return Output::Success;
	}

	// ***************************************************************************
	Output::ExitCodes InterfaceGenericMesh::InsertLinkedEdge(const vector<const GenericEdge*>& linkEdges, const unsigned int& position)
	{
		if (linkEdges.size() == 0 || position >= linkedEdges.size())
			return Output::GenericError;

		linkedEdges[position] = linkEdges;

		return Output::Success;
	}

	// ***************************************************************************
	Output::ExitCodes InterfaceGenericMesh::InsertLinkedPoint(const vector<const GenericPoint*>& linkPoints, const unsigned int& position)
	{
		if (linkPoints.size() == 0 || position >= linkedPoints.size())
			return Output::GenericError;

		linkedPoints[position] = linkPoints;

		return Output::Success;
	}

	// ***************************************************************************
	Output::ExitCodes InterfaceGenericMesh::AddLinkedCell(const vector<const GenericCell*>& linkCells)
	{
		if (linkCells.size() == 0)
			return Output::GenericError;

		linkedCells.push_back(linkCells);

		return Output::Success;
	}

	// ***************************************************************************
	Output::ExitCodes InterfaceGenericMesh::AddLinkedFace(const vector<const GenericFace*>& linkFaces)
	{
		if (linkFaces.size() == 0)
			return Output::GenericError;

		linkedFaces.push_back(linkFaces);

		return Output::Success;
	}

	// ***************************************************************************
	Output::ExitCodes InterfaceGenericMesh::AddLinkedEdge(const vector<const GenericEdge*>& linkEdges)
	{
		if (linkEdges.size() == 0)
			return Output::GenericError;

		linkedEdges.push_back(linkEdges);

		return Output::Success;
	}

	// ***************************************************************************
	Output::ExitCodes InterfaceGenericMesh::AddLinkedPoint(const vector<const GenericPoint*>& linkPoints)
	{
		if (linkPoints.size() == 0)
			return Output::GenericError;

		linkedPoints.push_back(linkPoints);

		return Output::Success;
	}

	// ***************************************************************************
	const Output::ExitCodes InterfaceGenericMesh::CleanInactiveTreeNode()
	{
		vector<GenericPoint*> pointsTemp;
		vector<GenericEdge*> edgesTemp;
		vector<GenericFace*> facesTemp;
		vector<GenericCell*> cellsTemp;

		vector< vector<const GenericCell*> > linkedCellsTemp;
		vector< vector<const GenericFace*> > linkedFacesTemp;
		vector< vector<const GenericEdge*> > linkedEdgesTemp;
		vector< vector<const GenericPoint*> > linkedPointsTemp;

		pointsTemp.reserve(NumberOfPoints());
		edgesTemp.reserve(NumberOfEdges());
		facesTemp.reserve(NumberOfFaces());
		cellsTemp.reserve(NumberOfCells());

		if(NumberOfLinkedCells() > 0)
			linkedCellsTemp.reserve(NumberOfCells());
		if(NumberOfLinkedFaces() > 0)
			linkedFacesTemp.reserve(NumberOfCells());
		if(NumberOfLinkedEdges() > 0)
			linkedEdgesTemp.reserve(NumberOfCells());
		if(NumberOfLinkedPoints() > 0)
			linkedPointsTemp.reserve(NumberOfPoints());

		unsigned int pointId = 0;
		unsigned int edgeId = 0;
		unsigned int faceId = 0;
		unsigned int cellId = 0;

		for (vector<GenericPoint*>::iterator pointPtr = points.begin(); pointPtr != points.end(); pointPtr++)
		{
			GenericPoint& point = **pointPtr;
			if(point.IsActive())
			{
				if(linkedPointsTemp.capacity() > 0)
					linkedPointsTemp.push_back(linkedPoints[point.Id()]);
				point.SetId(pointId++);
				pointsTemp.push_back(*pointPtr);
				for(int numEdge = (int)point.NumberOfEdges() - 1; numEdge >= 0; numEdge--)
					if(!point.Edge(numEdge)->IsActive())
						point.EraseEdge(numEdge);
				for(int numFace = (int)point.NumberOfFaces() - 1; numFace >= 0; numFace--)
					if(!point.Face(numFace)->IsActive())
						point.EraseFace(numFace);
				for(int numCell = (int)point.NumberOfCells() - 1; numCell >= 0 ; numCell--)
					if(!point.Cell(numCell)->IsActive())
						point.EraseCell(numCell);
			}
			else
				delete &point;
		}

		for (vector<GenericEdge*>::iterator edgePtr = edges.begin(); edgePtr != edges.end(); edgePtr++)
		{
			GenericEdge& edge = **edgePtr;
			if(edge.IsActive())
			{
				edge.SetId(edgeId++);
				edgesTemp.push_back(&edge);

				for(int numPnt = (int)edge.NumberOfPoints() - 1; numPnt >= 0; numPnt--)
					if(!edge.Point(numPnt)->IsActive())
						edge.ErasePoint(numPnt);
				for(int numEdge = (int)edge.NumberOfEdges() - 1; numEdge >= 0; numEdge--)
					if(!edge.Edge(numEdge)->IsActive())
						edge.EraseEdge(numEdge);
				for(int numFace = (int)edge.NumberOfFaces() - 1; numFace >= 0; numFace--)
					if(!edge.Face(numFace)->IsActive())
						edge.EraseFace(numFace);
				for(int numCell = (int)edge.NumberOfCells() - 1; numCell >= 0 ; numCell--)
					if((edge.Cell(numCell) != NULL) && !edge.Cell(numCell)->IsActive())
						edge.EraseCell(numCell);

				edge.SetFather(NULL);
			}
			else
				delete &edge;
		}

		for (vector<GenericFace*>::iterator facePtr = faces.begin(); facePtr != faces.end(); facePtr++)
		{
			GenericFace& face = **facePtr;
			if(face.IsActive())
			{
				face.SetId(faceId++);
				facesTemp.push_back(*facePtr);

				for(int numPnt = (int)face.NumberOfPoints() - 1; numPnt >= 0; numPnt--)
					if(!face.Point(numPnt)->IsActive())
						face.ErasePoint(numPnt);
				for(int numEdge = (int)face.NumberOfEdges() - 1; numEdge >= 0; numEdge--)
					if(!face.Edge(numEdge)->IsActive())
						face.EraseEdge(numEdge);
				for(int numFace = (int)face.NumberOfFaces() - 1; numFace >= 0; numFace--)
					if(!face.Face(numFace)->IsActive())
						face.EraseFace(numFace);
				for(int numCell = (int)face.NumberOfCells() - 1; numCell >= 0 ; numCell--)
					if((face.Cell(numCell) != NULL) && !face.Cell(numCell)->IsActive())
						face.EraseCell(numCell);

				face.SetFather(NULL);
			}
			else
				delete &face;
		}

		for (vector<GenericCell*>::iterator cellPtr = cells.begin(); cellPtr != cells.end(); cellPtr++)
		{
			GenericCell& cell = **cellPtr;
			if(cell.IsActive())
			{
				if(linkedCellsTemp.capacity() > 0)
					linkedCellsTemp.push_back(linkedCells[cell.Id()]);
				if(linkedFacesTemp.capacity() > 0)
					linkedFacesTemp.push_back(linkedFaces[cell.Id()]);
				if(linkedEdgesTemp.capacity() > 0)
					linkedEdgesTemp.push_back(linkedEdges[cell.Id()]);

				cell.SetId(cellId++);
				cellsTemp.push_back(*cellPtr);

				for(int numPnt = (int)cell.NumberOfPoints() - 1; numPnt >= 0; numPnt--)
					if(!cell.Point(numPnt)->IsActive())
						cell.ErasePoint(numPnt);
				for(int numEdge = (int)cell.NumberOfEdges() - 1; numEdge >= 0; numEdge--)
					if(!cell.Edge(numEdge)->IsActive())
						cell.EraseEdge(numEdge);
				for(int numFace = (int)cell.NumberOfFaces() - 1; numFace >= 0; numFace--)
					if(!cell.Face(numFace)->IsActive())
						cell.EraseFace(numFace);
				for(int numCell = (int)cell.NumberOfCells() - 1; numCell >= 0 ; numCell--)
					if((cell.Cell(numCell) != NULL) && !cell.Cell(numCell)->IsActive())
						cell.EraseCell(numCell);

				cell.SetFather(NULL);
			}
			else
				delete &cell;
		}

		unsigned int numPoints = pointsTemp.size();
		unsigned int numEdges = edgesTemp.size();
		unsigned int numFaces = facesTemp.size();
		unsigned int numCells = cellsTemp.size();

		points.resize(numPoints, NULL);
		edges.resize(numEdges, NULL);
		faces.resize(numFaces, NULL);
		cells.resize(numCells, NULL);

		AllocateLinkedCells(linkedCellsTemp.size());
		AllocateLinkedFaces(linkedFacesTemp.size());
		AllocateLinkedEdges(linkedEdgesTemp.size());
		AllocateLinkedPoints(linkedPointsTemp.size());

		for(unsigned int pnt = 0 ; pnt < numPoints; pnt++)
		{
			points[pnt] = pointsTemp[pnt];
			InsertLinkedPoint(linkedPointsTemp[pnt], pnt);
		}

		for(unsigned int edg = 0 ; edg < numEdges; edg++)
			edges[edg] = edgesTemp[edg];

		for(unsigned int fac = 0 ; fac < numFaces; fac++)
			faces[fac] = facesTemp[fac];

		for(unsigned int cel = 0 ; cel < numCells; cel++)
		{
			cells[cel] = cellsTemp[cel];
			if(linkedCellsTemp.size() > 0)
				InsertLinkedCell(linkedCellsTemp[cel],  cel);
			if(linkedFacesTemp.size() > 0)
				InsertLinkedFace(linkedFacesTemp[cel],  cel);
			if(linkedEdgesTemp.size() > 0)
				InsertLinkedEdge(linkedEdgesTemp[cel],  cel);
		}

		return Output::Success;

	}

	// ***************************************************************************
}
