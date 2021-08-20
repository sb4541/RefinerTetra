#ifndef GENERICMESH_H
#define GENERICMESH_H

#include <vector>
#include <map>
#include "Eigen/Eigen"
#include "Output.hpp"
#include "GenericDomain.hpp"

using namespace std;
using namespace Eigen;
using namespace MainApplication;

namespace GeDiM
{
	class GenericDomain;
	class GenericTreeNode;
	class GenericCell;
	class GenericFace;
	class GenericEdge;
	class GenericPoint;
	class GenericMesh;
	class GenericMeshImportInterface;

	class GenericTreeNode
	{
		protected:
			const GenericTreeNode* father; ///< Father node in the tree
			vector<const GenericTreeNode*> childs; ///< Child nodes in the tree.

			unsigned int id; ///< Internal id of the node
			unsigned int globalId; ///< Global id of the node, generally equal to id
			bool isActive; ///< Tells if the node is active for system matrix computation

			vector<unsigned int> markers; ///< Marker for Dirichlet and Neumann conditions

			map<string, void*> properties;

		public:
			GenericTreeNode(const unsigned int& _id);
			GenericTreeNode(const GenericTreeNode& treeNode);
			virtual ~GenericTreeNode();

			const unsigned int& Id() const { return id; }
			const unsigned int& GlobalId() const { return globalId; }
			const bool& IsActive() const { return isActive; }
			const unsigned int NumberOfMarkers() const {return markers.size();}
			const unsigned int& Marker(const unsigned int& position = 0) const { return markers[position]; }
			const vector<unsigned int>& Markers() const { return markers; }

			const bool IsDirichlet(const unsigned int& position = 0) const { return (markers[position] % 2 == 1); }
			const bool IsNeumann(const unsigned int& position = 0) const { return (markers[position] > 0 && markers[position] % 2 == 0); }
			const bool HasFather() const { return (father != NULL);}
			const bool HasChilds() const { return (childs.size() > 0);}
			const unsigned int NumberOfChilds() const { return childs.size();}

			const GenericTreeNode* Father() const{ return father;}
			const GenericTreeNode* Child(const unsigned int& position) const{ return childs[position];}

			void SetGlobalId(const unsigned int& _globalId) { globalId = _globalId; }
			void SetId(const unsigned int& _id) { id = _id; }
			void SetState(bool _isActive = true) { isActive = _isActive; }
			void SetMarker(const unsigned short& _marker, const unsigned int& position = 0) { markers[position] = _marker; }
			void SetSizeMarkers(const unsigned int& numberOfMarker) {markers.resize(numberOfMarker);}
			void SetMarkers(const vector<unsigned int>& _markers) { markers = _markers; }

			void SetFather(const GenericTreeNode* _father) {father = _father;}

			void AllocateChilds(const unsigned int& numChilds){ childs.resize( numChilds);}
			void InsertChild(const GenericTreeNode* child, const unsigned int& position) {childs[position] = child;}

			void InitializeChilds(const unsigned int& numChilds){ childs.reserve( numChilds);}
			void AddChild(const GenericTreeNode* child) {childs.push_back(child);}
			Output::ExitCodes InheritPropertiesByFather();

			const Output::ExitCodes InitializeProperty(const string& key);
			void AddProperty(const string& key, void* value) { properties[key] = value;}
			const map<string, void* >& GetAllProperties() const { return properties; }
			map<string, void* >& GetAllProperties() { return properties; }
			const void* GetProperty(const string& key) const { return properties.at(key); }
			void* GetProperty(const string& key) { return properties.at(key); }
			const bool HasProperty(const string& key) const { return (properties.find(key) != properties.end()); }
	};

	class GenericCell : public GenericTreeNode
	{
		protected:
			vector<const GenericCell*> cells; ///< Cells of the cell
			vector<const GenericFace*> faces; ///< Faces of the cell
			vector<const GenericEdge*> edges; ///< Edges of the cell
			vector<const GenericPoint*> points; ///< Points of the cell

		public:
			GenericCell(const unsigned int& _id);
			GenericCell(const GenericCell& cell);
			virtual ~GenericCell();

			const size_t NumberOfCells() const { return cells.size(); }
			const size_t NumberOfFaces() const { return faces.size(); }
			const size_t NumberOfEdges() const { return edges.size(); }
			const size_t NumberOfPoints() const { return points.size(); }
			const GenericCell* Cell(const unsigned int& position) const { return cells[position]; }
			const GenericFace* Face(const unsigned int& position) const { return faces[position]; }
			const GenericEdge* Edge(const unsigned int& position) const { return edges[position]; }
			const GenericPoint* Point(const unsigned int& position) const { return points[position]; }

			void InitializeCells(const size_t numberOfCells) { cells.reserve(numberOfCells); }
			void InitializeFaces(const size_t numberOfFaces) { faces.reserve(numberOfFaces); }
			void InitializeEdges(const size_t numberOfEdges) { edges.reserve(numberOfEdges); }
			void InitializePoints(const size_t numberOfPoints) { points.reserve(numberOfPoints); }

			void AllocateCells(const size_t numberOfCells) { cells.resize(numberOfCells, NULL); }
			void AllocateFaces(const size_t numberOfFaces) { faces.resize(numberOfFaces, NULL); }
			void AllocateEdges(const size_t numberOfEdges) { edges.resize(numberOfEdges, NULL); }
			void AllocatePoints(const size_t numberOfPoints) { points.resize(numberOfPoints, NULL); }

			void ReAllocateCells(const size_t numberOfCells) { cells.assign(numberOfCells, NULL); }
			void ReAllocateFaces(const size_t numberOfFaces) { faces.assign(numberOfFaces, NULL); }
			void ReAllocateEdges(const size_t numberOfEdges) { edges.assign(numberOfEdges, NULL); }
			void ReAllocatePoints(const size_t numberOfPoints) { points.assign(numberOfPoints, NULL); }

			void ErasePoint(const unsigned int& position) { points.erase(points.begin() + position); }
			void EraseEdge(const unsigned int& position) { edges.erase(edges.begin() + position); }
			void EraseFace(const unsigned int& position) { faces.erase(faces.begin() +position); }
			void EraseCell(const unsigned int& position) { cells.erase(cells.begin() + position); }

			void ShrinkCells() { cells.shrink_to_fit(); }
			void ShrinkFaces() { faces.shrink_to_fit(); }
			void ShrinkEdges() { edges.shrink_to_fit(); }
			void ShrinkPoints() { points.shrink_to_fit(); }

			Output::ExitCodes AddCell(const GenericCell* cell);
			Output::ExitCodes AddFace(const GenericFace* face);
			Output::ExitCodes AddEdge(const GenericEdge* edge);
			Output::ExitCodes AddPoint(const GenericPoint* point);

			Output::ExitCodes InsertCell(const GenericCell* cell, const unsigned int& position);
			Output::ExitCodes InsertFace(const GenericFace* face, const unsigned int& position);
			Output::ExitCodes InsertEdge(const GenericEdge* edge, const unsigned int& position);
			Output::ExitCodes InsertPoint(const GenericPoint* point, const unsigned int& position);

			bool PointInCellAndIdEdgeBoundary (const Vector3d& point, int& numEdges, const double& toll  = 1.0E-7) const;
			bool PointInCell(const Vector3d& point, const double& toll = 1.0E-7) const ;

			const unsigned int RatioMaxMinEdge() const;
			const Output::ExitCodes Centroid(Vector3d& centroid) const;
			const Output::ExitCodes Radius(double& radius, const Vector3d& centroid, const double squaredTolerance = 1.0E-12) const;
			const Output::ExitCodes Measure(double& measure) const;
			const Output::ExitCodes MeasureCentroid(double& measure, Vector3d& centroid) const;

			const Output::ExitCodes ComputeGeometricalProperties(Vector3d& normal, Matrix3d& rotationMatrix, const double& rotationTolerance = 1.0E-7) const;
			void AllignedEdgesPoints();
	};

	class GenericFace : public GenericTreeNode
	{
		protected:
			vector<const GenericCell*> cells; ///< Cells of the face
			vector<const GenericFace*> faces; ///< Faces of the face
			vector<const GenericEdge*> edges; ///< Edges of the face
			vector<const GenericPoint*> points; ///< Points of the face

			double measure;
			Vector3d* centroid;
			vector<Vector3d*> rotatedPoints;

			Vector3d* normal; ///The normal is congruent to the cell in position Zero
			Matrix3d* rotationMatrix;

		public:
			GenericFace(const unsigned int& _id);
			GenericFace(const GenericFace& face);
			virtual ~GenericFace();

			const size_t NumberOfCells() const { return cells.size(); }
			const size_t NumberOfFaces() const { return faces.size(); }
			const size_t NumberOfEdges() const { return edges.size(); }
			const size_t NumberOfPoints() const { return points.size(); }
			const GenericCell* Cell(const unsigned int& position) const { return cells[position]; }
			const GenericFace* Face(const unsigned int& position) const { return faces[position]; }
			const GenericEdge* Edge(const unsigned int& position) const { return edges[position]; }
			const GenericPoint* Point(const unsigned int& position) const { return points[position]; }

			const Matrix3d& RotationMatrix() const {return *rotationMatrix;}
			const double& Measure() const {return measure;}
			const Vector3d& RotatedPointsCoordinate(const unsigned int& position) const {return *rotatedPoints[position];}
			const Vector3d& Centroid() const {return *centroid;}
			const Vector3d& Normal() const  {return *normal;}
			Vector3d& Normal() {return *normal;}

			void SetRotationMatrix(const Matrix3d& rotMatrix){rotationMatrix = new Matrix3d(rotMatrix);}
			void SetMeasure(const double& _measure) {measure = _measure;}
			void SetRotatedPointsCoordinate(const Vector3d& coord, const unsigned int& position) {rotatedPoints[position] = new Vector3d(coord);}
			void SetCentroid(const Vector3d& coord) {centroid = new Vector3d(coord);}
			void SetNormal(const Vector3d& _normal) {normal = new Vector3d(_normal);}

			const int DirectionNormal(unsigned int idCell) const;

			void InitializeCells(const size_t numberOfCells) { cells.reserve(numberOfCells); }
			void InitializeFaces(const size_t numberOfFaces) { faces.reserve(numberOfFaces); }
			void InitializeEdges(const size_t numberOfEdges) { edges.reserve(numberOfEdges); }
			void InitializePoints(const size_t numberOfPoints) { points.reserve(numberOfPoints); }

			void AllocateCells(const size_t numberOfCells) { cells.resize(numberOfCells, NULL); }
			void AllocateFaces(const size_t numberOfFaces) { faces.resize(numberOfFaces, NULL); }
			void AllocateEdges(const size_t numberOfEdges) { edges.resize(numberOfEdges, NULL); }
			void AllocatePoints(const size_t numberOfPoints) { points.resize(numberOfPoints, NULL); }

			void ErasePoint(const unsigned int& position) {points.erase(points.begin() + position);}
			void EraseEdge(const unsigned int& position) {edges.erase(edges.begin() + position);}
			void EraseFace(const unsigned int& position) {faces.erase(faces.begin() +position);}
			void EraseCell(const unsigned int& position) {cells.erase(cells.begin() + position);}

			void ShrinkCells() { cells.shrink_to_fit(); }
			void ShrinkFaces() { faces.shrink_to_fit(); }
			void ShrinkEdges() { edges.shrink_to_fit(); }
			void ShrinkPoints() { points.shrink_to_fit(); }

			Output::ExitCodes AddCell(const GenericCell* cell);
			Output::ExitCodes AddFace(const GenericFace* face);
			Output::ExitCodes AddEdge(const GenericEdge* edge);
			Output::ExitCodes AddPoint(const GenericPoint* point);

			Output::ExitCodes InsertCell(const GenericCell* cell, const unsigned int& position);
			Output::ExitCodes InsertFace(const GenericFace* face, const unsigned int& position);
			Output::ExitCodes InsertEdge(const GenericEdge* edge, const unsigned int& position);
			Output::ExitCodes InsertPoint(const GenericPoint* point, const unsigned int& position);

			Output::ExitCodes InitializeGeometricalProperties();
			bool CheckGeometricalProperties() const {return (normal != NULL || centroid != NULL);}
			Output::ExitCodes ComputeNormal();
			Output::ExitCodes ComputeGeometricalProperties(const double rotationTolerance = 1.0E-7);

	};

	class GenericEdge : public GenericTreeNode
	{
		public:
			enum PositionPoint
			{
				AtTheLeft = 0,
				AtTheRight = 1,
				Beyond = 2,
				Behind = 3,
				Between = 4,
				AtTheOrigin = 5,
				AtTheEnd = 6
			};

		protected:
			vector<const GenericCell*> cells; ///< Cells of the edge
			vector<const GenericFace*> faces; ///< Faces of the edge
			vector<const GenericEdge*> edges; ///< Edges of the edge
			vector<const GenericPoint*> points; ///< Points of the edge

		public:
			GenericEdge(const unsigned int& _id);
			GenericEdge(const unsigned int& _id, const GenericPoint* origin, const GenericPoint* end);
			GenericEdge(const GenericEdge& edge);
			virtual ~GenericEdge();

			const size_t NumberOfCells() const { return cells.size(); }
			const size_t NumberOfFaces() const { return faces.size(); }
			const size_t NumberOfEdges() const { return edges.size(); }
			const size_t NumberOfPoints() const { return points.size(); }

			//position 0 = right position 1 = left
			const GenericCell* Cell(const unsigned int& position) const { return cells[position]; }
			const GenericCell* RightCell() const {return cells[0];}
			const GenericCell* LeftCell() const {return cells[1];}
			const GenericFace* Face(const unsigned int& position) const { return faces[position]; }
			const GenericEdge* Edge(const unsigned int& position) const { return edges[position]; }
			const GenericPoint* Point(const unsigned int& position) const { return points[position]; }

			const bool HasRightCell() const { return (cells[0] != NULL); }
			const bool HasLeftCell() const { return (cells[1] != NULL); }

			void InitializeCells(const size_t numberOfCells) { cells.reserve(numberOfCells); }
			void InitializeFaces(const size_t numberOfFaces) { faces.reserve(numberOfFaces); }
			void InitializeEdges(const size_t numberOfEdges) { edges.reserve(numberOfEdges); }

			void AllocateCells(const size_t numberOfCells) { cells.resize(numberOfCells, NULL); }
			void AllocateFaces(const size_t numberOfFaces) { faces.resize(numberOfFaces, NULL); }
			void AllocateEdges(const size_t numberOfEdges) { edges.resize(numberOfEdges, NULL); }
			void AllocatePoints(const size_t numberOfPoints) { points.resize(numberOfPoints, NULL); }

			void ErasePoint(const unsigned int& position) {points.erase(points.begin() + position);}
			void EraseEdge(const unsigned int& position) {edges.erase(edges.begin() + position);}
			void EraseFace(const unsigned int& position) {faces.erase(faces.begin() +position);}
			void EraseCell(const unsigned int& position) {cells.erase(cells.begin() + position);}

			void ShrinkCells() { cells.shrink_to_fit(); }
			void ShrinkFaces() { faces.shrink_to_fit(); }
			void ShrinkEdges() { edges.shrink_to_fit(); }
			void ShrinkPoints() { points.shrink_to_fit(); }

			Output::ExitCodes AddCell(const GenericCell* cell);
			Output::ExitCodes AddFace(const GenericFace* face);
			Output::ExitCodes AddEdge(const GenericEdge* edge);
			Output::ExitCodes AddPoint(const GenericPoint* point);

			Output::ExitCodes InsertCell(const GenericCell* cell, const unsigned int& position);
			Output::ExitCodes InsertFace(const GenericFace* face, const unsigned int& position);
			Output::ExitCodes InsertEdge(const GenericEdge* edge, const unsigned int& position);
			Output::ExitCodes InsertPoint(const GenericPoint* point, const unsigned int& position);

			Output::ExitCodes ChangeOrientation();
			GenericEdge::PositionPoint PointOnEdge (const Vector3d& point, const double& toll = 1.0E-7) const;

	};

	class GenericPoint : public GenericTreeNode
	{
		protected:
			Vector3d coordinates; ///< Geometry vertex of the degree of freedom. Its size depends on the dimension

			vector<const GenericCell*> cells; ///< Cells of the point
			vector<const GenericFace*> faces; ///< Faces of the point
			vector<const GenericEdge*> edges; ///< Edges of the point

		public:
			GenericPoint(const unsigned int& _id);
			GenericPoint(const GenericPoint& point);
			virtual ~GenericPoint();

			Vector3d& Coordinates() { return coordinates; }
			const Vector3d& Coordinates() const { return coordinates; }
			const double& X() const { return coordinates[0]; }
			const double& Y() const { return coordinates[1]; }
			const double& Z() const { return coordinates[2]; }
			const size_t NumberOfCells() const { return cells.size(); }
			const size_t NumberOfFaces() const { return faces.size(); }
			const size_t NumberOfEdges() const { return edges.size(); }
			const GenericCell* Cell(const unsigned int& position) const { return cells[position]; }
			const GenericFace* Face(const unsigned int& position) const { return faces[position]; }
			const GenericEdge* Edge(const unsigned int& position) const { return edges[position]; }

			Output::ExitCodes SetCoordinates(const Vector3d& _coordinates);
			Output::ExitCodes SetCoordinates(const double& x, const double& y = 0.0, const double& z = 0.0);

			void InitializeCells(const size_t numberOfCells) { cells.reserve(numberOfCells); }
			void InitializeFaces(const size_t numberOfFaces) { faces.reserve(numberOfFaces); }
			void InitializeEdges(const size_t numberOfEdges) { edges.reserve(numberOfEdges); }

			void AllocateCells(const size_t numberOfCells) { cells.resize(numberOfCells, NULL); }
			void AllocateFaces(const size_t numberOfFaces) { faces.resize(numberOfFaces, NULL); }
			void AllocateEdges(const size_t numberOfEdges) { edges.resize(numberOfEdges, NULL); }

			void EraseEdge(const unsigned int& position) {edges.erase(edges.begin() + position);}
			void EraseFace(const unsigned int& position) {faces.erase(faces.begin() +position);}
			void EraseCell(const unsigned int& position) {cells.erase(cells.begin() + position);}

			void ShrinkCells() { cells.shrink_to_fit(); }
			void ShrinkFaces() { faces.shrink_to_fit(); }
			void ShrinkEdges() { edges.shrink_to_fit(); }

			Output::ExitCodes AddCell(const GenericCell* cell);
			Output::ExitCodes AddFace(const GenericFace* face);
			Output::ExitCodes AddEdge(const GenericEdge* edge);

			Output::ExitCodes InsertCell(const GenericCell* cell, const unsigned int& position);
			Output::ExitCodes InsertFace(const GenericFace* face, const unsigned int& position);
			Output::ExitCodes InsertEdge(const GenericEdge* edge, const unsigned int& position);
	};

	class GenericMesh : public IRotation
	{
		protected:
			vector<GenericCell*> cells; ///< Cells of the mesh
			vector<GenericFace*> faces; ///< Faces of the mesh
			vector<GenericEdge*> edges; ///< Edges of the mesh
			vector<GenericPoint*> points; ///< Points of the mesh

		public:
			GenericMesh();
			GenericMesh(const GenericMesh& mesh);
			virtual ~GenericMesh();

			const size_t NumberOfCells() const { return cells.size(); }
			const size_t NumberOfFaces() const { return faces.size(); }
			const size_t NumberOfEdges() const { return edges.size(); }
			const size_t NumberOfPoints() const { return points.size(); }
			const GenericCell* Cell(const unsigned int& position) const { return cells[position]; }
			const GenericFace* Face(const unsigned int& position) const { return faces[position]; }
			const GenericEdge* Edge(const unsigned int& position) const { return edges[position]; }
			const GenericPoint* Point(const unsigned int& position) const { return points[position]; }
			GenericCell* Cell(const unsigned int& position) { return cells[position]; }
			GenericFace* Face(const unsigned int& position) { return faces[position]; }
			GenericEdge* Edge(const unsigned int& position) { return edges[position]; }
			GenericPoint* Point(const unsigned int& position) { return points[position]; }
			const unsigned short Dimension() const { if (faces.size() > 0) return 3; else if (edges.size() > 0) return 2; else return 1; }


			/// Rotate the mesh using the rotation in IRotation
			Output::ExitCodes Rotate(const bool& inverse = false);
			/// Rotate the mesh using matrix Rotation in Input
			Output::ExitCodes RotateWithInput(const Matrix3d& rotationMatrix, const Vector3d& translationVector, const bool& inverse = false);
			Output::ExitCodes RotateCellWithInput(const unsigned int& idCell, const Matrix3d& rotationMatrix, const Vector3d& translationVector, const bool& inverse = false);
			Output::ExitCodes RotateFaceWithInput(const unsigned int& idFace, const Matrix3d& rotationMatrix, const Vector3d& translationVector, const bool& inverse = false);

			virtual GenericCell* CreateCell() { return new GenericCell(cells.size()); }
			virtual GenericFace* CreateFace() { return new GenericFace(faces.size()); }
			virtual GenericEdge* CreateEdge() { return new GenericEdge(edges.size()); }
			virtual GenericPoint* CreatePoint() { return new GenericPoint(points.size()); }

			void InitializeCells(const size_t numberOfCells) { cells.reserve(numberOfCells); }
			void InitializeFaces(const size_t numberOfFaces) { faces.reserve(numberOfFaces); }
			void InitializeEdges(const size_t numberOfEdges) { edges.reserve(numberOfEdges); }
			void InitializePoints(const size_t numberOfPoints) { points.reserve(numberOfPoints); }

			Output::ExitCodes AddCell(GenericCell* cell);
			Output::ExitCodes AddFace(GenericFace* face);
			Output::ExitCodes AddEdge(GenericEdge* edge);
			Output::ExitCodes AddPoint(GenericPoint* point);

			const bool FindCoordinates(const Vector3d& coordinates, unsigned int& idPoint, const double& toll = 1.0E-5);

			const Output::ExitCodes CheckDoublePoints(const double& toll = 1.0E-7);
			const Output::ExitCodes CheckPointsInCells();
			const Output::ExitCodes CheckPointsInFaces();
			const Output::ExitCodes CheckPointsEqualsEdgesInCells2D();
			const Output::ExitCodes CheckPointsEqualsEdgesInFaces();
			const Output::ExitCodes CheckEdgesInFaces();
			const Output::ExitCodes CheckDoubleCells();
			const Output::ExitCodes CheckDoubleFaces();
			const Output::ExitCodes CheckDoubleEdges();
			const Output::ExitCodes CheckNeigs();

			void MovePointMesh(const unsigned int& position,const Vector3d& newCoordinate) { points[position]->SetCoordinates(newCoordinate); }
			const Output::ExitCodes CutEdgeWithPoints(const unsigned int& idEdge, const vector<Vector3d >& coordinatesPoints, const bool& inheritProperty = true);
			const Output::ExitCodes CutEdgeWithCoordinateCurvilinears(const unsigned int& idEdge, const vector<double>& coordinatesCurvilinear, const bool& inheritProperty = true);

			const Output::ExitCodes UpdateFace(const unsigned int& idFace, const int& idEdge = -1);
			const Output::ExitCodes CreateFaceChild(GenericFace& face, GenericFace& faceFather, const list<unsigned int>& idEdgesFace, const list<unsigned int>& idPointFace, const bool& property = true);

			const Output::ExitCodes UpdateCell(const unsigned int& idCell, const int& idEdge = -1);
			const Output::ExitCodes CreateCellChild2D(GenericCell& cell, GenericCell& cellFather, const list<unsigned int>& idEdgesCell, const list<unsigned int>& idPointCell, const bool& property = true);

			const Output::ExitCodes ComputeGeometricalProperties();

			const Output::ExitCodes ActivateFatherNodes();
			const Output::ExitCodes ActivateChildrenNodes();

			virtual const Output::ExitCodes CleanInactiveTreeNode();
	};

	class InterfaceGenericMesh : public GenericMesh
	{
		protected:
			vector< vector<const GenericCell*> > linkedCells;
			vector< vector<const GenericFace*> > linkedFaces;
			vector< vector<const GenericEdge*> > linkedEdges;
			vector< vector<const GenericPoint*> > linkedPoints;

		public:
			InterfaceGenericMesh();
			InterfaceGenericMesh(const InterfaceGenericMesh& mesh);
			virtual ~InterfaceGenericMesh();

			const size_t NumberOfLinkedCells() const { return linkedCells.size(); }
			const size_t NumberOfLinkedFaces() const { return linkedFaces.size(); }
			const size_t NumberOfLinkedEdges() const { return linkedEdges.size(); }
			const size_t NumberOfLinkedPoints() const { return linkedPoints.size(); }

			void InitializeLinkedCells(const size_t numberOfCells) { linkedCells.reserve(numberOfCells); }
			void InitializeLinkedFaces(const size_t numberOfFaces) { linkedFaces.reserve(numberOfFaces); }
			void InitializeLinkedEdges(const size_t numberOfEdges) { linkedEdges.reserve(numberOfEdges); }
			void InitializeLinkedPoints(const size_t numberOfPoints) { linkedPoints.reserve(numberOfPoints); }

			void AllocateLinkedCells(const size_t numberOfCells) { linkedCells.resize(numberOfCells, vector<const GenericCell*>()); }
			void AllocateLinkedFaces(const size_t numberOfFaces) { linkedFaces.resize(numberOfFaces, vector<const GenericFace*>()); }
			void AllocateLinkedEdges(const size_t numberOfEdges) { linkedEdges.resize(numberOfEdges, vector<const GenericEdge*>()); }
			void AllocateLinkedPoints(const size_t numberOfPoints) { linkedPoints.resize(numberOfPoints, vector<const GenericPoint*>()); }

			Output::ExitCodes InsertLinkedCell(const vector<const GenericCell*>& linkCells, const unsigned int& position);
			Output::ExitCodes InsertLinkedFace(const vector<const GenericFace*>& linkFaces, const unsigned int& position);
			Output::ExitCodes InsertLinkedEdge(const vector<const GenericEdge*>& linkEdges, const unsigned int& position);
			Output::ExitCodes InsertLinkedPoint(const vector<const GenericPoint*>& linkPoints, const unsigned int& position);

			Output::ExitCodes AddLinkedCell(const vector<const GenericCell*>& linkCells);
			Output::ExitCodes AddLinkedFace(const vector<const GenericFace*>& linkFaces);
			Output::ExitCodes AddLinkedEdge(const vector<const GenericEdge*>& linkEdges);
			Output::ExitCodes AddLinkedPoint(const vector<const GenericPoint*>& linkPoints);

			const vector<const GenericCell*>& LinkedCell(const unsigned int& position) const { return linkedCells[position]; }
			const vector<const GenericFace*>& LinkedFace(const unsigned int& position) const { return linkedFaces[position]; }
			const vector<const GenericEdge*>& LinkedEdge(const unsigned int& position) const { return linkedEdges[position]; }
			const vector<const GenericPoint*>& LinkedPoint(const unsigned int& position) const { return linkedPoints[position]; }

			vector<const GenericCell*>& LinkedCell(const unsigned int& position) { return linkedCells[position]; }
			vector<const GenericFace*>& LinkedFace(const unsigned int& position) { return linkedFaces[position]; }
			vector<const GenericEdge*>& LinkedEdge(const unsigned int& position) { return linkedEdges[position]; }
			vector<const GenericPoint*>& LinkedPoint(const unsigned int& position) { return linkedPoints[position]; }

			const Output::ExitCodes CleanInactiveTreeNode();
	};

	class GenericMeshImportInterface
	{
		protected:
			double maximumCellSize; ///< Size of the minimum cell of the mesh
			unsigned int minimumNumberOfCells; ///< Minimum number of cell of the mesh
			unsigned int markerDimension; ///< Minimum number of cell of the mesh

			vector<vector<unsigned int> > vertexMarkers; ///< Vector of boundary conditions of domain vertices per equation
			vector<vector<unsigned int> > edgeMarkers; ///< Vector of boundary conditions of domain edges per equation
			vector<vector<unsigned int> > faceMarkers; ///< Vector of boundary conditions of domain faces per equation

			vector<unsigned int> unidimensionalVertexMarkers;
			vector<unsigned int> unidimensionalEdgeMarkers;
			vector<unsigned int> unidimensionalFaceMarkers;

			vector< vector<unsigned int> > mappedMarkers;

			const Output::ExitCodes CreateUniDimMarkers();
			const Output::ExitCodes CreateMappedMarkers();
		public:
			GenericMeshImportInterface();
			virtual ~GenericMeshImportInterface();

			const double& MaximumCellSize() const { return maximumCellSize; }
			const unsigned int& MinimumNumberOfCells() const { return minimumNumberOfCells; }

			const vector<unsigned int>& VertexMarkers(const unsigned int& position = 0) const { return vertexMarkers[position]; }
			const vector<unsigned int>& EdgeMarkers(const unsigned int& position = 0) const { return edgeMarkers[position]; }
			const vector<unsigned int>& FaceMarkers(const unsigned int& position = 0) const { return faceMarkers[position]; }

			const vector<unsigned int> MarkersPerVertex(const unsigned int& position) const;
			const vector<unsigned int> MarkersPerEdge(const unsigned int& position) const;
			const vector<unsigned int> MarkersPerFace(const unsigned int& position) const;

			void SetMarkerDimension(const unsigned int& _markerDimension = 1) {markerDimension = _markerDimension; vertexMarkers.resize(_markerDimension); edgeMarkers.resize(_markerDimension); faceMarkers.resize(_markerDimension);}
			void SetMaximumCellSize(const double& _maximumCellSize) { maximumCellSize = _maximumCellSize; }
			void SetMinimumNumberOfCells(const unsigned int& _minimumNumberOfCells) { minimumNumberOfCells = _minimumNumberOfCells; }

			void SetBoundaryConditions(const vector<unsigned int>& _vertexMarkers, const unsigned int& position = 0);
			void SetBoundaryConditions(const vector<unsigned int>& _vertexMarkers, const vector<unsigned int>& _edgeMarkers, const unsigned int& position = 0);
			void SetBoundaryConditions(const vector<unsigned int>& _vertexMarkers, const vector<unsigned int>& _edgeMarkers, const vector<unsigned int>& _faceMarkers, const unsigned int& position = 0);

			virtual Output::ExitCodes CreateMesh(const GenericDomain& domain, GenericMesh& mesh) const = 0;
	};
}

#endif // GENERICMESH_H
