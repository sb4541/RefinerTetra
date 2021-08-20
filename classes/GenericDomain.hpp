#ifndef GENERICDOMAIN_H
#define GENERICDOMAIN_H

#include <vector>
#include <string>

#include "Eigen/Eigen"
#include "Output.hpp"

using namespace std;
using namespace Eigen;
using namespace MainApplication;

namespace GeDiM
{
	class IRotation;
	class GenericDomain;
	class GenericDomain1D;
	class GenericDomain2D;
	class GenericDomain3D;

	class IRotation
	{
		protected:
			Matrix3d* rotation; ///< rotation matrix
			Vector3d* originTranslation; ///< origin point of the translation

		public:
			IRotation();
			IRotation(const IRotation& irotation);
			virtual ~IRotation();

			const Matrix3d& RotationMatrix() const { return *rotation; }
			const Vector3d& Translation() const { return *originTranslation; }
			Matrix3d& RotationMatrix() { return *rotation; }
			Vector3d& Translation() { return *originTranslation; }

			void SetRotationMatrix(const Matrix3d& _rotationMatrix) { rotation = new Matrix3d(_rotationMatrix); }
			void SetOriginTranslation(const Vector3d& _originTranslation) { originTranslation = new Vector3d(_originTranslation); }

			Output::ExitCodes RotateVertices(const vector<Vector3d>& verticesToRotate, vector<Vector3d>& rotatedVertices) const;
			Vector3d RotatePoint(const Vector3d& point, const bool& translation = true, const bool& inverse = false) const;
			Vector3d RotatePointWithInput(const Vector3d& point, const Matrix3d& rotationMatrix, const Vector3d& translationVector, const bool& inverse = false) const;
	};


	class GenericDomain
	{
		protected:
			unsigned int globalId; ///< Id of the domain
			string description; ///< Description of the domain
			bool initialized; ///< True if the function Initialize finish with success
			bool toBeExport;

			unsigned int totalNumberVertices; ///< Number of vertices of the domain
			unsigned int totalNumberEdges; ///< Number of edges of the domain
			unsigned int totalNumberFaces; ///< Number of faces of the domain

			vector<Vector3d> vertices; ///< Array of domain vertices
			vector<unsigned int> edgeVertices; ///< Vertex indices for each edge of the domain

			vector< vector<unsigned int> > faceVertices; ///< Vertex indices for each face of the domain
			vector< vector<unsigned int> > faceEdges; ///< Edge indices for each face of the domain

			Matrix3d rotationMatrix; ///< rotation matrix
			vector<Vector3d> rotatedVertices; ///< Array of vertexes rotated

			double measure; ///< Measure of the domain

			Vector3d centroid; ///< Centroid of the domain
			double squaredRadius; ///< Squared radius of the domain (squared max distance between centroid and vertices)
			double radius; ///< Radius of the domain (max distance between centroid and vertices)

			map<string, void*> properties;

			void ConstructorBase(const unsigned int& _globalId);
		public:
			GenericDomain(const unsigned int& _globalId);
			GenericDomain(const unsigned int& _globalId, const unsigned int& _totalNumberVertices);
			GenericDomain(const unsigned int& _globalId, const unsigned int& _totalNumberVertices, const unsigned int& _totalNumberEdges);
			GenericDomain(const unsigned int& _globalId, const unsigned int& _totalNumberVertices, const unsigned int& _totalNumberEdges, const unsigned int& _totalNumberFaces);
			GenericDomain(const GenericDomain& domain);
			virtual ~GenericDomain();

			const unsigned int& GlobalId() const { return globalId; }
			const string& Description() const { return description; }

			const unsigned int& TotalNumberVertices() const { return totalNumberVertices; }
			const size_t NumberVertices() const { return vertices.size(); }
			const Vector3d& Vertex(const unsigned int& i) const { return vertices[i]; }

			const unsigned int& TotalNumberEdges() const { return totalNumberEdges; }
			const size_t NumberEdges() const { return edgeVertices.size() * 0.5; }
			const unsigned int& EdgeOriginIndex(const unsigned int& edgeNumber) const { return edgeVertices[2 * edgeNumber]; }
			const unsigned int& EdgeEndIndex(const unsigned int& edgeNumber) const { return edgeVertices[2 * edgeNumber + 1]; }

			const unsigned int& TotalNumberFaces() const { return totalNumberFaces; }
			const size_t NumberFaces() const { return faceVertices.size(); }
			const size_t NumberFacePoints(const unsigned int& faceNumber) const { return faceVertices[faceNumber].size(); }
			const size_t NumberFaceEdges(const unsigned int& faceNumber) const { return faceEdges[faceNumber].size(); }

			const vector<unsigned int>& StreamFacePoints(const unsigned int& faceNumber) {return faceVertices[faceNumber];}
			const unsigned int FacePointIndex(const unsigned int& faceNumber, const unsigned int& pointNumber) const { return faceVertices[faceNumber][pointNumber]; }
			const unsigned int FaceEdgeIndex(const unsigned int& faceNumber, const unsigned int& edgeNumber) const { return faceEdges[faceNumber][edgeNumber]; }

			const size_t NumberRotatedVertices() const { return rotatedVertices.size(); }
			const Vector3d& RotatedVertex(const unsigned int& i) const { return rotatedVertices[i]; }

			const double& Measure() const { return measure; }
			const Vector3d& Centroid() const { return centroid; }
			const double& Radius() const { return radius; }
			const double& SquaredRadius() const { return squaredRadius; }
			virtual const unsigned int Dimension() const = 0;
			const bool& IsInitialized() const { return initialized; }

			/// Set export options.
			const bool& IsToBeExport() const { return toBeExport; }
			void SetToBeExport(const bool& _toBeExport) { toBeExport = _toBeExport; }

			Output::ExitCodes ComputeRadius(const double& squaredTolerance = 1.0e-12);
			virtual Output::ExitCodes ComputeRotationMatrix(const double& rotationTolerance = 1.0e-12) = 0;

			void SetDescription(const string& _description) { description = _description; }
			virtual Output::ExitCodes Initialize() = 0;

			void InitializeVertices(const unsigned int& _totalNumberVertices) { totalNumberVertices =_totalNumberVertices; vertices.reserve(_totalNumberVertices); }
			Output::ExitCodes AddVertex(const Vector3d& vertex);

			void InitializeEdges(const unsigned int& _totalNumberEdges) { totalNumberEdges =_totalNumberEdges; edgeVertices.reserve(2 * _totalNumberEdges); }
			Output::ExitCodes AddEdge(const unsigned int& originNumber, const unsigned int& endNumber);

			void InitializeFaces(const unsigned int& _totalNumberFaces) { totalNumberFaces =_totalNumberFaces;faceVertices.reserve(_totalNumberFaces); faceEdges.reserve(_totalNumberFaces); }
			Output::ExitCodes AddFace(vector<unsigned int> _faceVertices, vector<unsigned int> _faceEdges);

			Output::ExitCodes InitializeProperty(const string& key);
			void AddProperty(const string& key, void* value) { properties[key] = value; }
			const map<string, void* >& GetAllProperties() const { return properties; }
			const void* GetProperty(const string& key) const { return properties.at(key); }
			void* GetProperty(const string& key) { return properties.at(key); }
			const bool HasProperty(const string& key) const { return (properties.find(key) != properties.end()); }

			Output::ExitCodes CopyPropertyFromDomain(const GenericDomain& domain);
	};

	class GenericDomain0D: public GenericDomain, public IRotation
	{
		public:

			GenericDomain0D(const unsigned int& _globalId, const Vector3d& coordinates);
			GenericDomain0D(const unsigned int& _globalId, const bool& initializeVertices = false);
			GenericDomain0D(const GenericDomain0D& domain);
			virtual ~GenericDomain0D();

			virtual const unsigned int Dimension() const { return 0; }

			const Vector3d& Coordinates() const { return vertices[0]; }

			virtual Output::ExitCodes ComputeRotationMatrix(const double& rotationTolerance = 1.0e-12) { return Output::Success;}

			virtual Output::ExitCodes Initialize() { return Output::Success;}
	};

	class GenericDomain1D : public GenericDomain, public IRotation
	{
		protected:
			Vector3d tangent; ///< Tangent of the domain in the space: it's norm is one

		public:
			GenericDomain1D(const unsigned int& _globalId, const bool& initializeVertices = false);
			GenericDomain1D(const GenericDomain1D& domain);
			virtual ~GenericDomain1D();

			const double& Length() const { return measure; }
			virtual const unsigned int Dimension() const { return 1; }

			const Vector3d& Origin() const { return vertices[0]; }
			const Vector3d& End() const { return vertices[1]; }
			const Vector3d Tangent() const { return tangent; }

			virtual Output::ExitCodes Initialize();
			virtual Output::ExitCodes ComputeRotationMatrix(const double& rotationTolerance = 1.0e-12);

			Output::ExitCodes ComputeTangent();
			Output::ExitCodes ComputeLengthAndCentroid(const double& lengthTolerance = 1e-16);

			/// Compute a point on the domain using the curvilinear abscissa
			void ComputePointOnLine(const double parametricCoordinate, Vector3d& point) const { point = tangent * parametricCoordinate * measure; point.noalias() += Origin(); }
	};

	class GenericDomain2D : public GenericDomain, public IRotation
	{
		public:
			enum PositionPoint
			{
				AllNegative = -1,
				Mixed = 0,
				AllPositive = 1,
			};

		protected:
			Vector3d planeNormal; ///< Normal of plane containing the domain in the space: comes from the equation ax+by+cz = d
			double planeTranslation; ///< plane translation in the space: comes from the equation ax+by+cz = d

		public:
			GenericDomain2D(const unsigned int& _globalId);
			GenericDomain2D(const unsigned int& _globalId, const unsigned int& _totalNumberVertices);
			GenericDomain2D(const GenericDomain2D& domain);
			virtual ~GenericDomain2D();

			const double& Area() const { return measure; }
			virtual const unsigned int Dimension() const { return 2; }
			const Vector3d& PlaneNormal() const { return planeNormal; }
			const double& PlaneTranslation() const { return planeTranslation; }

			virtual Output::ExitCodes Initialize();
			virtual Output::ExitCodes ComputeRotationMatrix(const double& rotationTolerance = 1.0e-12);

			Output::ExitCodes ComputePlane();
			Output::ExitCodes ComputeAreaAndCentroid(const double& areaTolerance = 1.0e-16);
			/// In the rotated domain if the point is inside or outside the domain. The border is considered outside.
			const bool PointInDomainRotatedNoBoundary (const Vector3d& point, const double& tolerance = 1.0e-7) const;
			const PositionPoint PositionPointDomain(const GenericDomain2D& domain, const double& tolerance = 1.0e-7) const;
	};

	class GenericDomain3D : public GenericDomain
	{
		protected:

		public:
			GenericDomain3D(const unsigned int& _globalId);
			GenericDomain3D(const unsigned int& _globalId, const unsigned int& _totalNumberVertices, const unsigned int& _totalNumberEdges, const unsigned int& _totalNumberFaces);
			virtual ~GenericDomain3D();

			const double& Volume() const { return measure; }
			virtual const unsigned int Dimension() const { return 3; }

			virtual Output::ExitCodes Initialize() { return ComputeMeausureAndCentroid(); }
			Output::ExitCodes ComputeMeausureAndCentroid(const double& rotationTolerance = 1.0e-12);
			virtual Output::ExitCodes ComputeRotationMatrix(const double& rotationTolerance = 1.0e-12) { return Output::UnimplementedMethod; }
	};
}

#endif // GENERICDOMAIN_H
