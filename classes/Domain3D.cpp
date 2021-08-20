#include "Domain3D.hpp"
#include "GenericDomain.hpp"
#include "Output.hpp"

namespace GeDiM
{
	// ***************************************************************************
	// Parallelepiped Implementation
	// ***************************************************************************
	Parallelepiped::Parallelepiped(const unsigned int& _globalId) : GenericDomain3D(_globalId, 8, 12, 6)
	{
		baseIndex = -1;
		heightIndex = -1;
	}
	Parallelepiped::~Parallelepiped()
	{
	}
	// ***************************************************************************
	Output::ExitCodes Parallelepiped::BuildParallelepiped(const Vector3d& origin, const Vector3d& lengthVector, const Vector3d& heightVector, const Vector3d& widthVector, const double& volumeTolerance)
	{
		/// <ul>
		if (abs(lengthVector.cross(widthVector).norm()) < 1e-12)
		{
			Output::PrintErrorMessage("Length and Width vectors are parallel, no parallelepiped can be created", false);
			return Output::GenericError;
		}

		if (abs(lengthVector.cross(heightVector).norm()) < 1e-12)
		{
			Output::PrintErrorMessage("Length and Height vectors are parallel, no parallelepiped can be created", false);
			return Output::GenericError;
		}

		if (abs(widthVector.cross(heightVector).norm()) < 1e-12)
		{
			Output::PrintErrorMessage("Width and Height vectors are parallel, no parallelepiped can be created", false);
			return Output::GenericError;
		}

		/// <li> Add base
		AddVertex(origin);
		AddVertex(origin + lengthVector);
		AddVertex(origin + lengthVector + widthVector);
		AddVertex(origin + widthVector);

		AddEdge(0, 1);
		AddEdge(1, 2);
		AddEdge(2, 3);
		AddEdge(3, 0);

		AddFace(vector<unsigned int> { 0, 3, 2, 1 }, vector<unsigned int> { 0, 3, 2, 1 });

		/// <li> Add heigth
		AddVertex(origin + heightVector);
		AddVertex(origin + heightVector + lengthVector);
		AddVertex(origin + heightVector + lengthVector + widthVector);
		AddVertex(origin + heightVector + widthVector);

		AddEdge(0, 4);
		AddEdge(1, 5);
		AddEdge(2, 6);
		AddEdge(3, 7);

		AddFace(vector<unsigned int> { 0, 1, 5, 4}, vector<unsigned int> { 0, 5, 8, 4 });
		AddFace(vector<unsigned int> { 1, 2, 6, 5 }, vector<unsigned int> { 1, 6, 9, 5 });
		AddFace(vector<unsigned int> { 2, 3, 7, 6 }, vector<unsigned int> { 2, 7, 10, 6 });
		AddFace(vector<unsigned int> { 0, 4, 7, 3 }, vector<unsigned int> { 3, 7, 11, 4 });

		AddEdge(4, 5);
		AddEdge(5, 6);
		AddEdge(6, 7);
		AddEdge(7, 4);

		AddFace(vector<unsigned int> { 4, 5, 6, 7 }, vector<unsigned int> { 8, 9, 10, 11 });

		baseIndex = 0;
		heightIndex = 4;

		/// <li> Compute Volume and Centroid
		measure = abs(heightVector.dot(lengthVector.cross(widthVector)));
		centroid = origin + lengthVector * 0.5 + heightVector * 0.5 + widthVector * 0.5;

		if (fabs(measure) < volumeTolerance)
		{
			Output::PrintErrorMessage("Parallelepiped %d has area too small: %.2e - tolerance %.2e", false, globalId, measure, volumeTolerance);
			return Output::GenericError;
		}

		/// <li> Compute Radius
		Output::ExitCodes result = ComputeRadius();
		if (result != Output::Success)
		{
			Output::PrintErrorMessage("Domain %d initialization. ComputeRadius failed.", false, globalId);
			return result;
		}

		return Output::Success;
		/// <ul>
	}
	// ***************************************************************************
}
