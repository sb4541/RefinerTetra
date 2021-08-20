#ifndef DOMAIN3D_H
#define DOMAIN3D_H

#include <vector>

#include "Eigen/Eigen"
#include "Output.hpp"
#include "GenericDomain.hpp"

using namespace std;
using namespace Eigen;
using namespace MainApplication;

namespace GeDiM
{
	class GenericDomain3D;
	class Parallelepiped;

	/// Parallelepiped as defined in https://en.wikipedia.org/wiki/Parallelepiped
	class Parallelepiped : public GenericDomain3D
	{
		protected:
			unsigned int baseIndex; ///< Indices of the base face
			unsigned int heightIndex; ///< Indices of the height edge

		public:
			Parallelepiped(const unsigned int& _globalId);
			virtual ~Parallelepiped();

			Output::ExitCodes BuildParallelepiped(const Vector3d& origin, const Vector3d& lengthVector, const Vector3d& heightVector, const Vector3d& widthVector, const double& volumeTolerance = 1.0e-16);
	};
}

#endif // DOMAIN3D_H
