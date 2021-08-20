#include "Intersector2D2D.hpp"

namespace GeDiM
{
    // ***************************************************************************
    Intersector2D2D::Intersector2D2D()
    {
        toleranceParallelism = 1.0e-07;
        toleranceIntersection = 1.0e-07;

        intersectionType = Intersector2D2D::NoInteresection;
        rightHandSide.setZero();
    }


    Intersector2D2D::~Intersector2D2D()
    {

    }
    // ***************************************************************************
    const Output::ExitCodes Intersector2D2D::SetFirstPlane(const Vector3d& planeNormal, const double& planeTranslation)
    {
        matrixNomalVector.row(0) = planeNormal;

        rightHandSide(0) = planeTranslation;

        return Output::Success;
    }

    // ***************************************************************************
    const Output::ExitCodes Intersector2D2D::SetSecondPlane(const Vector3d& planeNormal, const double& planeTranslation)
    {
        matrixNomalVector.row(1) = planeNormal;

        rightHandSide(1) = planeTranslation;

        return Output::Success;
    }

    // ***************************************************************************
    const Output::ExitCodes Intersector2D2D::ComputeIntersection()
    {
        tangentLine = matrixNomalVector.row(0).cross(matrixNomalVector.row(1));
        if(tangentLine.squaredNorm() < toleranceParallelism * toleranceParallelism)
        {

	  /// |d_0 - sign(N(0).N(1))*d_1| = |N(0).(X_0-X_1)|< tol
	  double samedirectionNormalVectors = (matrixNomalVector.row(0).dot(matrixNomalVector.row(1))>0.0 ? 1.0 : -1.0);
	    if( abs(rightHandSide(0) - samedirectionNormalVectors * rightHandSide(1)) < toleranceIntersection )
            {
                intersectionType = Intersector2D2D::Coplanar;
                Output::PrintGenericMessage("Coplanar domains", true);
                return Output::Success;
            }
            else
            {
                intersectionType = Intersector2D2D::NoInteresection;
                Output::PrintWarningMessage("Skew normal vector", true);
                return Output::Success;
            }
        }

	/// We look for the point that is common to the two given planes
	/// and the plane orthogonal to them and containing the origin
	/// N(0).X=N(0).X_0=d_0
	/// N(1).X=N(1).X_1=d_1
	/// T   .X=T*O=0
        matrixNomalVector.row(2) = tangentLine.normalized();
        pointLine = matrixNomalVector.colPivHouseholderQr().solve(rightHandSide);
        intersectionType = Intersector2D2D::LineIntersection;
        return Output::Success;
    }

    // ***************************************************************************
}
