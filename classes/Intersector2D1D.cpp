#include "Intersector2D1D.hpp"

namespace GeDiM
{
    // ***************************************************************************
    Intersector2D1D::Intersector2D1D()
    {
        planeTranslationPointer = NULL;
        planeNormalPointer = NULL;
        lineOriginPointer = NULL;
        lineTangentPointer = NULL;

        toleranceParallelism = 1.0e-07;
        toleranceIntersection = 1.0e-07;

        intersectionType = Intersector2D1D::NoInteresection;
        intersectionParametricCoordinate = 0.0;
        intersectionPoint.setZero();
    }
    Intersector2D1D::~Intersector2D1D()
    {
        planeTranslationPointer = NULL;
        planeNormalPointer = NULL;
        lineOriginPointer = NULL;
        lineTangentPointer = NULL;
    }
    // ***************************************************************************
    const Output::ExitCodes Intersector2D1D::SetPlane(const Vector3d& planeNormal, const double& planeTranslation)
    {
        Output::ExitCodes result = Output::Success;

        planeNormalPointer = &planeNormal;
        planeTranslationPointer = &planeTranslation;


        return result;
    }
    // ***************************************************************************
    const Output::ExitCodes Intersector2D1D::SetLine(const Vector3d& lineOrigin, const Vector3d& lineTangent)
    {
        Output::ExitCodes result = Output::Success;

        lineOriginPointer = &lineOrigin;
        lineTangentPointer = &lineTangent;


        return result;
    }
    // ***************************************************************************
    const Output::ExitCodes Intersector2D1D::ComputeIntersection()
    {
        Output::ExitCodes result = Output::Success;

        const Vector3d& planeNormal = *planeNormalPointer;
        const double& planeTranslation = *planeTranslationPointer;
        const Vector3d& lineTangent = *lineTangentPointer;
        const Vector3d& lineOrigin = *lineOriginPointer;

        double parallelTest = planeNormal.dot(lineTangent);
        double coplanarTest = planeNormal.dot(lineOrigin);

        if (abs(parallelTest) < toleranceParallelism)
        {
	  /// |*planeNormalPointer . (*lineOriginPointer - *planeOriginPointer)| < toleranceIntersection 
            if (abs(coplanarTest - planeTranslation) < toleranceIntersection)
            {
                intersectionType = Intersector2D1D::Coplanar;
                intersectionParametricCoordinate = 0.0;
            }
            else
            {
                intersectionType = Intersector2D1D::NoInteresection;
                intersectionParametricCoordinate = 0.0;
            }
        }
        else
        {
            intersectionType = Intersector2D1D::PointIntersection;
            intersectionParametricCoordinate = (planeTranslation - coplanarTest) / parallelTest;
            intersectionPoint = lineOrigin + intersectionParametricCoordinate * lineTangent;
        }

        return result;
    }
    // ***************************************************************************
}
