#ifndef INTERSECTOR2D1D_HPP
#define INTERSECTOR2D1D_HPP

#include "Output.hpp"
#include "Eigen/Eigen"
#include "list"

using namespace std;
using namespace MainApplication;
using namespace Eigen;

namespace GeDiM
{
    class Intersector2D1D;

    class Intersector2D1D
    {
        public:
            enum TypeIntersection
            {
                NoInteresection = 0,
                Coplanar = 1,
                PointIntersection = 2
            };

        protected:
            double toleranceParallelism;
            double toleranceIntersection;

            /// *planeTranslationPointer = d = *planeNormalPointer . *planeOriginPointer
            const double* planeTranslationPointer;
            const Vector3d* planeNormalPointer;
            const Vector3d* lineOriginPointer;
            const Vector3d* lineTangentPointer;

            TypeIntersection intersectionType;
            double intersectionParametricCoordinate;
            Vector3d intersectionPoint;

        public:
            Intersector2D1D();
            ~Intersector2D1D();

            void SetToleranceIntersection(const double& _tolerance) { toleranceIntersection = _tolerance; }
            void SetToleranceParallelism(const double& _tolerance) { toleranceParallelism = _tolerance; }

            const double& ToleranceIntersection() const {return toleranceIntersection; }
            const double& ToleranceParallelism() const {return toleranceParallelism; }
            const TypeIntersection& IntersectionType() const { return intersectionType; }
            const double& IntersectionParametricCoordinate() const { return intersectionParametricCoordinate; }
            const Vector3d& IntersectionPoint() const { return intersectionPoint; }

            const Output::ExitCodes SetPlane(const Vector3d& planeNormal, const double& planeTranslation);
            const Output::ExitCodes SetLine(const Vector3d& lineOrigin, const Vector3d& lineTangent);
            const Output::ExitCodes ComputeIntersection();
    };

}

#endif // INTERSECTOR2D1D_HPP

