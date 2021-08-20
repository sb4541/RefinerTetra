#ifndef INTERSECTOR1D1D_HPP
#define INTERSECTOR1D1D_HPP

#include "Output.hpp"
#include "Eigen/Eigen"
#include "list"

using namespace std;
using namespace MainApplication;
using namespace Eigen;

namespace GeDiM
{
    class Intersector1D1D;

    class Intersector1D1D
    {
        public:
            enum Type
            {
                NoIntersection = 0,
                IntersectionOnLine = 1,
                IntersectionOnSegment = 2,
                IntersectionParallelOnLine = 3,
                IntersectionParallelOnSegment = 4,
            };
            enum Position
            {
                Begin = 0,
                Inside = 1,
                End = 2,
                Outside = 3
            };

        private:

            double toleranceParallelism;
            double toleranceIntersection;

            Type type;
            Position positionIntersectionFirstEdge;
            Position positionIntersectionSecondEdge;

            Vector2d resultParametricCoordinates;
            Matrix2d matrixTangentVector;

        public:
            Intersector1D1D();
            ~Intersector1D1D();

            void SetToleranceParallelism(const double& _tolerance) { toleranceParallelism = _tolerance; }
            void SetToleranceIntersection(const double& _tolerance) { toleranceIntersection = _tolerance; }

            void SetMatrixTangentVector(const Vector3d& _vector0, const Vector3d& _vector1) {matrixTangentVector.col(0) = _vector0.head(2); matrixTangentVector.col(1) = -_vector1.head(2);}
            void SetFirstTangentVector(const Vector3d& _vector) {matrixTangentVector.col(0) = _vector.head(2);}
            void SetSecondTangentVector(const Vector3d& _vector) {matrixTangentVector.col(1) = -1.0 * _vector.head(2);}

            ///Compute the intersections between the line of the first edge and the second edge and returns the parametric coordinates and the intersection type
            /// line_1: X=X1_i+T1*s_0
            /// line_2: X=X2_i+T2*s_1
            /// tanVectorFirstEdg = T1
            /// tanVectorSecondEdg = T2
            /// tanVectorDifference = X2_i-X1_i
            ///The first parametric coordinate refers to the first tangentVector and ...
            ///@note In case of parallelism both the parametric coordinates refers to the first edge
            ///@note Before the call of this member the call of SetFirstTangentVector and  SetSecondTangentVector  or SetMatrixTangentVector are needed
            const Output::ExitCodes ComputeIntersectionEdges(const Vector3d& tanVectorFirstEdge, const Vector3d& tanVectorSecondEdge, const Vector3d& tanVectorDifference);

            const double& ToleranceIntersection() const {return toleranceIntersection; }
            const double& ToleranceParallelism() const {return toleranceParallelism; }

            const Vector2d& ParametricCoordinates() {return resultParametricCoordinates;}
            const double& FirstParametricCoordinate() {return resultParametricCoordinates(0);}
            const double& SecondParametricCoordinate() {return resultParametricCoordinates(1);}

            const Position& PositionIntersectionInFirstEdge(){return positionIntersectionFirstEdge;}
            const Position& PositionIntersectionInSecondEdge(){return positionIntersectionSecondEdge;}
            const Type& TypeIntersection() {return type;}

            //			const Vector3d IntersectionStartPointFirstVector(const Vector3d& startPoint) {return startPoint + resultParametricCoordinates(0) * matrixTangentVector.col(0);}
            //			const Vector3d IntersectionStartPointSecondVector(const Vector3d& startPoint) {return startPoint - resultParametricCoordinates(1) * matrixTangentVector.col(1);}
            const Vector3d IntersectionStartPointFirstVector(const Vector3d& startPoint,const Vector3d& tanVector) {return startPoint + resultParametricCoordinates(0) * tanVector;}
            const Vector3d IntersectionStartPointSecondVector(const Vector3d& startPoint,const Vector3d& tanVector) {return startPoint + resultParametricCoordinates(1) * tanVector;}

            static bool ParamCoordSetIntersection(vector<double>& firstParamCoord, vector<double>& secondParamCoord, vector<double>& set);
    };
}

#endif // INTERSECTOR1D1D_HPP

