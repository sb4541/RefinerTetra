#include "Intersector1D1D.hpp"

namespace GeDiM
{
    Intersector1D1D::Intersector1D1D()
    {
        toleranceParallelism = 1.0e-05;
        toleranceIntersection = 1.0e-07;
        type = Intersector1D1D::NoIntersection;
        resultParametricCoordinates.setZero(2);
    }
    Intersector1D1D::~Intersector1D1D()
    {

    }
    // ***************************************************************************
    bool Intersector1D1D::ParamCoordSetIntersection(vector<double>& firstParamCoord, vector<double>& secondParamCoord, vector<double>& set)
    {
        set[0] = MAX(firstParamCoord[0], secondParamCoord[0]);
        set[1] = MIN(firstParamCoord[1], secondParamCoord[1]);
        return set[0] < set[1];
    }

    // ***************************************************************************
    const Output::ExitCodes Intersector1D1D::ComputeIntersectionEdges(const Vector3d& tanVectorFirstEdge, const Vector3d& tanVectorSecondEdge, const Vector3d& tanVectorDifference)
    {
      /// parallelism = ||T1/\T2||
        double parallelism = fabs(tanVectorFirstEdge.x() * tanVectorSecondEdge.y() - tanVectorSecondEdge.x() * tanVectorFirstEdge.y());
        type = NoIntersection;

        double check = toleranceParallelism * toleranceParallelism * tanVectorFirstEdge.squaredNorm() * tanVectorSecondEdge.squaredNorm();

	/// ||T1/\T2||^2 >= tol^2 * ||T1||^2 * ||T2||^2
        if(parallelism * parallelism >= check)
        {
            /// <li> If the edge and the trace are not parallel look for the intersection with parametric coordinates
            FullPivHouseholderQR<Matrix<double, 2, 2>> qrFactorization = matrixTangentVector.fullPivHouseholderQr();
            resultParametricCoordinates = qrFactorization.solve(tanVectorDifference.head(2));

            if (resultParametricCoordinates(1) > -toleranceIntersection  && resultParametricCoordinates(1)-1.0 < toleranceIntersection)
            {
                type = IntersectionOnLine;
                if (resultParametricCoordinates(0) > -toleranceIntersection  && resultParametricCoordinates(0)-1.0 < toleranceIntersection)
                    type = IntersectionOnSegment;
            }
        }
        else
        {
	    /// parallelis2 = ||T1/\(X2_i-X1_i)||
            double parallelism2 = fabs(tanVectorFirstEdge.x() * tanVectorDifference.y() - tanVectorDifference.x() * tanVectorFirstEdge.y());
            /// <li> In case of parallelism check if the segment is the same with parametric coordinates

            double check2 = toleranceParallelism * toleranceParallelism * tanVectorFirstEdge.squaredNorm() * tanVectorDifference.squaredNorm();
	    /// ||T1/\(X2_i-X1_i)||^2 >= tol^2 * ||T1||^2 * ||X2_i-X1_i||^2
            if( parallelism2 * parallelism2 <= check2 )
            {
                double tempNorm = 1.0/(tanVectorFirstEdge.squaredNorm());
                // parametric coordinates on the trace of the starting point and end point
                resultParametricCoordinates(0) = tanVectorFirstEdge.dot(tanVectorDifference) * tempNorm ;
                resultParametricCoordinates(1) = resultParametricCoordinates(0) + tanVectorFirstEdge.dot(tanVectorSecondEdge) * tempNorm;

                type = IntersectionParallelOnLine;

                if(resultParametricCoordinates(1) < resultParametricCoordinates(0))
                {
                    double tmp = resultParametricCoordinates(0);
                    resultParametricCoordinates(0) = resultParametricCoordinates(1);
                    resultParametricCoordinates(1) = tmp;
                }
                // if one vertex is in the edge there is the intersection
		///  -toleranceIntersection < resultParametricCoordinates(0) < 1.0+toleranceIntersection
                if( (resultParametricCoordinates(0) > -toleranceIntersection && resultParametricCoordinates(0)-1.0 < toleranceIntersection) ||
                        (resultParametricCoordinates(1) > -toleranceIntersection && resultParametricCoordinates(1)-1.0 < toleranceIntersection)   )
                    type = IntersectionParallelOnSegment;
                else
                {
                    //IL PRIMO SEGMENTO DATO IN INPUT E' CONTENUTO NEL SECONDO
		  /// resultParametricCoordinates(0) < toleranceIntersection && resultParametricCoordinates(1) > 1.0-toleranceIntersection
                    if( ( resultParametricCoordinates(0) < toleranceIntersection && resultParametricCoordinates(1) - 1.0 > -toleranceIntersection) )
                        type = IntersectionParallelOnSegment;
                }
            }
        }
        /// </ul>
        if(resultParametricCoordinates(0) < -toleranceIntersection || resultParametricCoordinates(0) > 1.0 + toleranceIntersection)
            positionIntersectionFirstEdge =  Outside;
        else if((resultParametricCoordinates(0) > -toleranceIntersection) && (resultParametricCoordinates(0) < toleranceIntersection))
        {
            resultParametricCoordinates(0) = 0.0;
            positionIntersectionFirstEdge = Begin;
        }
        else if ((resultParametricCoordinates(0) > 1.0 - toleranceIntersection) && (resultParametricCoordinates(0) < 1.0 + toleranceIntersection))
        {
            resultParametricCoordinates(0) = 1.0;
            positionIntersectionFirstEdge = End;
        }
        else
            positionIntersectionFirstEdge = Inside;

        if(resultParametricCoordinates(1) < -toleranceIntersection || resultParametricCoordinates(1) > 1.0 + toleranceIntersection)
            positionIntersectionSecondEdge =  Outside;
        else if((resultParametricCoordinates(1) > -toleranceIntersection) && (resultParametricCoordinates(1) < toleranceIntersection))
        {
            resultParametricCoordinates(1) = 0.0;
            positionIntersectionSecondEdge = Begin;
        }
        else if ((resultParametricCoordinates(1) > 1.0 - toleranceIntersection) && (resultParametricCoordinates(1) <= 1.0 + toleranceIntersection))
        {
            resultParametricCoordinates(1) = 1.0;
            positionIntersectionSecondEdge = End;
        }
        else
            positionIntersectionSecondEdge = Inside;
        return Output::Success;
    }
    // ***************************************************************************
}
