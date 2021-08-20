//
// Created by Stefano Benedetto on 17/05/2020.
//
#include "Test.h"

namespace GeDiM
{
    void UnitTest::Test_MaxFace(GenericCell &cell, int id_from_FindMaxFace)
    {
        int id_test = 0;
        double maxAreaTest = 0.0;
        double currentArea= 0.0;
        Vector3d first_edge;
        Vector3d second_edge;

        for (int position_face = 0; position_face < cell.NumberOfFaces(); ++position_face)
        {
            const GenericFace& faceToTreat = *cell.Face(position_face);
            first_edge = faceToTreat.Point(0)->Coordinates();
            second_edge = faceToTreat.Point(faceToTreat.NumberOfPoints()-1)->Coordinates();
            currentArea = 0.5 * first_edge.cross(second_edge).norm();
            if (currentArea > maxAreaTest)
            {
                maxAreaTest = currentArea;
                id_test = position_face;
            }
        }

        if (id_from_FindMaxFace != id_test)
        {
            Output::PrintErrorMessage("Something went wrong in finding the biggest face\n", false);
        }
        else
        {
            Output::Assert(true, "%d is the correct id\n", id_from_FindMaxFace);
        }

        //oppure
        /*
         * const Output::ExitCodes RefinerTetra::FindMaxFace(GenericCell& cell, int& idMaxFace)
    {
        double MaxAreaFace = 0.0;
        int position = 0;
        int id = 2;
        int totalNumberofVertices = 3;
        GenericDomain2D* faceToTreatasDomain = new GenericDomain2D(id, totalNumberofVertices);
        const GenericFace* faceToTreat;

        for (position = 0; position < cell.NumberOfFaces(); ++position)
        {
            faceToTreat = cell.Face(position);
            for (int number_vertex = 0; number_vertex < faceToTreat->NumberOfPoints(); ++number_vertex)
            {
                faceToTreatasDomain->AddVertex(faceToTreat->Point(number_vertex)->Coordinates());
            }
            faceToTreatasDomain->Initialize();
            faceToTreatasDomain->ComputeAreaAndCentroid();
            if (faceToTreatasDomain->Area() > MaxAreaFace)
            {
                idMaxFace = position;
            }
        }

        return Output::Success;
    }
         */

    }

}