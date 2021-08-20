//
// Created by Stefano Benedetto on 17/05/2020.
//

#ifndef REFINERTETRA_TEST_H
#define REFINERTETRA_TEST_H

#include "Output.hpp"
#include "GenericMesh.hpp"

namespace GeDiM
{
    class UnitTest
    {
    public:
        static void Test_MaxFace(GenericCell &cell, int id_from_FindMaxFace);
    };
}

#endif //REFINERTETRA_TEST_H

