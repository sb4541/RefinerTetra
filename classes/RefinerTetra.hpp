#ifndef REFINERTETRA_HPP
#define REFINERTETRA_HPP

#include "GenericMesh.hpp"
#include "vector"
#include "list"
#include "Output.hpp"
#include "Cutter3D.hpp"

namespace GeDiM
{
    class RefinerTetra
    {
        protected:

            GenericMesh* meshPointer;
            vector<unsigned int> idCellToRefine;
            list<unsigned int> idEdgesCut;
            //vector<bool> check_presence;
            int refiner_counter;

            const Output::ExitCodes FindMaxFace(GenericCell& cell, unsigned int& idMaxFace);
            const Output::ExitCodes CutTetra(GenericCell& cell, bool& to_add);
            const Output::ExitCodes RecoverConformity();
            //const Output::ExitCodes FindLongestEdge(const GenericCell& cell, unsigned int& idLongestEdge);

        public:
            RefinerTetra();
            virtual ~RefinerTetra();

            const Output::ExitCodes FindLongestEdge(const GenericCell& cell, unsigned int& idLongestEdge);
            //void SetVector();
            const Output::ExitCodes CheckConformity();
            inline void SetCounter() {refiner_counter = 0;};
            const Output::ExitCodes PrintStatistics(unsigned int NumberOfCells = 10);
            void SetMesh(GenericMesh& mesh) {meshPointer = &mesh;}
            const Output::ExitCodes InitializeIdCells(const unsigned int& numberOfCells) { idCellToRefine.reserve(numberOfCells); return Output::Success;}
            const Output::ExitCodes	AddIdCell(const unsigned int& idCell);
            const Output::ExitCodes FindQualityCell(GenericCell& cell, double& quality, bool& check_problem);
            inline const Output::ExitCodes FindQualityCell_2(GenericCell& cell, double& quality);
            const Output::ExitCodes FindWorstQualityCells(multimap<double, unsigned int>& WorstQualityCells, unsigned int &num_Cells, double& mean);
            const Output::ExitCodes RefineMesh();
    };
}

#endif // REFINERTETRA_H
