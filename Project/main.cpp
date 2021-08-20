#include "GenericDomain.hpp"
#include "GenericMesh.hpp"
#include "Domain3D.hpp"
#include "MeshImport_Tetgen.hpp"
#include "RefinerTetra.hpp"
#include <time.h>
#include <iostream>

using namespace GeDiM;
using namespace Eigen;

int main(int argc, char** argv)
{
	/// CREATE DOMAIN
	unsigned int id = 1;

	clock_t start =clock();

    Parallelepiped domain(id);
    Vector3d origin(1.0, 1.0, 1.0), length(1.0, 0.0, 0.0), height(0.0, 0.0, 1.0), width(0.0, 1.0, 0.0);

    domain.BuildParallelepiped(origin, length, height, width);

    MeshImport_Tetgen meshImportTetgen;
    meshImportTetgen.SetMinimumNumberOfCells(100);
    GenericMesh mesh;
    meshImportTetgen.CreateTetgenInput(domain);
    meshImportTetgen.CreateTetgenOutput(domain);
    meshImportTetgen.CreateMesh(domain, mesh);

	/// REFINE MESH
	unsigned int numCellToRefiner = 10;
	RefinerTetra refiner;
	refiner.SetMesh(mesh);
	refiner.InitializeIdCells(numCellToRefiner);
	//refiner.SetVector();

	/*printf("I punti sono:");
	for(int i = 0; i < 4; ++i)
    {
	    printf(" %d ", mesh.Cell(347)->Point(i)->Id());
    }
	printf("\n");
	printf("%d %d %d\n", mesh.Face(492)->Point(0)->Id(), mesh.Face(492)->Point(1)->Id(), mesh.Face(492)->Point(2)->Id());
	GenericCell& cell = *mesh.Cell(20);
	unsigned int ids;
	refiner.FindLongestEdge(cell, ids);
	for(int i = 0; i < 6; ++i) {
        printf("%d: %lf\n", cell.Edge(i)->Id(),
               (cell.Edge(i)->Point(0)->Coordinates() - cell.Edge(i)->Point(1)->Coordinates()).squaredNorm());
    }
	printf("%d\n", ids);*/

	int refiner_rounds = 2;
	for(int round = 0; round < refiner_rounds; ++round)
    {
        multimap<double, unsigned int> WorstQualityCells;
        double mean = 0.0;
        refiner.FindWorstQualityCells(WorstQualityCells, numCellToRefiner, mean);
        std::multimap<double, unsigned int>::iterator it = WorstQualityCells.begin();
        for(; it != WorstQualityCells.end(); ++it)
        {
            refiner.AddIdCell(it->second);
        }
        refiner.PrintStatistics(numCellToRefiner);
        //Versione standard
        //for(unsigned int numCell = 0; numCell < numCellToRefiner; numCell++)
            //refiner.AddIdCell(numCell);
        refiner.RefineMesh();
        //refiner.PrintStatistics(numCellToRefiner);
    }
    refiner.PrintStatistics(numCellToRefiner);
	/*for(int cel = 0; cel < mesh.NumberOfCells(); ++cel)
    {
	    GenericCell& cell = *mesh.Cell(cel);
        for(int fac = 0; fac < cell.NumberOfFaces(); ++fac)
        {
            GenericFace& face = *mesh.Face(cell.Face(fac)->Id());
            for(int point = 0; point < face.NumberOfPoints(); ++point)
            {
                unsigned int idTemp = face.Point(point)->Id();
                int flag = 0;
                for(int vertex = 0; vertex < cell.NumberOfPoints(); ++vertex)
                {
                    if(idTemp == cell.Point(vertex)->Id())
                        ++flag;
                }
                if(flag != 1)
                {
                    printf("Error\n");
                }
            }
        }
        for(int edg = 0; edg < cell.NumberOfFaces(); ++edg)
        {
            GenericEdge& edge = *mesh.Edge(cell.Edge(edg)->Id());
            for(int point = 0; point < edge.NumberOfPoints(); ++point)
            {
                unsigned int idTemp = edge.Point(point)->Id();
                int flag = 0;
                for(int vertex = 0; vertex < cell.NumberOfPoints(); ++vertex)
                {
                    if(idTemp == cell.Point(vertex)->Id())
                        ++flag;
                }
                if(flag != 1)
                {
                    printf("Error\n");
                }
            }
        }

    }*/

    clock_t stop =clock();
    double time_elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
    printf("Time elapsed to carry out the refinement is %.5f\n", time_elapsed);

	/// OUTPUT MESH TO MATLAB SCRIPT FOR VISUALIZATION
    mesh.CleanInactiveTreeNode();
	ofstream file("plotTetrahedralMesh.m", ofstream::out );
	if(file.is_open())
	{
        file << "nodes = [";
        for(unsigned int i = 0; i < mesh.NumberOfPoints(); i++)
            file << mesh.Point(i)->Coordinates()(0) << "," <<  mesh.Point(i)->Coordinates()(1) << "," <<  mesh.Point(i)->Coordinates()(2) << ";" << endl;
        file << "];" << endl;

        file << "connectivity = [";
        for(unsigned int i = 0; i < mesh.NumberOfCells(); i++)
            {
                file << mesh.Cell(i)->Point(0)->Id()+1 << "," <<  mesh.Cell(i)->Point(1)->Id()+1 << "," << mesh.Cell(i)->Point(2)->Id()+1 << "," << mesh.Cell(i)->Point(3)->Id()+1 << ";" << endl;
            }
        file << "];" << endl;
        file << "tetramesh(connectivity, nodes);" << endl;
        file.close();
	}
	else
        Output::PrintErrorMessage("Unable to open the file", true);

}

