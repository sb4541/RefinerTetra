#ifndef MESHIMPORT_TETGEN_H
#define MESHIMPORT_TETGEN_H

#include "GenericMesh.hpp"
#include "GenericDomain.hpp"
#include "Eigen/Eigen"
#include "tetgen.h"

using namespace std;
using namespace Eigen;
using namespace MainApplication;

namespace GeDiM
{
	class MeshImport_Tetgen;

	class MeshImport_Tetgen : public GenericMeshImportInterface
	{
		protected:
			tetgenio* inputMeshPointer;
			tetgenio* outputMeshPointer;

			string tetgenOptions;

			MatrixXd constrainedPoints;
			vector< vector< int> > constrainedFacets;

		public:
			MeshImport_Tetgen();
			virtual ~MeshImport_Tetgen();

			// http://wias-berlin.de/software/tetgen/files/tetcall.cxx
			void SetInputMesh(tetgenio* _inputMeshPointer) { inputMeshPointer = _inputMeshPointer; }
			void SetOutputMesh(tetgenio* _outputMeshPointer) { outputMeshPointer = _outputMeshPointer; }
			void SetTetgenOptions(const string& _tetgenOptions) { tetgenOptions = _tetgenOptions; }
			void SetConstrainedPoints(const MatrixXd& points) {constrainedPoints = points;}
			void AllocateConstrainedFacets(const unsigned int& numFacets) {constrainedFacets.resize(numFacets);}
			void SetConstrainedFacets(const vector<int>& idPointsFacet, const unsigned int& position) {constrainedFacets[position] = idPointsFacet;}

			Output::ExitCodes CreateTetgenInput(const GenericDomain& domain);
			Output::ExitCodes CreateTetgenOutput(const GenericDomain& domain);
			virtual Output::ExitCodes CreateMesh(const GenericDomain& domain, GenericMesh& mesh) const;

			Output::ExitCodes ExportTetgenMesh(const string& nameFolder, const string& nameFile) const;
	};
}

#endif // MESHIMPORT_TETGEN_H
