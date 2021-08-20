// Copyright (C) 2016 Vicini Fabio
//
// This file is part_fracture of the dissertation of the author.
//
// This is a free program: you can redistribute it and/or modify
// it under the terms of the author.
//
// This program is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.
//
// Modified by Vicini Fabio 2016
//
// First added:  2016-05-13

#ifndef __CONFIGFILE_H
#define __CONFIGFILE_H

#include <sstream>

#include <vector>
#include <fstream>
#include <climits>
#include <cfloat>

using namespace std;

namespace Configuration
{
	class ConfigFile;
	class InputParameters;
	class CGParameters;
	class GeometricParameters;
	class DiscretizationParameters;
	class ParallelParameters;
	class GraphicParameters;
	class OutputParameters;
	class DebugParameters;
	class FileNameParameters;
	class LogParameters;
	class ModelParameters;

  class ConfigFile
  {
		private:
			static InputParameters inputParam; ///< Application Input parameters
			static CGParameters cGParam; ///< Conjugate Gradient parameters
			static GeometricParameters geomParam; ///< Geometric parameters
			static DiscretizationParameters discParam; ///< Discretization parameters
			static ParallelParameters paralParam; ///< Parallel parameters
			static GraphicParameters graphicParam; ///< Graphic parameters
			static OutputParameters outputParam; ///< Application Output parameters
			static DebugParameters debugParam; ///< Debug parameters
			static FileNameParameters filesParam; ///< Files name parameters
			static LogParameters logParam; ///< Log parameters
			static ModelParameters modelParam; ///< Model parameters
			static string defaultPropertyValue; ///< DefaultString value

		public:
			static const InputParameters& InParam() { return inputParam; }
			static const CGParameters& CG() { return cGParam; }
			static const GeometricParameters& Geom() { return geomParam; }
			static const DiscretizationParameters& Disc() { return discParam; }
			static const ParallelParameters& Paral() { return paralParam; }
			static const GraphicParameters& Graphic() { return graphicParam; }
			static const OutputParameters& OutParam() { return outputParam; }
			static const DebugParameters& DebugParam() { return debugParam; }
			static const FileNameParameters& FilesParam() { return filesParam; }
			static const LogParameters& Log() { return logParam; }
			static const ModelParameters& Model() { return modelParam; }

			static const string& DefaultPropertyValue() { return defaultPropertyValue; }
			static bool IsPropertyDefault(const string& propertyValue) { return propertyValue == defaultPropertyValue; }

			/// Check if the configuration file exists. If not then the file is created in folder.
			static bool ConfigurationFileExists(const string& nameFile, const string& nameFolder = "", const string& description = "");
			/// Split a string in more strings dived by character
			static vector<string> StringSplit(const string& stringToSplit, const char character);
			/// Find Property in Command Arguments. The function return true if the property was found and store the value in property value. If the property does not exist, it is written in the file with default value
			static bool FindStringProperty(const string& nameProperty, const string defaultValue, string& propertyValue, int argc, char** argv);
			/// Find the property on file. The function return true if the property was found and store the value in property value. If the property does not exist, it is written in the file with default value
			static bool ReadStringProperty(const string& nameProperty, const string defaultValue, string& propertyValue, const string& nameFile, const string& nameFolder = "", const string& description = "");
			/// Find Property in Command Arguments. The function return true if the property was found and store the value in property value. If the property does not exist, it is written in the file with default value
			static bool FindIntProperty(const string& nameProperty, const int defaultValue, int& propertyValue, int argc, char** argv);
			/// Find the property on file. The function return true if the property was found and store the value in property value. If the property does not exist, it is written in the file with default value
			static bool ReadIntProperty(const string& nameProperty, const int defaultValue, int& propertyValue, const string& nameFile, const string& nameFolder = "", const string& description = "");
			/// Find Property in Command Arguments. The function return true if the property was found and store the value in property value. If the property does not exist, it is written in the file with default value
			static bool FindIntProperty(const string& nameProperty, const vector<int> defaultValues, vector<int>& propertyValues, int argc, char** argv);
			/// Find the property on file. The function return true if the property was found and store the value in property value. If the property does not exist, it is written in the file with default value
			static bool ReadIntProperty(const string& nameProperty, const vector<int> defaultValues, vector<int>& propertyValues, const string& nameFile, const string& nameFolder = "", const string& description = "");
			/// Find Property in Command Arguments. The function return true if the property was found and store the value in property value. If the property does not exist, it is written in the file with default value
			static bool FindDoubleProperty(const string& nameProperty, const double defaultValue, double& propertyValue, int argc, char** argv);
			/// Find the property on file. The function return true if the property was found and store the value in property value. If the property does not exist, it is written in the file with default value
			static bool ReadDoubleProperty(const string& nameProperty, const double defaultValue, double& propertyValue, const string& nameFile, const string& nameFolder = "", const string& description = "");
			/// Find Property in Command Arguments. The function return true if the property was found and store the value in property value. If the property does not exist, it is written in the file with default value
			static bool FindDoubleProperty(const string& nameProperty, const vector<double> defaultValues, vector<double>& propertyValue, int argc, char** argv);
			/// Find the property on file. The function return true if the property was found and store the value in property value. If the property does not exist, it is written in the file with default value
			static bool ReadDoubleProperty(const string& nameProperty, const vector<double> defaultValues, vector<double>& propertyValue, const string& nameFile, const string& nameFolder = "", const string& description = "");
			/// Find Property in Command Arguments. The function return true if the property was found and store the value in property value. If the property does not exist, it is written in the file with default value
			static bool FindBoolProperty(const string& nameProperty, const bool defaultValue, bool& propertyValue, int argc, char** argv);
			/// Find the property on file. The function return true if the property was found and store the value in property value. If the property does not exist, it is written in the file with default value
			static bool ReadBoolProperty(const string& nameProperty, const bool defaultValue, bool& propertyValue, const string& nameFile, const string& nameFolder = "", const string& description = "");
			/// Check the property boundary. If the property is out of boundary then a message is print and the default value is set
			static bool CheckIntProperty(const string& nameFile, const string& nameProperty, const int defaultValue, int& propertyValue, const int& minValue = -INT_MAX, const int& maxValue = INT_MAX);
			/// Check the property boundary. If the property is out of boundary then a message is print and the default value is set
			static bool CheckIntProperty(const string& nameFile, const string& nameProperty, const vector<int> defaultValues, vector<int>& propertyValues, const int& minValue = -INT_MAX, const int& maxValue = INT_MAX);
			/// Check the property boundary. If the property is out of boundary then a message is print and the default value is set
			static bool CheckDoubleProperty(const string& nameFile, const string& nameProperty, const double defaultValue, double& propertyValue, const double& minValue = -DBL_MAX, const double& maxValue = DBL_MAX);
			/// Check the property boundary. If the property is out of boundary then a message is print and the default value is set
			static bool CheckDoubleProperty(const string& nameFile, const string& nameProperty, const vector<double> defaultValues, vector<double>& propertyValues, const double& minValue = -DBL_MAX, const double& maxValue = DBL_MAX);
			/// Read all the configuration for files
			static void InitializeConfigurations(const string& configurationFolder, int argc, char** argv);
	};

	class InputParameters
	{
		protected:
			string inputFolder; ///< Folder where all inputs are taken
			string inputData; ///< Folder name of input data taken from the ini file
			string runId; ///< Id of the program run
			int computeBoundaryConditions; ///< Used to compute the boundary conditions of DFN using tess fractures in .fab file
			bool computeConnectivity; ///< Compute DFN connectivity
			int discretizeDFN; ///< Discretize the DFN
			bool assembleSystem; ///< Assemble the system. Active only if discretizeDFN is active
			bool computeSolution; ///< Compute the solution of the system. Active only if assembleSystem is active
      bool activateProfiler; ///< Compute time for program profiling

		public:
			InputParameters();

			const string& InputData() const { return inputData; }
			const string InputFolder() const { return inputFolder + "/" + inputData + "/"; }
			const string& RunId() const { return runId; }
			const bool& ComputeConnectivity() const { return computeConnectivity; }
			const int& DiscretizeDFN() const { return discretizeDFN; }
			const bool& AssembleSystem() const { return assembleSystem; }
			const bool& ComputeSolution() const { return computeSolution; }
			const bool& ActivateProfiler() const { return activateProfiler; }
			const int& ComputeBoundaryConditions() const { return computeBoundaryConditions; }

			/// Read configuration from file ini
			void ReadConfiguration(const string& nameFile, const string& nameFolder = "", int argc = 0, char** argv = NULL);
	};
  class CGParameters
  {
		protected:
			double alpha; ///< The parameter used in the flow computation
			double penaltyFactorDirichlet; ///< The parameter used in the flow computation for boundary traces
			int diffusionTermAmplificationType; ///< Used to amplify the effect of the diffusion term
			double domainSize; ///< Size of the domain of the DFN
			double diffusionTermAmplificationConstant; ///< Used to amplify the effect of the diffusion term
			bool traceLengthScaling; ///< Activate the scaling factor by traces length
			int preconditioner; ///< Set Preconditioned CG or not
			bool restart; ///< Restart the approximation of the solution of the CG method
			int refine; ///< Used to refine the solution every \p refine multiple iterations
			vector<int> maxIterations; ///< Max number of iterations for the convergence of CG method
			vector<double> convergence; ///< Min value of the relative residual for the convergence of the CG method

		public:
			CGParameters();

      const double& Alpha() const { return alpha; }
      const double& PenaltyFactorDirichlet() const { return penaltyFactorDirichlet; }
      const double& DomainSize() const { return domainSize; }
      const int& DiffusionTermAmplificationType() const { return diffusionTermAmplificationType; }
      const double& DiffusionTermAmplificationConstant() const { return diffusionTermAmplificationConstant; }
			const bool& TraceLengthScaling() const { return traceLengthScaling; }
			const int& Preconditioner() const { return preconditioner; }
			const bool& Restart() const { return restart; }
			const int& Refine() const { return refine; }
			const int& MaxIterations() const { return maxIterations.back(); }
			const vector<int>& MaxIterationsVector() const { return maxIterations; }
			const double& Convergence() const { return convergence.back(); }
			const vector<double>& ConvergenceVector() const { return convergence; }

			/// Read configuration from file ini
			void ReadConfiguration(const string& nameFile, const string& nameFolder = "", int argc = 0, char** argv = NULL);
  };
  class GeometricParameters
  {
		protected:
			int numEdgesRoundedFractures; ///< Number of edges of the rounded fractures.
			double fractureIntersectionTolerance; ///< Tolerance used for trace research
			double fractureRotationTolerance; ///< Tolerance used for 2D points on fracture rotation
			double fractureCellTolerance; ///< Tolerance used to search points inside a fracture cell
			double traceFractureEdgeTolerance; ///< Tolerance used to search trace intersection with fracture edge
			double boundaryTraceFractureEdgeTolerance; ///< Tolerance used to search boundary trace intersection with fracture edge
			double pointCorrespondenceTolerance; ///< Tolerance used for points identification
			double lineParallelismTolerance; ///< Tolerance used for straight lines parallelism
			double smoothingTolerance; ///< Tolerance used for the smoothing process

		public:
			GeometricParameters();

			const int& NumEdgesRoundFractures() const { return numEdgesRoundedFractures; }
      const double& FractureIntersectionTolerance() const { return fractureIntersectionTolerance; }
      const double& FractureRotationTolerance() const { return fractureRotationTolerance; }
      const double& FractureCellTolerance() const { return fractureCellTolerance; }
      const double& TraceFractureEdgeTolerance() const { return traceFractureEdgeTolerance; }
			const double& BoundaryTraceFractureEdgeTolerance() const { return boundaryTraceFractureEdgeTolerance; }
      const double& PointCorrespondenceTolerance() const { return pointCorrespondenceTolerance; }
      const double& LineParallelismTolerance() const { return lineParallelismTolerance; }
      const double& SmoothingTolerance() const { return smoothingTolerance; }

			/// Read configuration from file ini
			void ReadConfiguration(const string& nameFile, const string& nameFolder = "", int argc = 0, char** argv = NULL);
  };
  class DiscretizationParameters
  {
		protected:
      int fractureMethodType; ///< Method type used for the fracture discretization
			int traceMethodType; ///< Method type used for the trace discretization
			int boundaryTraceMethodType; ///< Method type used for the boundary trace discretization
			int	removeTraceEndPoints; ///< Remove end points from the discretization on the trace
			int fracturePolynomialDegree; ///< Number of degrees of freedom on the fracture method
			int fractureQuadratureOrder; ///< Order of the quadrature formula on fractures
			vector<double> fractureMeshParameter; ///< Mesh parameter of fracture mesh
			int fractureMeshMinNumberElement; ///< Minimum number of element of fracture mesh
			int fractureMeshMaxNumberElement; ///< Maximum number of element of fracture mesh
			bool traceContinuosElements; ///< Activate the continous elements
			int tracePolynomialDegree; ///< Number of degrees of freedom on the trace method
			int traceQuadratureOrder; ///< Order of the quadrature formula on traces
			double traceMeshStep; ///< Mesh step of trace mesh. Use only if traceMethodType is higher than 1
			double traceMeshAmplificationConstant; ///< Mesh amplification constant applied to traceMeshStep. Use only if traceMethodType is higher than 2
			double boundaryTraceMeshStep; ///< Mesh step of boundary trace mesh. Use only if the boundaryTraceMethodType is 1
			int boundaryTracePolynomialDegree; ///< Number of degrees of freedom on the boundary trace method
			int boundaryTraceQuadratureOrder; ///< Order of the quadrature formula on boundary traces
			bool smoothFractureMesh; ///< Smooth the fractures discretization
			int smoothFractureMeshType; ///< Method type used to do the smoothing of fracture mesh
			double fractureMeshRefinementThreshold; ///< Activate the refining process
			int fractureMeshRefinementMaxStep; ///< Max number of steps of refining process
			int fractureMeshRefinementThresholdType; ///< Refining process threshold computation type
			int fractureMeshRefinementAreaDivision; ///< Number of division of area cell in refinement

		public:
			DiscretizationParameters();

      const int& FractureMethodType() const { return fractureMethodType; }
      const int& TraceMethodType() const { return traceMethodType; }
      const int& BoundaryTraceMethodType() const { return boundaryTraceMethodType; }
      const int& RemoveTraceEndPoints() const { return removeTraceEndPoints; }
      const int& FracturePolynomialDegree() const { return fracturePolynomialDegree; }
      const int& FractureQuadratureOrder() const { return fractureQuadratureOrder; }
			const bool& TraceContinuosElements() const { return traceContinuosElements; }
      const int& TracePolynomialDegree() const { return tracePolynomialDegree; }
      const int& TraceQuadratureOrder() const { return traceQuadratureOrder; }
      const vector<double>& FractureMeshParameter() const { return fractureMeshParameter; }
      const int& FractureMeshMinNumberElement() const { return fractureMeshMinNumberElement; }
      const int& FractureMeshMaxNumberElement() const { return fractureMeshMaxNumberElement; }
      const double& TraceMeshStep() const { return traceMeshStep; }
      const double& TraceMeshAmplificationConstant() const { return traceMeshAmplificationConstant; }
      const double& BoundaryTraceMeshStep() const { return boundaryTraceMeshStep; }
      const int& BoundaryTraceQuadratureOrder() const { return boundaryTraceQuadratureOrder; }
      const int& BoundaryTracePolynomialDegree() const { return boundaryTracePolynomialDegree; }
			const bool& SmoothFractureMesh() const { return smoothFractureMesh; }
			const int& SmoothFractureMeshType() const { return smoothFractureMeshType; }
			const double& FractureMeshRefinementThreshold() const { return fractureMeshRefinementThreshold; }
			const int& FractureMeshRefinementMaxStep() const { return fractureMeshRefinementMaxStep; }
			const int& FractureMeshRefinementThresholdType() const { return fractureMeshRefinementThresholdType; }
			const int& FractureMeshRefinementAreaDivision() const { return fractureMeshRefinementAreaDivision; }

			/// Read configuration from file ini
			void ReadConfiguration(const string& nameFile, const string& nameFolder = "", int argc = 0, char** argv = NULL);
  };
  class ParallelParameters
  {
		protected:
      int dfnConnectivityDivisionType; ///< Type of division of fractures among processes on DFN Connectivity computation
			double masterFracturePartition; ///< Partition weight on Master for division, from 0.0 to 1.0
			int partitionType; ///< Partition type
			int metisDivisionType; ///< METIS division type. Corresponds to METIS_OPTION_OBJTYPE
			int cudaStreamNumber; ///< Number of the Streams in GPU
			int cudaVersionSolver; ///< Version of the triangular solver in GPU

		public:
			ParallelParameters();

      const int& DfnConnectivityDivisionType() const { return dfnConnectivityDivisionType; }
      const double& MasterFracturePartition() const { return masterFracturePartition; }
			const int& PartitionType() const { return partitionType; }
			const int& MetisDivisionType() const { return metisDivisionType; }
			const int& CudaStreamNumber() const  { return cudaStreamNumber; }
			const int& CudaVersionSolver() const  { return cudaVersionSolver; }

			/// Read configuration from file ini
			void ReadConfiguration(const string& nameFile, const string& nameFolder = "", int argc = 0, char** argv = NULL);
  };
  class GraphicParameters
  {
		protected:
      bool plotGlobalDFN; ///< Plot Global DFN
      bool plotLocalDFN; ///< Plot Local DFN on process
      bool plotFractures; ///< Plot Fractures
      vector<int> fracturesToPlot; ///< Filter on Fractures to plot
			bool plotTraces; ///< Plot Trace
			vector<int> tracesToPlot; ///< Filter on Traces to plot
			bool plotBoundaryTraces; ///< Plot Trace
			vector<int> boundaryTracesToPlot; ///< Filter on Traces to plot
			double scaleFactor; ///< Scale factor of the output images
			vector<double> rotationGraphic; ///< Angle values for graphics rotation

		public:
			GraphicParameters();

      const bool& PlotGlobalDFN() const { return plotGlobalDFN; }
      const bool& PlotLocalDFN() const { return plotLocalDFN; }
			const bool& PlotFractures() const { return plotFractures; }
			const vector<int>& FracturesToPlot() const { return fracturesToPlot; }
			const bool& PlotTraces() const { return plotTraces; }
			const vector<int>& TracesToPlot() const { return tracesToPlot; }
			const bool& PlotBoundaryTraces() const { return plotBoundaryTraces; }
			const vector<int>& BoundaryTracesToPlot() const { return boundaryTracesToPlot; }
			const double& ScaleFactor() const { return scaleFactor; }
			const vector<double>& RotationGraphic() const { return rotationGraphic; }

			/// Read configuration from file ini
			void ReadConfiguration(const string& nameFile, const string& nameFolder = "", int argc = 0, char** argv = NULL);
  };
	class OutputParameters
  {
		protected:
			int printSolution; ///< Print the solution
			vector<int> fracturesPrintSolution; ///< Filter on Fractures to print solution
			int exportFormat; ///< Format of the output. If 0 then no output is made
			int exportAdjacencyMatrix; ///< Export the Adjacency Matrix in folder ./Export
			bool exportConnectedNetwork; ///< Export the connected network into folder ./Export
			bool exportMesh; ///< Export mesh into folder ./Export
			bool exportModel; ///< Export the model matrices into folder ./Export
			int exportToParaview; ///< Activate export to Paraview
			int exportToParaviewFormat; ///< Format of the Paraview export
			int exportToParaviewEachIteration; ///< Export to paraview each iteration (Active in non stationary mode)
			int postProcessSolutionMarker; ///< Marker used to identify post process solution
			bool exportCudaSystem; ///< Export Matrices in a binary file

		public:
			OutputParameters();

			const int& PrintSolution() const { return printSolution; }
			const vector<int>& FracturesPrintSolution() const { return fracturesPrintSolution; }
			const int& ExportFormat() const { return exportFormat; }
			const int& ExportAdjacencyMatrix() const { return exportAdjacencyMatrix; }
			const bool& ExportConnectedNetwork() const { return exportConnectedNetwork; }
			const bool& ExportMesh() const { return exportMesh; }
			const bool& ExportModel() const { return exportModel; }
			const int& ExportToParaview() const { return exportToParaview; }
			const int& ExportToParaviewFormat() const { return exportToParaviewFormat; }
			const int& ExportToParaviewEachIteration() const { return exportToParaviewEachIteration; }
			const int& PostProcessSolutionMarker() const { return postProcessSolutionMarker; }
			const bool& ExportCudaSystem() const { return exportCudaSystem; }

			/// Read configuration from file ini
			void ReadConfiguration(const string& nameFile, const string& nameFolder = "", int argc = 0, char** argv = NULL);
  };
	class DebugParameters
	{
		protected:
			int startingIndexToPrint; ///< Index where to start to print vectors, lists and so on
			int maxNumElementToPrint; ///< Max number of elements to print for vectors, lists and so on
      bool printApplicationInfo; ///< Print file with application general information
      bool printProcessInfo; ///< Print file with process general information
      vector<int> processesToPrint; ///< Filter on process to plot
      bool printFractureInfo; ///< Print file with fracture general information
      vector<int> fracturesToPrint; ///< Filter on fracture to plot
			bool printTraceInfo; ///< Print file with trace general information
			vector<int> tracesToPrint; ///< Filter on trace to plot
			bool printBoundTraceInfo; ///< Print file with boundary trace general information
			vector<int> boundTracesToPrint; ///< Filter on boundary trace to print
			bool mantainAlteredFiles; ///< Mantain DFN and BoundaryConditions altered files
			int printEachIteration; ///< Print solver iterations every printEachIteration

		public:
			DebugParameters();

			const int& StartingIndexToPrint() const { return startingIndexToPrint; }
			const int& MaxNumElementToPrint() const { return maxNumElementToPrint; }
			const bool& PrintApplicationInfo() const { return printApplicationInfo; }
			const bool& PrintProcessInfo() const { return printProcessInfo; }
			const vector<int>& ProcessesToPrint() const { return processesToPrint; }
			const bool& PrintFractureInfo() const { return printFractureInfo; }
			const vector<int>& FracturesToPrint() const { return fracturesToPrint; }
			const bool& PrintTraceInfo() const { return printTraceInfo; }
			const vector<int>& TracesToPrint() const { return tracesToPrint; }
			const bool& PrintBoundTraceInfo() const { return printBoundTraceInfo; }
			const vector<int>& BoundTracesToPrint() const { return boundTracesToPrint; }
			const bool& MantainAlteredFiles() const { return mantainAlteredFiles; }
			const int& PrintEachIteration() const { return printEachIteration; }

			/// Read configuration from file ini
			void ReadConfiguration(const string& nameFile, const string& nameFolder = "", int argc = 0, char** argv = NULL);
	};
	class FileNameParameters
	{
		protected:
			int dfnInputFormat; ///< Format for DFN file
      string dfnFileName; ///< Input DFN file name
      string outputFolder; ///< Output folder name
      string outputData; ///< Output folder data name
      string discretizationFileName; ///< Input discretization condition file name
      string boundaryConditionFileName; ///< Input boundary condition file name
      string dirichletFileName; ///< Input Dirichlet file name
      string neumannFileName; ///< Input Neumann file name
      string forcingTermFileName; ///< Input forcing term file name
      string diffusionTermFileName; ///< Input diffusion term file name
      string exactSolutionFileName; ///< Input exact solution file name
			string exactSolutionOnTraceFileName; ///< Input exact solution on Trace file name
      string dfnStructFile; ///< Input File containing the DFN structure
      string dfnDofFile; ///< Input File containing the DFN dofs
      string dfnFracFile; ///< Input File containing the fracture informations
			string triangleMeshOptions; ///< Option string for Triangle mesh generator
			string debugFolder; ///< Folder name where debug outputs are put
      string graphicFolder; ///< Folder name where graphic outputs are put
      string processFolder; ///< Folder name where process data are put
			string fracturesFolder; ///< Folder name where fractures data are put
			string tracesFolder; ///< Folder name where traces data are put
			string boundTracesFolder; ///< Folder name where boundary traces data are put
			string solutionsFolder; ///< Folder name where solution files of the program are put
			string errorsFile; ///< File name of Solution Summary
			string iterationFolder; ///< Folder name where solution outputs for iterations are put
			string residualFolder; ///< Folder name where solution outputs for residuals are put
			string finalSolutionFolder; ///< Folder name where last solution outputs are put
			string timeFolder; ///< Folder name where time data are put
			string cGLogFile; ///< File name of the CG output
			string jFile; ///< File with the values of the functional J in solutions
			string resFile; ///< File with the values of the residuals in solutions
			string exportFolder; ///< Folder name where export data are put
			string exportMeshFile; ///< File name for the mesh export
			string importFolder; ///< Folder name where external input data are taken

		public:
			FileNameParameters();

			const int& DfnInputFormat() const { return dfnInputFormat; }
			const string& DfnFileName() const { return dfnFileName; }
			const string OutputFolder() const { return outputFolder + "/"; }
			const string& OutputData() const { return outputData; }
			const string& DfnStructFile() const { return dfnStructFile; }
			const string& DfnDofFile() const { return dfnDofFile; }
			const string DfnFracFile(const int& id) const { ostringstream nameFolder;  nameFolder << dfnFracFile << id; return nameFolder.str(); }
			const string& DiscretizationFileName() const { return discretizationFileName; }
			const string& BoundaryConditionFileName() const { return boundaryConditionFileName; }
			const string& DirichletFileName() const { return dirichletFileName; }
			const string& NeumannFileName() const { return neumannFileName; }
			const string& ForcingTermFileName() const { return forcingTermFileName; }
			const string& DiffusionTermFileName() const { return diffusionTermFileName; }
			const string& ExactSolutionFileName() const { return exactSolutionFileName; }
			const string& ExactSolutionOnTraceFileName() const { return exactSolutionOnTraceFileName; }
			const string& TriangleMeshOptions() const { return triangleMeshOptions; }
			const string DebugFolder() const { return debugFolder + "/"; }
      const string GraphicFolder() const { return graphicFolder + "/"; }
      const string ProcessFolder(const int& rank) const { ostringstream nameFolder;  nameFolder << processFolder << rank << "/"; return nameFolder.str(); }
			const string FracturesFolder() const { return fracturesFolder + "/"; }
			const string TracesFolder() const { return tracesFolder + "/"; }
			const string BoundaryTracesFolder() const { return boundTracesFolder + "/"; }
			const string SolutionsFolder() const { return solutionsFolder + "/"; }
			const string& ErrorsFile() const { return errorsFile; }
			const string IterationFolder(const int& iterationIndex) const { ostringstream nameFolder;  nameFolder << iterationFolder << iterationIndex << "/"; return nameFolder.str(); }
			const string ResidualFolder(const int& residualIndex) const { ostringstream nameFolder;  nameFolder << residualFolder << residualIndex << "/"; return nameFolder.str(); }
			const string FinalSolutionFolder() const { return finalSolutionFolder + "/"; }
			const string TimeFolder() const { return timeFolder + "/"; }
			const string& CGLogFile() const { return cGLogFile; }
			const string& JFile() const { return jFile; }
			const string& ResFile() const { return resFile; }
			const string ExportFolder() const { return exportFolder + "/"; }
			const string& ExportMeshFile() const { return exportMeshFile; }
			const string ImportFolder() const { return importFolder + "/"; }

			/// Read configuration from file ini
			void ReadConfiguration(const string& nameFile, const string& nameFolder = "", int argc = 0, char** argv = NULL);
	};
	class LogParameters
	{
		protected:
			string logFile; ///< File name of Log
			bool createAnalysisFile; ///< Save Analysis file
			string analysisFile; ///< File name of Summary
			int logMaxFileSize; ///< Max file log size (MB)
      int logMaxNumFiles; ///< Max number of backup log files to maintain
      string timeFile; ///< File where time profiling data are stored
			string memoryFile; ///< File where memory profiling data are stored
			string logFolder; ///< Folder name where log files of the program are put

		public:
			LogParameters();

			const string& LogFile() const { return logFile; }
			const bool& CreateAnalysisFile() const { return createAnalysisFile; }
			const string& AnalysisFile() const { return analysisFile; }
			const int& LogMaxFileSize() const { return logMaxFileSize; }
			const int& LogMaxNumFiles() const { return logMaxNumFiles; }
			const string TimeFile() const { return timeFile; }
			const string MemoryFile() const { return memoryFile; }
			const string LogFolder() const { return logFolder + "/"; }

			/// Read configuration from file ini
			void ReadConfiguration(const string& nameFile, const string& nameFolder = "", int argc = 0, char** argv = NULL);
	};
	class ModelParameters
	{
		protected:
			double maxTime; ///< Maximum time for non stationary problem
			double timeStep; ///< Time step used for the non stationary problem discretization
			bool steadyState; ///< If true computes the steady state of the time evolving problem

		public:
			ModelParameters();

			const double& MaxTime() const { return maxTime; }
			const double& TimeStep() const { return timeStep; }
			const bool& SteadyState() const { return steadyState; }

			/// Read configuration from file ini
			void ReadConfiguration(const string& nameFile, const string& nameFolder = "", int argc = 0, char** argv = NULL);
	};
}

#endif
