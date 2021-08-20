#ifndef __MACROPARAMETERS_H
#define __MACROPARAMETERS_H

/// @file macroDefinitions.cpp
/// @brief Used to set special Macro such as error code for the program or Particular Functions
/// @warning If the file is modified then the program \e must be re-compiled with the command \c make \c recompile


/// Use MPI - Supported only for the new Solver
/// - 0 false
/// - 1 true
#define USE_MPI 0

/// Use PETSC - Supported only for the new Solver
/// - 0 false
/// - 1 true
#define USE_PETSC 0

/// Verbose Levels
/// - 0 No output
/// - 1 Only command
/// - 2 Only Files
/// - 3 Command and Files
#define VERBOSE 1

/// Warning Outputs
/// - 0 false
/// - 1 true
#define CONTROL_OUTPUT 0

/// Graphic with mathgl library
/// - 0 Graphic not enabled
/// - 1 Graphic enabled - png files are printed in \c Graphic folder
#define ENABLE_GRAPHIC 0

/// Use Paraview library with XDMF & HDF5 export
/// - 0 Paraview not enabled
/// - 1 Paraview enabled
#define ENABLE_PARAVIEW 0

/// CUDA DEBUG
/// - 0 CUDA DEBUG not enabled
/// - 1 CUDA DEBUG enabled
/// - 2 CUDA save some info about matrices into a file
#define DEBUG_CUDA_SOLVER 0

/// Solver to use to solve Fracture System
/// - 0 Eigen with separate fracture systems
/// - 1 Eigen with agglomerate fracture systems
/// - 2 GPU (CUDA)
#define FRACTURE_SYSTEMS_SOLVER 1

/// SOLVER FORMAT (compatible only with Solver 5)
/// - 0 Solver 5 with CSR Format
/// - 1 Solver 5 with HYB Format
#define CUDA_SOLVER_FORMAT 0

/// MATRIX FORMAT
/// - 0 Matrices with CSR Format
/// - 1 Matrices with HYB Format
#define CUDA_MATRIX_FORMAT 0

/// Non Stationary Problem
#define NON_STATIONARY 0

/// Use VTK library
/// - 0 VTK not enabled
/// - 1 VTK enabled
#define ENABLE_VTK 0

/// Use METIS library
/// - 0 METIS not enabled
/// - 1 METIS enabled
#define ENABLE_METIS 0


/// @name Code Simplifications
///@{
#ifndef MIN
#define MIN(a,b) (a < b) ? a : b
#endif

#ifndef MAX
#define MAX(a,b) (a > b) ? a : b
#endif
///@}

#endif
