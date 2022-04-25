// ##################################################
// #   main.c - finally revised on Apr 2022         #
// #   coded by Taeyoon Kim                         #
// #   Copyright (C) 2005 - 2022, Taeyoon Kim,      #
// #   All rights reserved.                         #
// ##################################################
// This file contains the main process which calls many functions to perform
// a simulation.

#include "common.c"

int main (int argc, char *argv[]) { 
  // The general initialization for MP
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nCpu);

  // Check whether 'condition' and 'parallel' exist.
  InitFileCheck();
  // Initialize the number of cells in x,y,z directions based on 'parallel'
  InitCellNumber();
  // Load initial parameters from two files, "condition" and "Config"
  LoadInitParameter(); 
  // Initialize a run - initialize the random-number generator, 
  // assign values to variables, define arrays,and assign values to the arrays
  InitRun();    
  PrepareStateWoNetworkData();
  // If necessary, add free actins and ABPs
  AddFreeActin();
  AddFreeAbp();
  // Update the neighboring list
  UpdateNeighborList();
  // Record the loaded parameter values
  if (rank == 0) { RecordInitParameter(); }
  // Main process
  MainProcess();
  return 0;
}
