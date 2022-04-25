// ##################################################
// #   process.c - finally revised on Apr 2022      #
// #   coded by Taeyoon Kim                         #
// #   Copyright (C) 2005 - 2022, Taeyoon Kim,      #
// #   All rights reserved.                         #
// ##################################################
// This file contains detailed processes.

void InitSingleProcessSubroutine(double *arr, int cnt) {
  double *ip;

  ip = arr;
  while (ip < arr + cnt) {
    *ip = 0.;
    ip++;
  }
}

// Initiate each single step 
void InitSingleProcess(void) {
  currTimeStep++;
  // Forces are set to zero before SingleProcess().
  // forces of actins.
  InitSingleProcessSubroutine(act.f, (nActMe + nActCp) * NDIM);
  // Brownian forces of actins.
  InitSingleProcessSubroutine(act.fBr, (nActMe + nActCp) * NDIM);
  // forces of ACPs and motors.
  InitSingleProcessSubroutine(abp.f, (nAbpMe + nAbpCp) * NDIM);
  // Brownian forces of ACPs and motors.
  InitSingleProcessSubroutine(abp.fBr, (nAbpMe + nAbpCp) * NDIM);
  if (tglActFormDyn != 0) {
	cntNucAss = 0;
  }

  if (currTimeStep == netForm.dur && rheo.dur > 0) {
	Printf0("=============================== Main measurement "
			"===============================\n");
  }
}

void SingleProcess(void) {

  // Initialize each step
  InitSingleProcess ();
  // Adjust size of subdomains depending on the number of particles in each 
  // Measure displacements of all particles to determine whether
  // neighboring list should be updated or not at a current time step.
  MeasureDisplacementForNL();
  // Repulsive forces between actin segments
  if (actF.rep.facStf > 0. || abpF.rep.facStf > 0.) {
	CalcRepulsiveForces(); 
  }
  // Spring forces for actin, ACP, and motor
  CalcSpringForces();
  // Bending forces for actin filaments
  if (actF.bend.facStf > 0.) {
	CalcFilaBendForces();
  }
  // Brownian forces for actin, ACP, and motor
  CalcBrownForces();
  // Two kinds of bending forces for ACPs
  if (!(motSA.gTgl != 0 && nAbp - nMot == 0)) {
	CalcAbpBendForces();
  }
  // Spring and bending forces acting on the backbone of multimerized motors
  if (motSA.gTgl != 0 && nMot > 0) { 
	CalcMotorBackboneForces();
  }

  /*--------------------------- Dynamics of actins ---------------------------*/
  if (rho.tgl != 0) {
	if (currTimeStep >= netForm.dur) {
		UpdateRhoActivation();
	}
  }
  // Polymerization of actin monomers on filaments
  if (actAss.tgl != 0 && actM.c > 0) {
	UpdateActinAssembly();
  }
  // Depolymerization of actin monomers on filaments
  if (actDis.tgl != 0) {
	UpdateActinDisassembly();
  }
  // Nucleation of actin
  if (actNuc.tgl != 0 && actM.c > 1) {
	UpdateActinNucleation();
  }
  // Eliminate expired elements in noActDyn.l
  if (tglActMoDyn != 0) {
	UpdateNoActinDynamicsList();
  }
  // If actin disassembly/severing is allowed with bound ABPs, this function 
  // should be executed here to let other adjacent subdomains know about ABP 
  // unbinding events caused by disassembly/severing.
  if (nCpu > 1 && actDis.tgl != 0) {
	UpdateActinAbpDynamicsEvents(0); 
  }
  /*----------------------- Dynamics of ACPs or motors -----------------------*/
  // Binding event of monomeric ACP or motor
  if ((acpMoBind.tgl != 0 && acpM.c > 0 && gTglImpAcpM != 0) 
		|| (motMoBind.tgl != 0 && motM.c > 0 && gTglImpMotM != 0)) {
	UpdateAbpMonomerBinding();
  }
  if (motSA.gTgl != 0) {
	if (motSA.to.gTgl != 0) {
		UpdateMotorTurnover();
	}
	if (currTimeStep > netForm.dur) {
		if ((motM.c > 0 && gTglImpMotM != 0) || gTglImpMotM == 0) {
			UpdateMotorAssembly();
		}
	}
  }
  // Unbinding event of inactive ACP or motor to monomer
  if (acpInaUnb.tgl != 0 || motInaUnb.tgl != 0) {
	UpdateInactAbpUnbinding();
  }
  // Binding event of inactive ACP or motor
  if (acpReb.tgl != 0 || motReb.tgl != 0) {
	UpdateAbpBinding();
  }
  // Unbinding event of active ACP or motor
  if (acpUnb.tgl != 0 || motUnb.tgl != 0) { 
	UpdateActiveAbpUnbinding();
  }
  // Walking event of active motor
  if (motWalk.tgl != 0) {
	UpdateMotorWalking();
  }
  // Eliminate expired elements in noAbpDyn.l 
  if (tglAbpAcInaDyn != 0 || acpMoBind.tgl != 0 || motMoBind.tgl != 0 
		|| (actDis.tgl != 0 && actDis.facKWA > 0.)) { 
	UpdateNoAbpUnbindList();
  }
  // transfer actin severing/disassembly and ABP unbinding/binding events 
  // stored in sendAbpDyn.l to adjacent CPUs
  if (nCpu > 1 && tglAbpAcInaDyn != 0) { 
	UpdateActinAbpDynamicsEvents(1);   
  }
  /*--------------------------------------------------------------------------*/
  // Calculate new positions of all particles based on forces 
  UpdateNewLocation();
  // Apply the periodic boundary condition to all particles
  ApplyBoundCondAll();
  // Only for a case with more than 1 CPU
  if (nCpu > 1) { 
	// Process longCh.l in either normal or Plympton way
	UpdateLongChainNormal(); 
	// Move particles escaping from a current subdomain to adjacent subdomains
	MoveParticles();
	// Copy particles located near boundaries to adjacent subdomains 
	// for synchronization
	CopyParticles();
  }
}

void InitToggleParameters(int mode) {
  // Actin assembly or disassembly is on or off 
  // depending on parameters in 'condition'.
  if (nAct > 0 && (mode == 1 || (mode == 0 && gTglActDynNF != 0))) {
	actAss.tgl = actAss.gTgl;
	actDis.tgl = actDis.gTgl;
	actNuc.tgl = actNuc.gTgl;
  }
  else {
	actAss.tgl = 0;
	actDis.tgl = 0;
	actNuc.tgl = 0;
  }
  // Binding of monomeric ACPs or motors is on or off.
  if ((mode == 1 && nAbp > 0) 
		|| (mode == 0 && ((gTglAcpDynNF != 0 && nAbp - nMot > 0) 
		|| (gTglMotUnbRebNF != 0 && nMot > 0)))) { 
	acpMoBind.tgl = acpMoBind.gTgl;
	motMoBind.tgl = motMoBind.gTgl;
  }
  else { 
	acpMoBind.tgl = 0; 
	motMoBind.tgl = 0; 
  }
  // Unbinding of active ACPs or motors and 
  // unbinding and binding of inactive ACPs of motors are on or off.
  if (nAbp - nMot > 0 && (mode == 1 || (mode == 0 && gTglAcpDynNF != 0))) {
	acpUnb.tgl = acpUnb.gTgl;
	acpReb.tgl = acpReb.gTgl;
	acpInaUnb.tgl = acpInaUnb.gTgl;
  }
  else {
	acpUnb.tgl = 0;
	acpReb.tgl = 0;
	acpInaUnb.tgl = 0;
  }
  // Unbinding of active motors and the unbinding or binding of inactive
  // motors are on or off. 
  if (nMot > 0 && (mode == 1 || (mode == 0 && gTglMotUnbRebNF != 0))) {
	motUnb.tgl = motUnb.gTgl;
	motReb.tgl = motReb.gTgl;
	motInaUnb.tgl = motInaUnb.gTgl;
  }
  else {
	motUnb.tgl = 0;
	motReb.tgl = 0;
	motInaUnb.tgl = 0;
  }
  // Walking of active motors is on or off.
  if (nMot > 0 && (mode == 1 || (mode == 0 && gTglMotWalkNF != 0))) {
	motWalk.tgl = motWalk.gTgl;
  }
  else { 
	motWalk.tgl = 0; 
  }
  // Toggle parameters representing multiple toggles are on or off.
  tglAcpAcInaDyn = (acpUnb.tgl != 0 || acpReb.tgl != 0 
		|| acpInaUnb.tgl != 0) ? 1 : 0;
  tglMotAcInaDyn = (motUnb.tgl != 0 || motReb.tgl != 0 || motWalk.tgl != 0 
		|| motInaUnb.tgl != 0) ? 1 : 0;
  tglAbpAcInaDyn = (tglAcpAcInaDyn != 0 || tglMotAcInaDyn != 0) ? 1 : 0;
  tglAbpInaMoDyn = (acpMoBind.tgl != 0 || motMoBind.tgl != 0
		|| acpInaUnb.tgl != 0 || motInaUnb.tgl != 0) ? 1 : 0;
  tglActFormDyn = (actAss.tgl != 0 || actNuc.tgl != 0) 
		? 1 : 0;
  tglActMoDyn = (actAss.tgl != 0 || actDis.tgl != 0	|| actNuc.tgl != 0) ? 1 : 0;
  tglNeiAbpSC = (acpReb.gTgl != 0 || motReb.gTgl != 0 
		|| abpF.rep.facStf > 0.) ? 1 : 0;
  tglNeiAbpDC = (abpF.rep.facStf > 0.) ? 1 : 0;
}

// Main process
// There are two separate parts for each rheology measurement.
void MainProcess(void) {
  // Synchronization
  MPI_Barrier(MPI_COMM_WORLD);
  // Initial output.
  if (rank == 0) { 
	InitOutput(); 
  }
  InitToggleParameters(0);
  // This loop is iterated until currTimeStep reaches a last time step.
  // "netForm.dur" is the time step required for network formation.
  // "rheo.dur" is the time of rheological measurement in "condition".
  while (currTimeStep <= netForm.dur + rheo.dur) {
	if (currTimeStep == netForm.dur) {
		InitToggleParameters(1);
	}
	// Excute the main part of the process run every time step.
	SingleProcess();
	// Check and record things after SingleProcess().
	TaskAfterSingleProcess();
  }
  Printf0("All process is done!\n"); 
  // Finalize MPI processes
  MPI_Finalize();
}

void TaskAfterSingleProcess(void) {
  /*------------------------------ Record progress ---------------------------*/
  // Record progress information and the mechanical energy of networks.
  if (currTimeStep % recProg.prd == 0 || currTimeStep == 1) {
    RecordProgress(); 
  }
  /*--------------------------- Network information --------------------------*/
  if (recConf.tgl != 0) {
	// Record network information including the chains and positions of 
	// particles.
	if (currTimeStep % recConf.prd == 0 || currTimeStep == 1) {
		RecordConfig(GenFileName("AccuConf"), recConf.prd);
	}
  }
  if (recConfVmd.tgl != 0) {
	// Record networks for visualization via VMD. 
	if (currTimeStep % recConfVmd.prd == 0 || currTimeStep == 1) {
		RecordConfigVmd(0);
	}
  }
  /*--------------------------- Energy calculation ---------------------------*/
  if (recE.tgl != 0) {
	if (currTimeStep % recE.prd == 0 || currTimeStep == 1) {
	    RecordMechEnergy(0);
	}
  }
  /*-------- Information in unit of elements, segments, and filaments --------*/
  if (recInfo.tgl != 0) { 
	if (currTimeStep % recInfo.prd == 0) {
		RecordIndvSegFilaInformation(recInfo.prd);
	}
  }
  if (confVmdInfo.tgl != 0) {
	// If a specific coloring method is used, the accumulated values of chain  
	// lengths and forces have to recorded.
    RecordAccuLengthForces();
	if (currTimeStep % confVmdInfo.prd == 0) {
		ResetAccuLengthForces();
	}
  }
  if (recLongF.tgl != 0) {
	RecordAccuLongSpringForces(recLongF.prd);
  }
  /*------------------ Dynamics behaviors of actin and ABPs ------------------*/

  if (recAbpUnb.tgl != 0) {
	if (currTimeStep % recAbpUnb.prd == 0) {
  		RecordAbpUnbindEvent(-1, -1, 2);
	}
  }
  if (recAbpBind.tgl != 0) {
	if (currTimeStep % recAbpBind.prd == 0) {
  		RecordAbpBindEvent(-1, -1, 2);
	}
  }
  if (recAbpTurn.tgl != 0) {
	if (currTimeStep % recAbpTurn.prd == 0) {
  		RecordAbpTurnover(-1, -1, -1, 2);
	}
  }
  /*-------------------- Structural properties of network --------------------*/
  if (recFilaL.tgl != 0) {
	if (currTimeStep % recFilaL.prd == 0) {
		RecordFilamentLength(recFilaL.prd);
	}
  }

  // If an error signal appears, the run shoulde be terminated.
  if (stopSig != 0) {
	printf("Due to the large force, the code is stopped at currTimeStep = "
		"%lld!\n", currTimeStep);
	while(1);
		exit(-1);
	}
}

void PrepareStateWoNetworkData(void) {
  int n, nIFilaP;

  nAct = 0;
  nActMe = 0;
  nActCp = 0;
  actM.c = 0; 
  nActFilaMe = 0;

  nAbp = 0;
  nAbpMe = 0;
  nAbpCp = 0; 
  if (gTglImpAcpM != 0) { acpM.c = 0; }
  if (gTglImpMotM != 0) { motM.c = 0; }

  neigh.c = 0;
  if (actNuc.gTgl != 0) {
	nIFilaP = (int)(nActGoal / nCpu);
	for(n = 0; n < nIFilaP; n++) {
		iFilaP.l[n] = rank * nIFilaP + n;
	}
	iFilaP.c = nIFilaP;
  }
  UpdateChainList();
  
}
