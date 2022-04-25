// ##################################################
// #   dataMgr.c - finally revised on Apr 2022      #
// #   coded by Taeyoon Kim                         #
// #   Copyright (C) 2005 - 2022, Taeyoon Kim,      #
// #   All rights reserved.                         #
// ##################################################
// This file contains functions performing the modification or loading of 
// data.

/*----------------------- Loading initial parameters -------------------------*/

// Subroutine for LoadInitParameter(), which judges "yes" or "no" from 
// the file, "condition".
void CheckAnswerYesNo(char *tag, const char *which, int *tgl) {
  if (strcmp(tag, "yes") == 0) { *tgl = 1; }
  else if (strcmp(tag, "no") == 0) { *tgl = 0; } 
  else { 
	Printf0("Error: only yes or no is allowed (%s).\n\n", which); 
	exit(-1); 
  }
}

void LoadInitParameterSubroutine2(char *str, int ind, int *tgl1, int *tgl2) {
  *tgl1 = (ind > 0) ? 1 : 0;
  *tgl2 = (ind == 2) ? 1 : 0;
  if (ind < 0 || ind > 2) {
	Printf0("Error: choice of %s should be either of 0, 1, or 2.\n\n", str);
  	MPI_Barrier(MPI_COMM_WORLD);
	exit(-1);
  }
}

void LoadInitParameterSubroutine(FILE *fIn, const char *which, int *tgl) {
  char tag[80];

  fscanf(fIn, "%s", tag);
  if (tag[0] == 'y' && tag[1] == 'e' && tag[2] == 's') 
  { *tgl = 1; }
  else if (tag[0] == 'n' && tag[1] == 'o') 
  { *tgl = 0; }
  else {
	Printf0("Error: only yes or no is allowed (%s).\n\n", which); 
	exit(-1); 
  }
}

// Load parameters from "Config" and "condition" at the very beginning.
// Those parameters are needed to define arrays, and so on.
void LoadInitParameter(void) {
  int n, tempInt[4];
  double tempDbl[3];
  char tag[200], direc[4]= "xyz";
  FILE *fIn;

  nChAb = 5;
  MALLOC(abpF.a90,Force,NK_ABP);
  MALLOC(abpF.cr,Force,NK_ABP);
  MALLOC(abpF.bend,Force,NK_ABP);
  MALLOC(abpF.spr,Force,NK_ABP);

  // From "condition", several types of informatio are loaded: MPI method,
  // rheological methods, dynamic behaviors of ABPs
  if ((fIn = fopen("condition", "r")) == NULL) {
	Printf0("File doesn't exist: condition\n");
	exit(-1); 
  }
  fgets(tag, 200, fIn);
  fscanf(fIn, "======================= Parameters for parallel processing "
		"=====================\n");
  // Coarse-graining using cylindrical segments
  fscanf(fIn, "========== Parameters for coarse-graining using cylindrical "
		"segments ===========\n");
  fscanf(fIn, "Length of cylindrical segments for actins(=i*7nm) = %d\n", 
		&nActPerSeg);
  fscanf(fIn, "Number of binding sites on each segment in longitudinal and "
		"transverse directions = %d, %d\n", &nChAcX, &nChAcY);
  // Network condition
  fscanf(fIn, "======================= Parameters for network condition "
		"=======================\n");
  fscanf(fIn, "Size of domain(x, y, z in um) = %lf, %lf, %lf\n", 
		&tempDbl[0], &tempDbl[1], &tempDbl[2]);
  FOR_NDIM(n) { dimDom[n] = L_UM2S(tempDbl[n]); }
  fscanf(fIn, "Periodic boundary condition(yes/no in x, y, z) = ");
  FOR_NDIM(n) {
	LoadInitParameterSubroutine(fIn, "PBC", &pbc[n]);
  }
  fscanf(fIn, "\n");
  fscanf(fIn, "If 2D network, specifiy the normal direction(no, x, y, or z) "
		"= %s\n", tag);
  if (strcmp(tag, "x") == 0) { dir2D = 0; }
  else if (strcmp(tag, "y") == 0) { dir2D = 1; } 
  else if (strcmp(tag, "z") == 0) { dir2D = 2; } 
  else if (strcmp(tag, "no") == 0) { dir2D = -1; } 
  else {
	Printf0("Error: Normal direction specifying 2D geometry should be "
			"either no, x, y, or z.\n\n"); 
	exit(-1); 
  }
  if (dir2D > -1) {
	if (pbc[dir2D] == 1) {
		Printf0("Error: periodic boundary condition should not exist "
				"in %c direction with the 2D geometry.\n\n", direc[dir2D]); 
		exit(-1); 
	}
  }
  fscanf(fIn, "Duration of network formation(s) = %lf\n", &netForm.durR);

  fscanf(fIn, "Actin concentration(in uM) = %lf\n", &cAct);
  // R values of ACPC, ACPB, and motor
  fscanf(fIn, "ACPC density(R value) = %lf\n", &RAbp[0]);
  fscanf(fIn, "ACPB density(R value) = %lf\n", &RAbp[1]);
  fscanf(fIn, "Motor density(R value) = %lf\n", &RAbp[2]);

  fscanf(fIn, "========================= Parameters for measurement "
		"===========================\n");
  fscanf(fIn, "Duration of simulation(s) = %lf\n", &rheo.durR);

  // Toggle dynamic behaviors of actin and ABP
  fscanf(fIn, "==================== Toggle the dynamic behaviors of actin "
		"=====================\n");
  fscanf(fIn, "Allow the thermal fluctuation of actin filaments(yes/no) "
		"= %s\n", tag);
  CheckAnswerYesNo(tag, "actin thermal fluctuation", &gTglActTherm);

  fscanf(fIn, "Allow actin nucleation/assembly/disassembly during "
		"network formation(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "actin nucleation/assembly/disassembly during network "
		"formation", &gTglActDynNF);

  fscanf(fIn, "Allow nucleation of actins(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "actin nucleation", &actNuc.gTgl);
  fscanf(fIn, "Allow assembly of actins(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "actin assembly", &actAss.gTgl);
  fscanf(fIn, "Allow disassembly of actins(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "actin disassembly", &actDis.gTgl);

  fscanf(fIn, "===================== Toggle the dynamic behaviors of ACP "
		"======================\n");
  fscanf(fIn, "Allow the thermal fluctuation of ACPs(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "ACP thermal fluctuation", &gTglAcpTherm);

  fscanf(fIn, "Allow ACP unbinding/binding during "
		"network formation(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "ACP unbinding/binding during network formation", 
		&gTglAcpDynNF);

  fscanf(fIn, "Allow unbinding of inactive ACPs(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "inactive ACP unbinding", &acpInaUnb.gTgl);
  fscanf(fIn, "Allow binding of monomeric ACPs(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "monomeric ACP binding", &acpMoBind.gTgl);
  fscanf(fIn, "Allow unbinding of active ACPs(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "active ACP unbinding", &acpUnb.gTgl);
  fscanf(fIn, "Allow binding of inactive ACPs(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "inactive ACP binding", &acpReb.gTgl);
  fscanf(fIn, "Allow implicit consideration of monomeric ACPs(yes/no) = %s\n", 
		tag);
  CheckAnswerYesNo(tag, "implicit consideration of monomeric ACPs", 
		&gTglImpAcpM);

  fscanf(fIn, "==================== Toggle the dynamic behaviors of motor "
		"=====================\n");
  fscanf(fIn, "Allow the thermal fluctuation of motors(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "motor thermal fluctuation", &gTglMotTherm);

  fscanf(fIn, "Allow motor unbinding/binding during network formation(yes/no) "
		"= %s\n", tag);
  CheckAnswerYesNo(tag, "motor unbinding/binding during "
		"network formation", &gTglMotUnbRebNF);

  fscanf(fIn, "Allow motor walking during network formation(yes/no) = %s\n",
		tag);
  CheckAnswerYesNo(tag, "motor walking during network formation", 
		&gTglMotWalkNF);

  fscanf(fIn, "Allow walking of motors(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "motor walking", &motWalk.gTgl);
  fscanf(fIn, "Allow self-assembly of motors into filaments(yes/no) = %s\n", 
		tag);
  CheckAnswerYesNo(tag, "motor assembly", &motSA.gTgl);

  fscanf(fIn, "----------------------- If motors are independent dimers "
		"-----------------------\n");
  fscanf(fIn, "Allow unbinding of inactive motors(yes/no) = %s\n", tag);
  if (motSA.gTgl == 0) {
	CheckAnswerYesNo(tag, "inactive motor unbinding", &motInaUnb.gTgl);
  }
  fscanf(fIn, "Allow binding of monomeric motors(yes/no) = %s\n", tag);
  if (motSA.gTgl == 0) {
	CheckAnswerYesNo(tag, "monomeric motor binding", &motMoBind.gTgl);
  }
  fscanf(fIn, "Allow unbinding of active motors(yes/no) = %s\n", tag);
  if (motSA.gTgl == 0) {
	CheckAnswerYesNo(tag, "active motor unbinding", &motUnb.gTgl);
  }
  fscanf(fIn, "Allow binding of inactive motors(yes/no) = %s\n", tag);
  if (motSA.gTgl == 0) {
	CheckAnswerYesNo(tag, "inactive motor binding", &motReb.gTgl);
  }
  fscanf(fIn, "Allow implicit consideration of monomeric motors(yes/no) = %s\n",
		tag);
  CheckAnswerYesNo(tag, "implicit consideration of monomeric motor", 
		&gTglImpMotM);
  if (motSA.gTgl != 0 && gTglImpMotM == 0) {
	gTglImpMotM = 1;
  }

  fscanf(fIn, "------------------------- If motors are self-assembled "
		"-------------------------\n");
  fscanf(fIn, "Allow unbinding of motors(yes/no) = %s\n", tag);
  if (motSA.gTgl != 0) {
	CheckAnswerYesNo(tag, "motor unbinding", &motInaUnb.gTgl);
	motUnb.gTgl = motInaUnb.gTgl;
  }
  fscanf(fIn, "Allow binding of motors(yes/no) = %s\n", tag);
  if (motSA.gTgl != 0) {
	CheckAnswerYesNo(tag, "motor binding", &motMoBind.gTgl);
	motReb.gTgl = motMoBind.gTgl;
  }
  fscanf(fIn, "If yes, allow binding of motors only when actin filaments are "
		"favorably aligned with motor filaments(yes/no) = %s\n", tag);
  if (motSA.gTgl != 0) {
	CheckAnswerYesNo(tag, "actin-motor-filament alignment", &motReb.gTglCrsAng);
  }
  fscanf(fIn, "Allow turnover of motor filaments(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "motor turnover", &motSA.to.gTgl);

  // Parameters for dynamic behaviors of actin and ABP
  fscanf(fIn, "================ Parameters for the dynamic behaviors of "
		"actin =================\n");
  fscanf(fIn, "k for actin nucleation(1/uM s) = %lf\n", &actNuc.k);
  fscanf(fIn, "k for actin assembly at barbed and pointed ends(1/uM s) "
		"= %lf, %lf\n", &actAss.k[0], &actAss.k[1]);
  fscanf(fIn, "k for actin disassembly at barbed and pointed ends(1/s) "
		"= %lf, %lf\n", &actDis.k[0], &actDis.k[1]);
  fscanf(fIn, "Factor for varying disassembly rate of actin with ACPs or "
		"motors = %lf\n", &actDis.facKWA);

  fscanf(fIn, "================= Parameters for the dynamic behaviors of "
		"ACP ==================\n");
  fscanf(fIn, "k for ACP binding(k=i*k_0) = %lf\n", &acpReb.facK);
  acpMoBind.facK = acpReb.facK;

  fscanf(fIn, "k0 and x for ACP slip unbinding = %lf, %lf\n", 
		&acpUnb.facK0, &acpUnb.facX);
  fscanf(fIn, "k0 and x for ACP catch unbinding = %lf, %lf\n",
		&acpUnb.facK0c, &acpUnb.facXc);

  fscanf(fIn, "================ Parameters for the dynamic behaviors of "
		"motor =================\n");
  fscanf(fIn, "k for motor binding(k=i*k_0) = %lf\n", 
		&motReb.facK);
  motMoBind.facK = motReb.facK;
  fscanf(fIn, "k0 for motor unbinding and walking(k0=i*k0_0) = %lf\n", 
		&motUnb.facK0);
  motWalk.facK0 = motUnb.facK0;
  fscanf(fIn, "Number of heads which each motor arm represents = %d\n", 
		&motMC.nHead);
  fscanf(fIn, "Transition rates between mechanochemical states in each head"
		"(k01, k10, k12, k21, k20 in 1/s) = %lf, %lf, %lf, %lf, %lf\n", 
		&motMC.k01, &motMC.k10, &motMC.k12, &motMC.k21, &motMC.k20);
  fscanf(fIn, "------------------------- If motors are self-assembled "
		"-------------------------\n");
  fscanf(fIn, "Average number of motors per each self-assembled structure = "
		"%d\n", &motSA.nMotPerTF);
  fscanf(fIn, "k for motor assembly(in 1/s) = %lf\n", &motSA.kAss);
  fscanf(fIn, "k for motor turnover(in 1/s) = %lf\n", &motSA.to.k);
  fscanf(fIn, "======== Adjustment for rates of the dynamic behaviors of "
		"actin and ABP ========\n");
  fscanf(fIn, "During network formation(i*k) = %lf\n", &netForm.facK);
  fscanf(fIn, "=============== Parameters for the mechanical behaviors of "
		"actin ===============\n");
  fscanf(fIn, "Strength of repulsive forces between actin filaments(Kr=i*Kr_0)"
		" = %lf\n", &actF.rep.facStf);
  fscanf(fIn, "Bending stiffness of actin filaments(Kb=i*Kb_0) = %lf\n", 
		&actF.bend.facStf);
  fscanf(fIn, "Extensional stiffness of actin filaments(Ks=i*Ks_0) = %lf\n", 
		&actF.spr.facStf);

  fscanf(fIn, "================ Parameters for the mechanical behaviors of "
		"ABP ================\n");
  fscanf(fIn, "Strength of repulsive forces between actin and ABP and between "
		"ABPs(Kr=i*Kr_0) = %lf\n", &abpF.rep.facStf);
  fscanf(fIn, "Bending stiffness which maintains an assigned angle between two "
		"arms of ACPC, ACPB, and motor(Kb=i*Kb_0) = %lf, %lf, %lf\n", 
		&abpF.bend[0].facStf, &abpF.bend[1].facStf, &abpF.bend[2].facStf);
  fscanf(fIn, "Bending stiffness which maintains 90 deg between the axis of "
		"a filament and the arm of ACPC, ACPB, and motor(Kb=i*Kb_0) = %lf, "
		"%lf, %lf\n", &abpF.a90[0].facStf, &abpF.a90[1].facStf, 
		&abpF.a90[2].facStf);
  fscanf(fIn, "Extensional stiffness of ACPC, ACPB, and motor(Ks=i*Ks_0) = "
		"%lf, %lf, %lf\n", &abpF.spr[0].facStf, &abpF.spr[1].facStf, 
		&abpF.spr[2].facStf);

  fscanf(fIn, "======================== Parameters for data recording ======"
		"===================\n");
  fscanf(fIn, "Period of recording Output and Progress(s) = %lf\n", 
		&recProg.prdR);

  fscanf(fIn, "----------------------------- Structural information "
		"---------------------------\n");
  fscanf(fIn, "Period of recording Config(in s, -1 for deactivation) = %lf\n",
		&recConf.prdR);
  fscanf(fIn, "Period of recording structural information for visualization "
		"via VMD(in s, -1 for deactivation) = %lf\n", &recConfVmd.prdR);
  fscanf(fIn, "Number of VMD files(multiple/single) = %s\n", tag);
  if (strcmp(tag, "multiple") == 0) { recConfVmd.mode = 0; }
  else if (strcmp(tag, "single") == 0) { recConfVmd.mode = 1; } 
  else {
	Printf0("Error: the number of VMD files should be multiple or single\n");
	exit(-1); 
  }
  fscanf(fIn, "Show the boundaries of a network drawn by VMD(yes/no) = %s\n", 
		tag);
  CheckAnswerYesNo(tag, "show boundaries in VMD", &recConfVmd.gTglBnd);
  fscanf(fIn, "Record information for coloring a network drawn by "
		"VMD(yes/no) = %s\n", tag);
  CheckAnswerYesNo(tag, "information for coloring VMD", &recConfVmd.gTglInfo);
  fscanf(fIn, "Lower and upper limits of forces for coloring(in pN) = %lf, "
		"%lf\n", &recConfVmd.minF, &recConfVmd.maxF);
  fscanf(fIn, "Period of recording length of actin filaments(in s, -1 for "
		"deactivation) = %lf\n", &recFilaL.prdR);

  fscanf(fIn, "--------------------------- Force, stress, and energy "
		"--------------------------\n");
  fscanf(fIn, "Period of recording longitudinal forces of ABPs(in s, -1 for "
		"deactivation) = %lf\n", &recLongF.prdR);
  fscanf(fIn, "Period of recording mechanical energy(in s, -1 for "
		"deactivation) = %lf\n", &recE.prdR);

  fscanf(fIn, "--------------------- Dynamic behaviors of actin and ABPs "
		"----------------------\n");
  fscanf(fIn, "Period of recording ABP unbinding(in s, -1 for deactivation) "
		"= %lf\n", &recAbpUnb.prdR);
  fscanf(fIn, "Period of recording ABP binding(in s, -1 for deactivation) "
		"= %lf\n", &recAbpBind.prdR);
  fscanf(fIn, "Period of recording ABP turnover(in s, -1 for deactivation) "
		"= %lf\n", &recAbpTurn.prdR);

  fscanf(fIn, "---------------------------------- Miscellany "
		"----------------------------------\n");
  fscanf(fIn, "Period of recording information in unit of individual "
		"elements, filament segments, and filaments(in s, -1 for deactivation)"
		" = %lf\n", &recInfo.prdR);

  fscanf(fIn, "========================== Parameters for rho dynamics "
		"=========================\n");
  fscanf(fIn, "Allow rho activation(in s, -1 for deactivation) = %lf\n", 
		&rho.prdR);
  fscanf(fIn, "Number of regions for rho activation in x, y, and z "
		"directions = %d, %d, %d\n", &rho.nReg[0], &rho.nReg[1], &rho.nReg[2]);
  fscanf(fIn, "Number of regions activated by rho at once = %d\n", &rho.nActiv);

  fscanf(fIn, "Delay and duration of motor activation(in s) = "
		"%lf, %lf\n", &rho.dlyMotR, &rho.durMotR);
  fscanf(fIn, "Delay and duration of altered actin polymerization(in s) = "
		"%lf, %lf\n", &rho.dlyActR, &rho.durActR);
  fscanf(fIn, "Density of barbed ends for altered actin polymerization(0-1) "
		"= %lf\n", &rho.RFor);
  fscanf(fIn, "Factor for actin polymerization due to rho activation = %lf, "
		"%lf\n", &rho.facAct[0], &rho.facAct[1]);
  fscanf(fIn, "Portion of activable motors = %lf\n", &rho.fracM);
  fscanf(fIn, "Enhancement of assembly rate of activated motors = %lf\n", 
		&rho.facM);
  fclose(fIn);

  if (cAct >= 0. && actAss.gTgl == 0 && actDis.gTgl == 0) { 
	cAct = -1.; 
  }
  if (actAss.gTgl == 0 && actNuc.gTgl == 0) {
	Printf0("Error: actin assembly and nucleation should be allowed "
			"at least.\n\n"); 
	exit(-1); 
  }

  V3SET_ALL(nAbpDet, 0);

  nAct = 0;
  nAbp = 0;

}

/*----------------------- Loading initial parameters -------------------------*/

/*----------------------- Adding and deleting elements -----------------------*/

void AddFreeActinAbpSubroutine(int nMall, int *nAddMe, int *nAddMeSum) {
  int k, mSum, *mAll;
  double volSum, *volAll, dimDomC[NDIM], volMe;

  V3SUB(dimDomC, &P2A(bnd.r,1,0,NDIM), &P2A(bnd.r,0,0,NDIM));
  volMe = V3PROD(dimDomC);
  
  MALLOC(volAll, double, nCpu);
  MALLOC(mAll, int, nCpu);
  MPI_Gather(&volMe, 1, MPI_DOUBLE, volAll, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	volSum = SumArrDbl(volAll, nCpu);
	for(k = 0; k < nCpu; k++) {
		mAll[k] = (int)(volAll[k] / volSum * nMall);
	}
	mSum = SumArrInt(mAll, nCpu);
	mAll[nCpu - 1] += nMall - mSum;
  }
  MPI_Bcast(mAll, nCpu, MPI_INT, 0, MPI_COMM_WORLD);
  *nAddMe = mAll[rank];
  *nAddMeSum = SumArrInt(mAll, rank);

  free(volAll);
  free(mAll);
}

void AddFreeAbp(void) {
  int n, kind, CS, sum, nAbpAddSumMe, nAbpMallSum;
  int nAbpAddMe[NK_ABP], *nAbpAddSumAll;
  int diff, ind, *nNucAll;
  int cnt;

  Printf0("%d ACPC, %d ACPB, %d motor are added.\n",
		nAbpMall[0], nAbpMall[1], nAbpMall[2]);
  nMot += nAbpMall[2];
  nAbpMallSum = V3SUM(nAbpMall);

  for(n = 0; n < NK_ABP; n++) {
	nAbpAddMe[n] = nAbpMall[n] / nCpu
			+ (rank < (nAbpMall[n] % nCpu) ? 1 : 0);
  }
  if (motSA.gTgl != 0) {
	MALLOC(nNucAll, int, nCpu);
	if (rank == 0) {
		SetAllValue1dArrayInt(nNucAll, nCpu, motSA.nNucMe / nCpu);
		for(n = 0; n < motSA.nNucMe % nCpu; n++) {
			while(1) {
				ind = GenRandIntIndex(nCpu);
				BREAK(nNucAll[ind] == motSA.nNucMe / nCpu);
			}
			nNucAll[ind]++;
		}
		if (SumArrInt(nNucAll, nCpu) == 0 && nAbpMall[2] > 0) {
			ind = GenRandIntIndex(nCpu);
			nNucAll[ind] = 1;
		}
	}
	MPI_Bcast(nNucAll, nCpu, MPI_INT, 0, MPI_COMM_WORLD);
	motSA.nNucMe = nNucAll[rank];
	sum = SumArrInt(nNucAll, rank);
	diff = nAbpMall[2] - sum * motSA.nMotPerSide * 2;
	if (diff <= 0) {
		nAbpAddMe[2] = 0;
	}
	else {
		nAbpAddMe[2] = motSA.nNucMe * motSA.nMotPerSide * 2;
		if (diff - motSA.nNucMe * motSA.nMotPerSide * 2
				< motSA.nMotPerSide * 2 && motSA.nNucMe > 0) {
			nAbpAddMe[2] += diff - motSA.nNucMe * motSA.nMotPerSide * 2;
		}
	}
	free(nNucAll);
  }
  nAbpAddSumMe = SumArrInt(nAbpAddMe, NK_ABP);
  MALLOC(nAbpAddSumAll, int, nCpu);
  MPI_Allgather(&nAbpAddSumMe, 1, MPI_INT, nAbpAddSumAll, 1, MPI_INT,
 		MPI_COMM_WORLD);
  sum = SumArrInt(nAbpAddSumAll, rank);
  free(nAbpAddSumAll);

  for(n = 0; n < NK_ABP; n++) {
    nAbpDet[n] += nAbpMall[n];
    nAbpMall[n] = 0;
  }
  nAcpMme += nAbpAddMe[0] + nAbpAddMe[1];
  nMotMme += nAbpAddMe[2];
  // Shift the copied ABPs
  for(n = nAbpCp - 1; n >= 0; n--) {
	UpdateActinAbpMonomerListSubroutine(nAbpMe + nAbpAddSumMe + n,
			nAbpMe + n, 1);
  }

  cnt = 0;
  // Determine the kind of free ABP to add
  for(n = 0; n < nAbpAddSumMe; n++) {
	CS = 0;
	while(CS == 0) {
		kind = GenRandIntIndex(NK_ABP);
		if (nAbpAddMe[kind] > 0) {
			nAbpAddMe[kind]--;
			CS = 1;
		}
	}
	if (kind == 2) {
		abp.kind[nAbpMe + n] =  (cnt < (int)(rho.fracM * (double)nAbpAddMe[2])) 
				? 3 : 2;
		cnt++;
	}
	// Add the free ABP
	UpdateActinAbpMonomerListSubroutine2(nAbpMe + n, nAbp + sum + n, kind + 1);


	if ((kind != 2 && gTglImpAcpM == 0) || (kind == 2 && gTglImpMotM == 0)) {
		GenRandPosSubdom(&P2(abp.r,nAbpMe + n,0));
	}
	if (kind != 2 && gTglImpAcpM != 0) {
		InsertElement1dArrayWoChk(acpM.l, &acpM.c, nAbp + sum + n);
	}
	if (kind == 2 && gTglImpMotM != 0) {
		InsertElement1dArrayWoChk(motM.l, &motM.c, nAbp + sum + n);
	}
  }
  // All other free ABPs belonging to other submains should be not here
  for (n = 0; n < nAbpMallSum; n++) {
	if (!(n >= sum && n < sum + nAbpAddSumMe)) {
		iAbp[nAbp + n] = -1;
	}
  }
  nAbpMe += nAbpAddSumMe;
  nAbp += nAbpMallSum;

  motM.cM = (int)(rho.fracM * (double)motM.c);
  motM.cOri = motM.c - motM.cM;
  motM.cMOri = motM.cM;
  
  motSA.nNucMeOri = motSA.nNucMe;
  motSA.nNucMeM = (int)(rho.fracM * (double)motSA.nNucMe);
  motSA.nNucMeMOri = motSA.nNucMeM;

}

void AddFreeActin(void) {
  int n, nActAddMe, sum;

  if (nActMall > 0) {
	nActAddMe = nActMall / nCpu;
	sum = rank * nActAddMe;
	for (n = nActCp - 1; n >= 0; n--) {
		UpdateActinAbpMonomerListSubroutine(nActMe + nActAddMe + n,
				nActMe + n, 0);
	}
	for (n = 0; n < nActAddMe; n++) {
		UpdateActinAbpMonomerListSubroutine2(nActMe + n,
				nAct + sum + n, 0);
		actM.l[n + actM.c] = nAct + sum + n;
	}
	for (n = 0; n < nActMall; n++) {
		CONT(n >= sum && n < sum + nActAddMe);
		iAct[nAct + n] = -1;
	}
	actM.c += nActAddMe;
	nActMe += nActAddMe;
	nAct += nActMall;
	if (rank == 0) {
  		Printf0("%d actins are added.\n", nActMall);
	} 
  } 
}

/*----------------------- Adding and deleting elements -----------------------*/

/*--------------------------- Pack information -------------------------------*/

// 0: actin0, 2: stalling at barbed end, 4: stalling by traffic jam; 
// 6: stalling by forces
// 1, 3, 5, 7: same for actin 1
// 8: the kind of ABP
int BinaryPackAbpChainArray(int abpInd) {
  int k, ch, side, actInd, locActInd, locAbpInd, locNextActInd;
  double pMotW, critPMotW;

  critPMotW = 1.0e-10;
  ch = 0;
  locAbpInd = iAbp[abpInd];
  for(k = 0; k < 2; k++) {
	actInd = P2A(abp.ch,locAbpInd,k,nChAb);
	// Check whether any actin is bound on ABP
	if (actInd > -1) { 
		ch += (int)(pow(2, k)); 
		locActInd = iAct[actInd];
		CONT(locActInd < 0);
	}
	else { continue; }
	locNextActInd = iAct[P2A(act.ch,locActInd,0,nChAc)];
	// Check traffic jam
	side = FindAbpActinChain(locActInd, abpInd, 0);
	side = (side - 2) / nChAcY;
	if (side < nChAcX - 1) {
		side = (((side - 2) / nChAcY) + 1) * nChAcY + 2;
		side = FindElementArray(&P2A(act.ch,locActInd,side,nChAc), nChAcY, 
				-1, 0, 1);
	}
	else {
		if (locNextActInd > -1) {
			side = FindElementArray(&P2A(act.ch,locNextActInd,2,nChAc), nChAcY, 
					-1, 0, 1);
			// Check whether ABP is located at the barbed ends
			if (P2A(act.ch,locNextActInd,0,nChAc) < 0) {
				ch += (int)pow(2, k + 2);
			}
		}
	}
	if (side < 0) { 
		ch += (int)pow(2., k + 4);
	}
	CONT(K_ABP(locAbpInd) != 2 || motWalk.gTgl == 0);
	// Check stalling by forces
	pMotW = UpdateMotorWalkingSubroutine(abpInd, k);
	if (pMotW < critPMotW) {
		ch += (int)pow(2., k + 6);
	}
  }
  ch += (int)pow(2., K_ABP(locAbpInd) + 8);
  if (ISMTF(K_ABP(locAbpInd))) {
	if (abp.mId[locAbpInd] > -1) {
		ch += (int)pow(2., abp.mId[locAbpInd] + 11);
	}
  }
  return ch;
}

int BinaryPackActinChainArray(int actInd) {
  int k, ch, cntAbp;

  ch = 0;
  cntAbp = 0;
  for(k = 0; k < nChAc; k++) {
	CONT(!(P2A(act.ch,iAct[actInd],k,nChAc) > -1));
	if (k < 2) { ch += (int)(pow(2, k + 1)); }
	else { cntAbp++; }
  }
  k = 3;
  while(cntAbp != 0) {
	if (cntAbp % 2 == 1) { ch += (int)(pow(2, k)); }
	cntAbp /= 2;
	k++;
  }
  return ch;
}

/*--------------------------- Pack information -------------------------------*/
