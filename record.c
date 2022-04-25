// ##################################################
// #   record.c - finally revised on Apr 2022       #
// #   coded by Taeyoon Kim                         #
// #   Copyright (C) 2005 - 2022, Taeyoon Kim,      #
// #   All rights reserved.                         #
// ##################################################
// This file contains functions which record data.

/*------------------- Recording the progress of simulations ------------------*/

void RecordProgressSubroutine(int *cntMe, FILE *fOut) {
  int *cntAll, sumCnt;

  MALLOC(cntAll,int,nCpu);
  MPI_Gather(cntMe, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	sumCnt = SumArrInt(cntAll, nCpu);
	fprintf(fOut, "%d\t", sumCnt);
	Printf("%5d ", sumCnt);
  }
  free(cntAll);
}

void RecordProgressSubroutine2(int *nPar1, int *nPar2, int nPar3, 
		FILE *fOut, int mode) {
  int n, *cntAll, sumCnt[2];

  MALLOC(cntAll,int,nCpu);
  V2SET_ALL(sumCnt, 0);
  for(n = 0; n < 2; n++) {
	CONT(!((n == 0 && (int)(mode / 2) == 1) || (n == 1 && mode % 2 == 1)));
	MPI_Gather(((n == 0) ? nPar1 : nPar2), 1, MPI_INT, cntAll, 
			1, MPI_INT, 0, MPI_COMM_WORLD);
	CONT(rank != 0);
	sumCnt[n] = SumArrInt(cntAll, nCpu);
	fprintf(fOut, "%d\t", sumCnt[n]);
	Printf("%5d/", sumCnt[n]);
  }
  if (rank == 0) {
	if (mode == 0) { 
		fprintf(fOut, "%d\t", nPar3);
		Printf("%5d ", nPar3);
	}
	else {
		fprintf(fOut, "%d\t%d\t", nPar3 - sumCnt[0] - sumCnt[1], nPar3);
		Printf("%5d/%5d ", nPar3 - sumCnt[0] - sumCnt[1], nPar3);
	}
  }
  free(cntAll);
}

void RecordProgressSubroutine3(int *cntMe, int *cntMe2, FILE *fOut) {
  int *cntAll, sumCnt[2];

  MALLOC(cntAll,int,nCpu);
  MPI_Gather(cntMe, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	sumCnt[0] = SumArrInt(cntAll, nCpu);
  }
  MPI_Gather(cntMe2, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	sumCnt[1] = SumArrInt(cntAll, nCpu);
	fprintf(fOut, "%d\t", sumCnt[0] + sumCnt[1]);
	Printf("%5d ", sumCnt[0] + sumCnt[1]);
  }
  free(cntAll);
}

// Record the progress of the simulation.
void RecordProgress(void) {
  int n, nAcp, mode, nMot2[3], *pArr;
  time_t now;
  FILE *fOut;

  if (rank == 0) {
	fOut = fopen(GenFileName("Progress"), "a");
	fprintf(fOut, "%lld\t%g\t", currTimeStep, (double)currTimeStep * dtReal);
	Printf("%12d ", currTimeStep);
  }
  nAcp = nAbp - nMot;
  // Gather information from all subdomains first.
  if (nAct > 0) {
	mode = 0;
	if (actAss.gTgl != 0 || actDis.gTgl != 0) { mode++; }
	RecordProgressSubroutine2(&actM.c, &actM.c, nAct, fOut, mode);
	if (actNuc.gTgl != 0 || actDis.gTgl != 0) 
	{ RecordProgressSubroutine(&nActFilaMe, fOut); }
	if (actNuc.gTgl != 0) { RecordProgressSubroutine(&actNuc.cntMe, fOut); }
	if (actAss.gTgl != 0) { RecordProgressSubroutine(&actAss.cntMe, fOut); }
	if (actDis.gTgl != 0) 
	{ RecordProgressSubroutine(&actDis.cntMe, fOut); }
  }
  if (nAcp > 0) { 
	mode = 0;
	if (acpInaUnb.gTgl != 0 || acpMoBind.gTgl != 0) { mode += 2; }
	if (acpUnb.gTgl != 0 || acpReb.gTgl != 0) { mode++; }
	RecordProgressSubroutine2(&nAcpMme, &nAcpInaMe, nAcp, fOut, mode); 
	if (acpInaUnb.gTgl != 0) 
	{ RecordProgressSubroutine(&acpInaUnb.cntMe, fOut); }
	if (acpMoBind.gTgl != 0) 
	{ RecordProgressSubroutine(&acpMoBind.cntMe, fOut); }
	if (acpUnb.gTgl != 0) { RecordProgressSubroutine(&acpUnb.cntMe, fOut); }
	if (acpReb.gTgl != 0) { RecordProgressSubroutine(&acpReb.cntMe, fOut); }
  }
  if (nMot > 0) { 
	mode = 0;
	if (motInaUnb.gTgl != 0 || motMoBind.gTgl != 0) { mode += 2; }
	if (motUnb.gTgl != 0 || motReb.gTgl != 0) { mode++; }
	if (motSA.gTgl == 0) { 
		RecordProgressSubroutine2(&nMotMme, &nMotInaMe, nMot, fOut, mode); 
		if (motInaUnb.gTgl != 0) 
		{ RecordProgressSubroutine(&motInaUnb.cntMe, fOut); }
		if (motMoBind.gTgl != 0) 
		{ RecordProgressSubroutine(&motMoBind.cntMe, fOut); }
		if (motUnb.gTgl != 0) { RecordProgressSubroutine(&motUnb.cntMe, fOut); }
		if (motReb.gTgl != 0) { RecordProgressSubroutine(&motReb.cntMe, fOut); }
	}
	else {
		nMot2[0] = 0;
		FOR_ABPME(n) {
			pArr = &P2A(abp.ch,n,0,nChAb);
			CONT(pArr[2] != 2);
			CONT(!(pArr[0] < 0 && pArr[1] < 0 && pArr[3] < 0 && pArr[4] < 0));
			nMot2[0] += 2;
		}
		nMot2[1] = 2 * nMotMme - nMot2[0] + nMotInaMe;
		nMot2[2] = nMot * 2;
		RecordProgressSubroutine2(&nMot2[0], &nMot2[1], nMot2[2], fOut, mode); 
		if (motUnb.gTgl != 0) {
			RecordProgressSubroutine3(&motUnb.cntMe, &motInaUnb.cntMe, fOut); 
		}
		if (motReb.gTgl != 0) { 
			RecordProgressSubroutine3(&motReb.cntMe, &motMoBind.cntMe, fOut); 
		}
	}
	if (motWalk.gTgl != 0) { RecordProgressSubroutine(&motWalk.cntMe, fOut); }
    if (motSA.gTgl != 0) {
		RecordProgressSubroutine(&motSA.cntNucMe, fOut);
		RecordProgressSubroutine(&motSA.cntAssMe, fOut);
		if (motSA.to.gTgl != 0) {
			RecordProgressSubroutine(&motSA.cntTurnMe, fOut);
		}
    }
  }
  if (rho.tgl != 0) {
	if (nAct > 0) { RecordProgressSubroutine(&rho.cntActFor, fOut); }
	if (nAct > 0) { RecordProgressSubroutine(&rho.cntActAss, fOut); }
	if (nAct > 0) { RecordProgressSubroutine(&rho.cntMotNuc, fOut); }
	if (nMot > 0) { RecordProgressSubroutine(&rho.cntMotAss, fOut); }
  }

  if (rank == 0) {
	now = time(NULL);
	fprintf(fOut, "%s", ctime(&now));
	fclose(fOut);
	Printf("\n");
  }

}

/*------------------- Recording the progress of simulations ------------------*/

/*------------------------- Recording parameter values -----------------------*/

void RecordInitParameterSubroutine(const char *str, int tgl1,
		int tgl3, FILE *fOut) {
  fprintf(fOut, "%s = %s / %s\n", str, 
		((tgl1 != 0 && tgl3 != 0) ? "on" : "off"),
		(tgl1 != 0 ? "on": "off"));
}

void RecordInitParameterSubroutine2(const char *str, 
		FuncCont *cont, FILE *fOut) {
  fprintf(fOut, "%s = %s\n", str, (cont->tgl != 0) ? "on" : "off");
  if (cont->tgl != 0) {
	fprintf(fOut, "- Period = %g [s] (%d timesteps)\n", cont->prdR, cont->prd);
  }
}

// Record the parameter values which were initially loaded from "condition" 
// and "Config".
void RecordInitParameter(void) {
  int chk[NDIM][2], chk2, k, nAcp;
  double dimDomC[NDIM] , vol, ratio, ratio2, cAbp, rAbp;
  char direc[4] = "xyz";
  FILE *fOut;  

  nAcp = nAbp - nMot;

  V3DIV(dimDomC, dimDom, nCell);
  vol = V3PROD(dimDom) * CUBE(L_SCALE_IN_M) * 1e3;
  cAbp = (double)nAbp * 2. / N_AVO / vol * 1e6;
  rAbp = (double)nAbp / ((double)nAct * nActPerSeg);
  fOut = fopen(GenFileName("Parameter"), "w");

  fprintf(fOut, "Duration of initial network formation = %g [s] "
		"(%lld timesteps)\n", netForm.durR, netForm.dur);
  fprintf(fOut, "Time step = %g [s] (%g)\n", dtReal, dt);
  fprintf(fOut, "Threshold of unstable force = %g [pN] (%g)\n",
		MAG_UNSTABLE_FORCE * 1e12, magUnstF);

  fprintf(fOut, "\n================== Number and concentration of actins and "
		"ABPs =================\n");
  fprintf(fOut, "Number of actin segments = %d (C_A = %g uM)\n", nAct, 
		(double)nAct * 2. * nActPerSeg / N_AVO / vol * 1e6);
  fprintf(fOut, "Length of actin segments = %g [nm]\n", L_SCALE_IN_NM);
  fprintf(fOut, "Number of ABPs = %d (C_ABP = %g uM, R_ABP = %g)\n", nAbp, 
		cAbp, rAbp);
  if (nAbp > 0) {
	ratio = (nAbp > 0) ? (double)nMot / (double)nAbp : 0.;
	ratio2 = (nAcp > 0) ? (double)nAbpDet[0] / (double)nAcp : 0.;
	fprintf(fOut, "(N_ACP^C = %d, C_ACP^C = %g uM, R_ACP^C = %g)\n", nAbpDet[0],
			cAbp * (1. - ratio) * ratio2, rAbp * (1. - ratio) * ratio2);
	fprintf(fOut, "(N_ACP^B = %d, C_ACP^B = %g uM, R_ACP^B = %g)\n", nAbpDet[1],
			cAbp * (1. - ratio) * (1. - ratio2), rAbp * (1. - ratio) 
			* (1. - ratio2));
	fprintf(fOut, "(N_M = %d, C_M = %g uM, R_M = %g)\n", nMot, 
			cAbp * ratio, rAbp * ratio);
  }

  fprintf(fOut, "\n===================== Conditions for domain and boundaries "
		"=====================\n");
  fprintf(fOut, "Domain width: x = %g [um] (%g), y = %g [um] (%g),"
		" z = %g [um] (%g)\n", L_S2UM(dimDom[0]), dimDom[0], 
		L_S2UM(dimDom[1]), dimDom[1], L_S2UM(dimDom[2]), dimDom[2]);
  fprintf(fOut, "Periodic boundary condition: x = %s, y = %s, z = %s\n", 
		(pbc[0] == 1 ? "on" : "off"), (pbc[1] == 1 ? "on" : "off"),
		(pbc[2] == 1 ? "on" : "off"));
  fprintf(fOut, "Geometry of a domain = ");
  if (dir2D > -1) {
	fprintf(fOut, "2-D with a normal direction in %c\n", direc[dir2D]);
  }
  else {
	fprintf(fOut, "3-D\n");
  }

  fprintf(fOut, "\n======================== Information about measurement "
			"=========================\n");
  fprintf(fOut, "Net duration of simulation = %g [s] "
		"(%lld timesteps)\n", rheo.durR, rheo.dur);

  if (nAct > 0) {
	fprintf(fOut, "\n================== Conditions for dynamic behaviors of "
			"actins ==================\n");
	RecordInitParameterSubroutine("Nucleation of actins", 
			actNuc.gTgl, gTglActDynNF, fOut);
	RecordInitParameterSubroutine("Assembly of actins", 
			actAss.gTgl, gTglActDynNF, fOut);
	RecordInitParameterSubroutine("Disassembly of actins", 
			actDis.gTgl, gTglActDynNF, fOut);
  }
  if (nAcp > 0) {
	fprintf(fOut, "\n=================== Conditions for dynamic behaviors "
			"of ACPs ===================\n");
	RecordInitParameterSubroutine("Unbinding of inactive ACPs", 
			acpInaUnb.gTgl, gTglAcpDynNF, fOut);
	RecordInitParameterSubroutine("Binding of monomeric ACPs", 
			acpMoBind.gTgl, gTglAcpDynNF, fOut);
	RecordInitParameterSubroutine("Unbinding of active ACPs", 
			acpUnb.gTgl, gTglAcpDynNF, fOut);
	RecordInitParameterSubroutine("Binding of inactive ACPs", 
			acpReb.gTgl, gTglAcpDynNF, fOut);
	fprintf(fOut, "Implicit consideration of monomeric ACPs = %s\n", 
			(gTglImpAcpM != 0 ? "on" : "off"));
  }
  if (nMot > 0) {
	fprintf(fOut, "\n================== Conditions for dynamic behaviors "
			"of motors ==================\n");
	RecordInitParameterSubroutine("Walking of active motors", 
			motWalk.gTgl, gTglMotWalkNF, fOut);
	fprintf(fOut, "Self-assembly of motors = %s\n",
			(motSA.gTgl != 0 ? "on" : "off"));
	if (motSA.gTgl != 0) {
		RecordInitParameterSubroutine("Unbinding of motors", 
				motUnb.gTgl, gTglMotUnbRebNF, fOut);
		RecordInitParameterSubroutine("Binding of motors", 
				motReb.gTgl, gTglMotUnbRebNF, fOut);
		fprintf(fOut, "Check the alignment between actin and motor filaments "
				"for binding of motors = %s\n", (motReb.gTglCrsAng != 0 
				? "on" : "off"));
		fprintf(fOut, "Allow the turnover of motor filaments = %s\n", 
				(motSA.to.gTgl != 0) ? "on" : "off");
	}
	else {
		RecordInitParameterSubroutine("Unbinding of inactive motors", 
				motInaUnb.gTgl, gTglMotUnbRebNF, fOut);
		RecordInitParameterSubroutine("Binding of monomeric motors", 
				motMoBind.gTgl, gTglMotUnbRebNF, fOut);
		RecordInitParameterSubroutine("Unbinding of active motors", 
				motUnb.gTgl, gTglMotUnbRebNF, fOut);
		RecordInitParameterSubroutine("Binding of inactive motors", 
				motReb.gTgl, gTglMotUnbRebNF, fOut);
		fprintf(fOut, "Implicit consideration of monomeric motors = %s\n", 
				(gTglImpMotM != 0 ? "on" : "off"));
	}
  }
  fprintf(fOut, "\n======== Adjustment for rates of the dynamic behaviors of "
		"actin and ABP ========\n");
  fprintf(fOut, "During network formation(i*k) = %g\n", netForm.facK);

  if (nAct > 0) {
	if (actNuc.gTgl != 0 || actAss.gTgl != 0 || actDis.gTgl != 0) {
		fprintf(fOut, "\n======================== Parameters of actin dynamics "
				"=========================\n");
	}
	if (actNuc.gTgl != 0) {
		fprintf(fOut, "Nucleation rate = %g [1/uM s]\n", actNuc.k);
	}
	if (actAss.gTgl != 0) {
		fprintf(fOut, "Assembly rate = %g [1/uM s] at barbed, %g [1/uM s] at "
				"pointed\n", actAss.k[0], actAss.k[1]);
	}
	if (actDis.gTgl != 0) {
		fprintf(fOut, "Disassembly rate = %g [1/s] at barbed, %g [1/s] at "
				"pointed\n", actDis.k[0], actDis.k[1]);
		fprintf(fOut, "Factor for varying disassembly rate of actin with "
				"ABPs = %g\n", actDis.facKWA);
	}
  }

  if (nAcp > 0) {
	if (acpMoBind.gTgl != 0 || acpReb.gTgl != 0 || acpInaUnb.gTgl != 0 
			|| acpUnb.gTgl != 0) {
		fprintf(fOut, "\n================= Parameters for unbinding and "
				"binding of ACPs =================\n");
	}
	if (acpMoBind.gTgl != 0 || acpReb.gTgl != 0) {
		fprintf(fOut, "Binding and binding rate(k) = %g [1/uM s]\n", 
				K_ACP_BIND);
		fprintf(fOut, "Factor that adjusts k = %g\n", acpReb.facK);
	}
  }

  if (nMot > 0) {
	if (motSA.gTgl != 0) { 
		fprintf(fOut, "\n====== Parameters for motor self-assembly, "
				"unbinding, binding, and walking =====\n");
	}
	else {
		fprintf(fOut, "\n============== Parameters for motor unbinding, "
				"binding, and walking ============\n");
	}
	if (motMoBind.gTgl != 0 || motReb.gTgl != 0) {  
		fprintf(fOut, "Binding rate(k) = %g [1/uM s]\n", 
				40. * motMC.nHead);
		fprintf(fOut, "Factor that adjusts k = %g\n", motReb.facK);
	}
	if (motInaUnb.gTgl != 0 || motUnb.gTgl != 0) {
		fprintf(fOut, "Unbinding rate without applied force(k0) = %g [1/s]\n", 
				(motUnb.gTgl != 0 || motInaUnb.gTgl != 0) ? 
				log(1 - motUnb.p[motUnb.maxF[0] - 1]) / REVSIGN(dtReal) : 0.);
	}
	if (motWalk.gTgl != 0) {
		fprintf(fOut, "Walking rate without applied force(k0) = %g [1/s]\n", 
				(motWalk.gTgl != 0) ? log(1 - motWalk.p[motWalk.maxF[0] - 1]) 
				/ REVSIGN(dtReal) : 0.);
		fprintf(fOut, "Factor that adjusts k0 of unbinding and walking = %g\n", 
				motUnb.facK0);
		fprintf(fOut, "Each walking distance = %g [nm]\n", 
				L_SCALE_IN_NM / (double)nChAcX);
	}
	if (motMoBind.gTgl != 0 || motReb.gTgl != 0 || motInaUnb.gTgl != 0 
			|| motUnb.gTgl != 0 || motWalk.gTgl != 0) {
		fprintf(fOut, "Transition rates between mechanochemical states in each "
				"head(k01, k10, k12, k21, k20) = %g, %g, %g, %g, %g [1/s]\n", 
				motMC.k01, motMC.k10, motMC.k12, motMC.k21, motMC.k20);
		fprintf(fOut, "Number of heads which each motor arm represents = %d\n", 
				motMC.nHead);
	}
	if (motSA.gTgl != 0) {
		fprintf(fOut, "Average number of motors per each self-assembled "
				"structure = %d\n", motSA.nMotPerTF);
		fprintf(fOut, "k for motor assembly = %g [1/s]\n", motSA.kAss);
		fprintf(fOut, "Turnover of motor filaments = %s\n", 
				(motSA.to.gTgl != 0) ? "on" : "off");
		if (motSA.to.gTgl != 0) {
			fprintf(fOut, "k for motor turnover = %g [1/s]\n", motSA.to.k);
		}
	}
  }

  if (nAbp > 0) {
	fprintf(fOut, "\n============================== Geometry of ABP"
			" =================================\n");
	fprintf(fOut, "Number of binding sites for ABPs on each actin segment in "
			"longitudinal and transverse directions = %d, %d\n", nChAcX, 
			nChAcY);
  }
  if (nAbpDet[0] > 0) {
	fprintf(fOut, "Length of ACP^C arm = %g [nm] (%g)\n", 
			L_ACPC_ARM * 1.0e9, L_M2S(L_ACPC_ARM));
  }
  if (nAbpDet[1] > 0) {
	fprintf(fOut, "Length of ACP^B arm = %g [nm] (%g)\n", 
			L_ACPB_ARM * 1.0e9, L_M2S(L_ACPB_ARM));
  }
  if (nMot > 0) {
	fprintf(fOut, "Length of motor arm = %g [nm] (%g)\n", 
			L_MOT_ARM * 1.0e9, L_M2S(L_MOT_ARM));
	if (motSA.gTgl != 0) {
		fprintf(fOut, "Length of bare zone in self-assembled motors = %g [nm] "
				"(%g)\n", L_MOTBACK_CEN_DIST * 1.0e9, 
				L_M2S(L_MOTBACK_CEN_DIST));
		fprintf(fOut, "Spacing between motors in self-assembled motors = %g "
				"[nm] (%g)\n", L_MOTBACK_DIST * 1.0e9, L_M2S(L_MOTBACK_DIST));
	}
  }
  fprintf(fOut, "\n=============================== Drag coefficients "
		"==============================\n");
  fprintf(fOut, "Viscosity of medium = %g [Pa s] (%g)\n", VISCOSITY, 
		VISCOSITY * L_SCALE_IN_M / actF.dragR);
  if (nAct > 0) {
	fprintf(fOut, "Drag coefficient of actin segment = %g [kg/s] (%g)\n", 
			actF.dragR, 1.0);
  }
  if (nAbpDet[0] > 0) {
	fprintf(fOut, "Drag coefficient of ACP^C = %g [kg/s] (%g)\n", 
			abpF.drag[0].n * actF.dragR, abpF.drag[0].n);
  }
  if (nAbpDet[1] > 0) {
	fprintf(fOut, "Drag coefficient of ACP^B = %g [kg/s] (%g)\n", 
			abpF.drag[1].n * actF.dragR, abpF.drag[1].n);
  }
  if (nMot > 0) {
	fprintf(fOut, "Drag coefficient of motor = %g [kg/s] (%g)\n", 
			abpF.drag[2].n * actF.dragR, abpF.drag[2].n);
  }

  fprintf(fOut, "\n============================== Thermal fluctuation "
		"=============================\n");
  if (nAct > 0) {
	fprintf(fOut, "Thermal fluctuation of actin filaments = %s\n",
			(gTglActTherm != 0 ? "on" : "off"));
  }
  if (nAcp > 0) {
	fprintf(fOut, "Thermal fluctuation of ACPs = %s\n", 
			(gTglAcpTherm != 0 ? "on" : "off"));
  }
  if (nMot > 0) {
	fprintf(fOut, "Thermal fluctuation of motors = %s\n", 
			(gTglMotTherm != 0 ? "on" : "off"));
  }

  if (nAct > 0 || nAbp > 0) {
	fprintf(fOut, "\n========================= Stiffness of force models "
			"============================\n");
  }
  if (nAbp > 0) {
  fprintf(fOut, "Strength of repulsive forces between ABPs and between ABP "
		"and actin = %g [N/m] (%g)\n",	KS_S2NPM(abpF.rep.stf), 
		abpF.rep.stf);
  }
  if (nAct > 0) {
	fprintf(fOut, "------------------------------------- actin "
			"------------------------------------\n");
	fprintf(fOut, "Strength of repulsive forces between actin filaments = %g "
			"[N/m] (%g)\n",	KS_S2NPM(actF.rep.stf), actF.rep.stf);
	fprintf(fOut, "Bending stiffness of actin filaments = %g [Nm] (%g)\n",
			KB_S2NM(actF.bend.stf), actF.bend.stf);
	fprintf(fOut, "Persistence length of actin filaments = %g [um]\n", 
			KB_S2NM(actF.bend.stf) * L_SCALE_IN_M / KT_IN_J * 1e6);
	fprintf(fOut, "Extensional stiffness of actin filaments = %g [N/m] (%g)\n",
			KS_S2NPM(actF.spr.stf), actF.spr.stf);
  }
  if (nAbpDet[0] > 0) {
	fprintf(fOut, "------------------------------------- ACP^C "
			"------------------------------------\n");
	fprintf(fOut, "Bending stiffness of ACP^C (theta 1) = %g [Nm] (%g)\n",
			KB_S2NM(abpF.bend[0].stf), abpF.bend[0].stf);
	fprintf(fOut, "Bending stiffness of ACP^C (theta 2) = %g [Nm] (%g)\n",
			KB_S2NM(abpF.a90[0].stf), abpF.a90[0].stf);
	fprintf(fOut, "Extensional stiffness of ACP^C = %g [N/m] (%g)\n",
			KS_S2NPM(abpF.spr[0].stf), abpF.spr[0].stf);
  }
  if (nAbpDet[1] > 0) {
	fprintf(fOut, "------------------------------------- ACP^B "
			"------------------------------------\n");
	fprintf(fOut, "Bending stiffness of ACP^B (theta 1) = %g [Nm] (%g)\n",
			KB_S2NM(abpF.bend[1].stf), abpF.bend[1].stf);
	fprintf(fOut, "Bending stiffness of ACP^B (theta 2) = %g [Nm] (%g)\n",
			KB_S2NM(abpF.a90[1].stf), abpF.a90[1].stf);
	fprintf(fOut, "Extensional stiffness of ACP^B = %g [N/m] (%g)\n",
			KS_S2NPM(abpF.spr[1].stf), abpF.spr[1].stf);
  }
  if (nMot > 0) {
	fprintf(fOut, "------------------------------------- motor "
			"------------------------------------\n");
	fprintf(fOut, "Bending stiffness of motor (theta 1) = %g [Nm] (%g)\n",
			KB_S2NM(abpF.bend[2].stf), abpF.bend[2].stf);
	fprintf(fOut, "Bending stiffness of motor (theta 2) = %g [Nm] (%g)\n",
			KB_S2NM(abpF.a90[2].stf), abpF.a90[2].stf);
	fprintf(fOut, "Extensional stiffness 1 of motor = %g [N/m] (%g)\n",
			KS_S2NPM(abpF.spr[2].stf), abpF.spr[2].stf);
	fprintf(fOut, "Extensional stiffness 2 of motor = %g [N/m] (%g)\n",
			KS_S2NPM(abpF.spr[2].stf2), abpF.spr[2].stf2);
	if (motSA.gTgl != 0) {
	  	fprintf(fOut, "Bending stiffness of motor backbone = %g [Nm] (%g)\n",
				KB_S2NM(motSA.bend.stf), motSA.bend.stf);
	  	fprintf(fOut, "Extensional stiffness of motor backbone = %g [N/m] "
				"(%g)\n", KS_S2NPM(motSA.spr.stf), motSA.spr.stf);
	}
  }
  fprintf(fOut, "\n============================ Parallel processing "
			"===============================\n");
  fprintf(fOut, "Number of cores used for the job = %d\n", nCpu);
  fprintf(fOut, "Number of subdomains: %d in x, %d in y, %d in z\n", nCell[0], 
		nCell[1], nCell[2]);
  fprintf(fOut, "Initial Subdomain size: x = %g [um] (%g), y = %g "
		"[um] (%g), z = %g [um] (%g)\n", L_S2UM(dimDomC[0]), 
		dimDomC[0], L_S2UM(dimDomC[1]), dimDomC[1], L_S2UM(dimDomC[2]), 
		dimDomC[2]);
  fprintf(fOut, "Thickness of overlapping regions = %g [um] (%g)\n",
		L_S2UM(neiEdge), neiEdge);
  FOR_NDIM(k) {
	chk[k][0] = (nCell[k] > 1) ? 1 : 0; 
	chk[k][1] = (nCell[(k+1)%NDIM] > 1 && nCell[(k+2)%NDIM] > 1) ? 1 : 0;
  }
  chk2 = (nCell[0] > 1 && nCell[1] > 1 && nCell[2] > 1)  ? 1 : 0;
  fprintf(fOut, "The ratio of an overlapping region per to a subdomain = %g "
		"[percent]\n", (2. * neiEdge * (dimDomC[0] * dimDomC[1] * chk[2][0]
		+ dimDomC[1] * dimDomC[2] * chk[0][0] + dimDomC[2] * dimDomC[0] * 
		chk[1][0]) - 4. * SQR(neiEdge) * (dimDomC[0] * chk[0][1] + dimDomC[1] 
		* chk[1][1] + dimDomC[2] * chk[2][1]) + 8. * CUBE(neiEdge) * chk2) 
		/ V3PROD(dimDomC) * 100.);	

  fprintf(fOut, "\n=============================== Data recording "
		"=================================\n");
  fprintf(fOut, "Directory for saving data = %s\n", dataFold);
  fprintf(fOut, "Deletion of pre-existing data files = %s\n", 
		(DELETE_FILE != 0 ? "on" : "off"));
  RecordInitParameterSubroutine2("Record Output and Progress", &recProg, fOut);

  fprintf(fOut, "----------------------------- Structural information "
		"---------------------------\n");
  RecordInitParameterSubroutine2("Record Config", &recConf, fOut);
  fprintf(fOut, "Recording structural information for visualization via VMD "
		"= %s\n", (recConfVmd.tgl != 0) ? "on" : "off");
  if (recConfVmd.tgl != 0) {
	fprintf(fOut, "- Period = %g [s] (%d timesteps)\n", recConfVmd.prdR, 
			recConfVmd.prd);
	if (recConfVmd.mode == 0) { fprintf(fOut, "- Multiple files\n"); }
	else { fprintf(fOut, "- Single file\n"); }
	fprintf(fOut, "- Show boundaries in the network drawn via VMD = %s\n",
			(recConfVmd.gTglBnd != 0 ? "on" : "off"));
	fprintf(fOut, "- Record information for coloring a network in VMD = %s\n",
			(recConfVmd.gTglInfo != 0 ? "on" : "off"));
  }
  RecordInitParameterSubroutine2("Record the length of actin filaments", 
		&recFilaL, fOut);

  fprintf(fOut, "--------------------------- Force, stress, and energy "
		"--------------------------\n");
  RecordInitParameterSubroutine2("Record longitudinal forces acting on ABPs", 
		&recLongF, fOut);
  RecordInitParameterSubroutine2("Record the mechanical energy of a network", 
		&recE, fOut);
  
fprintf(fOut, "--------------------- Dynamic behaviors of actin and ABPs "
		"----------------------\n");
  RecordInitParameterSubroutine2("Record ABP unbinding", &recAbpUnb, fOut);
  RecordInitParameterSubroutine2("Record ABP binding", &recAbpBind, fOut);
  RecordInitParameterSubroutine2("Record ABP turnover", &recAbpTurn, fOut);

  fprintf(fOut, "---------------------------------- Miscellany "
		"----------------------------------\n");
  RecordInitParameterSubroutine2("Record the information in unit of "
		"filaments", &recInfo, fOut);
  fclose(fOut);
}

/*------------------------- Recording parameter values -----------------------*/

/*-------------------------- Recording network structure ---------------------*/

// Record network configuration. Most of the data are related to the positions
// and chain information of particles.
void RecordConfig(char *fn, int period) {
  int n, k, *nActAll, *nAbpAll;
  double *rActOff, *rAbpOff;
  FILE *fOut;

  MALLOC(nActAll,int,nCpu);
  MALLOC(nAbpAll,int,nCpu);
  MALLOC(rActOff,double,nActC*NDIM);
  MALLOC(rAbpOff,double,nAbpC*NDIM);

  MPI_Gather(&nActMe, 1, MPI_INT, nActAll, 1, MPI_INT, 0,
        MPI_COMM_WORLD);
  MPI_Gather(&nAbpMe, 1, MPI_INT, nAbpAll, 1, MPI_INT, 0,
        MPI_COMM_WORLD);
  if (rank == 0) {
	if (SumArrInt(nActAll, nCpu) != nAct) {
		Printf("Error: the sum of collected nActMe is different from "
				"nAct!\n");
		exit(-1);
	}
	if (SumArrInt(nAbpAll, nCpu) != nAbp) {
		Printf("Error: the sum of collected nAbpMe is different from "
				"nAbp!\n");
		exit(-1);
	}
	
	fOut = fopen(fn, ((period > 0) ? "a" : "w")); 
	if (period > 0) {
		fprintf(fOut, "------------------------------------ %6d ------------"
				"------------------------\n", (int)(currTimeStep / period));
	}
	fprintf(fOut, "nAct : %d\n", nAct);
	fprintf(fOut, "nAbp : %d\n", nAbp);
	fprintf(fOut, "dimDom : %g, %g, %g\n", dimDom[0], 
			dimDom[1], dimDom[2]);
	fprintf(fOut, "nActPerSeg, nChAcX, nChAcY : %d, %d, %d\n\n", nActPerSeg, 
			nChAcX, nChAcY);
  }
  // Offset the positions of particles. This can be recovered later using the
  // information of rGrid below.
  FOR_ACTME(n) { 
	FOR_NDIM(k) {
		P2(rActOff,n,k) = P2(act.r,n,k) - rGrid[k][0];
	}
  }
  FOR_ABPME(n) { 
	FOR_NDIM(k) {
		P2(rAbpOff,n,k) = P2(abp.r,n,k) - rGrid[k][0];
	}
  }
  // Positions and chain information
  if (rank == 0) { fprintf(fOut, "## Position for actin ##\n"); }
  RecordGather2dArrayDouble(nActMe, nActAll, 3, act.id, rActOff, fOut, 0);
  if (rank == 0) { fprintf(fOut, "## Position for ABP ##\n"); }
  RecordGather2dArrayDouble(nAbpMe, nAbpAll, 3, abp.id, rAbpOff, fOut, 0);
  if (rank == 0) { fprintf(fOut, "## Chain for actin ##\n"); }
  RecordGather2dArrayInt(nActMe, nActAll, nChAc, act.id, act.ch, fOut, 0);
  if (rank == 0) { fprintf(fOut, "## Chain for ABP ##\n"); }
  RecordGather2dArrayInt(nAbpMe, nAbpAll, nChAb, abp.id, abp.ch, fOut, 0);
  if (rank == 0) {
	fprintf(fOut, "currTimeStep = %lld\n", currTimeStep);
	fprintf(fOut, "rGrid = %g, %g, %g, %g, %g, %g\n", rGrid[0][0], rGrid[1][0],
			rGrid[2][0], rGrid[0][nGrid[0] - 1], rGrid[1][nGrid[1] - 1], 
			rGrid[2][nGrid[2] - 1]);
	fprintf(fOut, "pbc = %d, %d, %d\n", pbc[0], pbc[1], pbc[2]);
  }
  free(rActOff);
  free(rAbpOff);
  free(nActAll);
  free(nAbpAll);
}

int RecordConfigVmdSubroutine(double *r1, double *r2, int *ind, int *mul) {
  int k, adjInd, fac[NDIM];
  double dr[NDIM];

  V3COPY(fac, ind);
  V3SUB(dr, r1, r2);
  FOR_NDIM(k) {
	CONT(!(fabs(dr[k]) > dimDomH[k]));
	fac[k] += (dr[k] > dimDomH[k] ? 1 : -1);
	if (fac[k] >= mul[k]) { fac[k] = 0; }
	else if (fac[k] < 0) { fac[k] = mul[k] - 1; } 
  }
  V3IND_BACK_INT(adjInd, fac, mul);

  return adjInd;
}

// Record configuration for visualization via VMD. Positions and chain 
// information are recorded in the PDB and PSF formats. Actin, ACP, and motor
// have different SEGNAMEs for distinction. 
// mode = 0: bonds between actins and ABPs are drawn to the end points of 
//           actin segments. 
//        1: actin segments are divided to sub-segments, so the bonds are
//           drawn to actual points. This will cost more memory in VMD.
void RecordConfigVmd(int mode) {
  int m, n, k, ind, ind2[NDIM], ind3, ind4, ind5, cnt, curr, amp = 1;
  int nRec, nD, ch, *pArr, *pArr2, loc[2], nCh;
  int nActVmd, nAbpVmd, mul[NDIM];
  int facInfoConfVmd, facBndConfVmd;
  int *nActAll, *nAbpAll, *chActVmd, *chAbpVmd, *iFila;
  int *chActAlt;
  int *actRho, *abpRho;
  int nBeadVmd;
  int bndCh[] = {0, 4, 1, 5, 2, 6, 3, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 
				2, 1, 3, 4, 6, 5, 7};
  double temp, dr[NDIM], domCen[NDIM], rPos[NDIM], sft[NDIM];
  double *rActVmd, *rAbpVmd, *occu, *mass;
  double *recActG, *recAbpG, *sendRecAct, *sendRecAbp, *pntRec;
  char **segType, fn[80], fnIn[80];
  FILE *fOut;
  ListInt chain;

  if (recConfVmd.mode == 0) {
	sprintf(fnIn, "ConfVmd_%d", (int)(currTimeStep / recConfVmd.prd));
  }
  else {
	sprintf(fnIn, "ConfVmd");
  }

  facInfoConfVmd = (recConfVmd.gTglInfo > 0) ? 1 : 0;
  facBndConfVmd = (recConfVmd.gTglBnd != 0) ? 8 : 0;

  MALLOC(chActAlt,int,nChAc);
  // If a direction has periodic boundary condition, the configuration is
  // automatically duplicated in the direction for a clearer view.
  FOR_NDIM(k) {
	mul[k] = (pbc[k] == 1) ? 2 : 1;
  }
  V3SET_ALL(mul, 1);

  nActVmd = nAct * V3PROD(mul) * ((mode != 0) ? nChAcX : 1);
  nAbpVmd = nAbp * V3PROD(mul);

  MALLOC(nActAll,int,nCpu);
  MALLOC(nAbpAll,int,nCpu);

  if (rank == 0) {
	MALLOC(rAct,double,nAct*NDIM);
	MALLOC(rAbp,double,nAbp*NDIM);
	MALLOC(chAct,int,nAct*nChAc);
	MALLOC(chAbp,int,nAbp*nChAb);
	MALLOC(iFila,int,nAct);
if (rho.tgl != 0) {
	MALLOC(actRho,int,nAct);
	MALLOC(abpRho,int,nAbp);
}
  } 
  // Gather the positions and chain information of particles.
  GatherActChainPosition(nActAll, chAct, rAct);
  GatherAbpChainPosition(nAbpAll, chAbp, rAbp);
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, act.iF, iFila);  

if (rho.tgl != 0) {
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, act.rho, actRho);  
  Gather2dArrayInt(nAbpMe, nAbpAll, 1, abp.id, abp.rho, abpRho);  
}

  // Record a separate data file (.dat) containing information for coloring 
  // methods via VMD.
  nD = 6;
  if (recConfVmd.gTglInfo > 0) {
	MALLOC(sendRecAct,double,nActC * nD);
	MALLOC(sendRecAbp,double,nAbpC * nD);
	if (rank == 0) {
		MALLOC(recActG,double,nAct * nD);
		MALLOC(recAbpG,double,nAbp * nD);
	}
	FOR_ACTME(n) { 
		ch = BinaryPackActinChainArray(act.id[n]);
		V6SET(&P2A(sendRecAct,n,0,nD), recAct.len[n], 
				P2A(recAct.allF,n,NDIM,NDIM + 1), recAct.sprF[n], 
				recAct.bendF[n], ch, (double)recAct.cnt[n]);
	}
	FOR_ABPME(n) { 
		ch = BinaryPackAbpChainArray(abp.id[n]);
		V6SET(&P2A(sendRecAbp,n,0,nD), P2A(recAbp.len,n,0,recAbp.nL), 
				P2A(recAbp.allF,n,NDIM,NDIM + 1), recAbp.sprF[n], 
				recAbp.bendF[n], ch, (double)recAbp.cnt[n]);
	}
	Gather2dArrayDouble(nActMe, nActAll, nD, act.id, sendRecAct, recActG);	
	Gather2dArrayDouble(nAbpMe, nAbpAll, nD, abp.id, sendRecAbp, recAbpG);	

	if (rank == 0) {
		for(m = 0; m < 2; m++) {
			nRec = (m == 0) ? nAct : nAbp;
			pntRec = (m == 0) ? recActG : recAbpG;
			for(n = 0; n < nRec; n++) {
				for(k = 0; k < 4; k++) {
					if (P2A(pntRec,n,nD - 1,nD) == 0) { 
						V4SET_ALL(&P2A(pntRec,n,0,nD), 0.);
					}
					else {
						P2A(pntRec,n,k,nD) /= P2A(pntRec,n,nD - 1,nD);
		 				if (k > 0) {
							P2A(pntRec,n,k,nD) = F_S2PN(P2A(pntRec,n,k,nD));
						}
					}
				}
			}
		}
	}
	free(sendRecAct);
	free(sendRecAbp);
  }

  // Prepare the (duplicated) positions and chain information of particles.
  if (rank == 0) {
	MALLOC(rActVmd, double, nActVmd * NDIM);
	MALLOC(rAbpVmd, double, nAbpVmd * NDIM);
	MALLOC(chActVmd, int, nActVmd * ((mode != 0) ? nChAcY + 2 : nChAc));
	MALLOC(chAbpVmd, int, nAbpVmd * nChAb);
	MALLOC(occu, double, nAct + nAbp);
	MALLOC(mass, double, nAct + nAbp);
	MALLOC2(segType, char, nActVmd + nAbp); 
	MALLOC(chain.l, int, (nActVmd + nAbpVmd * 4 
			+ ((recConfVmd.gTglBnd != 0) ? 12 : 0)) * 2);
	for (k = 0; k < nActVmd + nAbp; k++) {
		MALLOC(segType[k], char, 4);
	}
	memset(chActVmd, -1, sizeof(int) * nActVmd 
			* ((mode != 0) ? nChAcY + 2 : nChAc));
	memset(chAbpVmd, -1, sizeof(int) * nAbpVmd * nChAb);
	// Copy position and chain information to multiplied array
	for(ind = 0; ind < V3PROD(mul); ind++) {
		V3IND_ASSIGN_INT(ind, mul, 1, ind2);
		FOR_ACT(n) {
			pArr = &P2A(chAct,n,0,nChAc);
			memset(chActAlt, -1, sizeof(int) * nChAc);
			for (m = 0; m < nChAc; m++) {
				ind5 = pArr[m];
				CONT(ind5 < 0);
				ind4 = RecordConfigVmdSubroutine(&P2(rAct,n,0), (m < 2 ? 
						&P2(rAct,ind5,0) : &P2(rAbp,ind5,0)), ind2, mul);
				chActAlt[m] = ind5 + ind4 * ((m < 2) ? nAct : nAbp);
			}
			if (mode == 0) {
				ind3 = n + ind * nAct;
				VSV3ADD(&P2(rActVmd,ind3,0), &P2(rAct,n,0), dimDom, ind2);	
				for (m = 0; m < nChAc; m++) {
					P2A(chActVmd,ind3,m,nChAc) = chActAlt[m];
				}
			}
			else {
				ind3 = (n + ind * nAct) * nChAcX;
				VSV3ADD(&P2(rActVmd,ind3,0), &P2(rAct,n,0), dimDom, ind2);	
				pArr2 = &P2A(chActVmd,ind3,0,nChAcY + 2);
				if (pArr[1] > -1) { 
					pArr2[1] = (chActAlt[1] + ind4 * nAct + 1) * nChAcX - 1;
				}
				if (pArr[0] > -1) {
					for(m = 1; m < nChAcX; m++) {
						// Position
						CalcPosOnActSeg(&P2(rAct,n,0), &P2(rAct,pArr[0],0),
								rPos, (double)m / (double)nChAcX);
						VSV3ADD(&P2(rActVmd,ind3 + m,0), rPos, dimDom, ind2);
						// Chain
						V2SET(&P2A(pArr2,m,0,nChAcY + 2), ind3 + m + 1, 
								ind3 + m - 1);
					}
					// Chain at ends
					pArr2[0] = ind3 + 1;
					P2A(pArr2,nChAcX - 1,0,nChAcY + 2) 
							= (chActAlt[0] + ind4 * nAct) * nChAcX;

					for(m = 2; m < nChAc; m++) {
						CONT(pArr[m] < 0);
						loc[0] = (m - 2) / nChAcY;
						loc[1] = (m - 2) % nChAcY;
						P2A(pArr2,loc[0],loc[1] + 2,nChAcY + 2) 
								= chActAlt[m] + ind4 * nAbp;
					}
				}
				else {
					for(m = 1; m < nChAcX; m++) {
						V3SET_ALL(&P2(rActVmd,ind3 + m,0), 0.);
					}
				}

			}
		}
		FOR_ABP(n) {
			ind3 = n + ind * nAbp;
			VSV3ADD(&P2(rAbpVmd,ind3,0),&P2(rAbp,n,0),dimDom,ind2);	
			for (m = 0; m < 2; m++) {
				ind5 = P2A(chAbp,n,m,nChAb);
				CONT(ind5 < 0);
				ind4 = RecordConfigVmdSubroutine(&P2(rAbp,n,0), 
						&P2(rAct,ind5,0), ind2, mul);
				if (mode == 0) { 
					P2A(chAbpVmd,ind3,m,nChAb) = ind5 + ind4 * nAct;
				}
				else {
					loc[0] = FindAbpActinChain(ind5, n, 1);
					loc[0] = (loc[0] - 2) / nChAcY;
					P2A(chAbpVmd,ind3,m,nChAb) 
							= (ind4 * nAct + ind5) * nChAcX + loc[0];
				}
			}
			P2A(chAbpVmd,ind3,2,nChAb) = P2A(chAbp,n,2,nChAb);
			for (m = 3; m < 4; m++) {
				ind5 = P2A(chAbp,n,m,nChAb);
				CONT(ind5 < 0);
				ind4 = RecordConfigVmdSubroutine(&P2(rAbp,n,0), 
						&P2(rAbp,ind5,0), ind2, mul);
				P2A(chAbpVmd,ind3,m,nChAb) = ind5 + ind4 * nAbp;
			}
		}
	}
	// Fill in segment type of actin, ACP, and motor.
	for(n = 0; n < nActVmd; n++) { 
		ind = (int)(n / ((mode != 0) ? nChAcX : 1)) % nAct;
		pArr = &P2A(chAct,ind,0,nChAc);
		if (pArr[0] < 0 && pArr[1] < 0) {
			sprintf(segType[n], "G");
		}
		else if (pArr[0] > -1 && pArr[1] < 0) {
			sprintf(segType[n], "FP");
			if (mode == 1) {
				if (P2A(chActVmd,n,1,nChAcY + 2) > -1) {
					sprintf(segType[n], "F");
				}
			}
		}
		else if (pArr[0] < 0 && pArr[1] > -1) {
			sprintf(segType[n], "FB");
			if (mode == 1) {
				if (n % nChAcX > 0) {
					sprintf(segType[n], "G");
				}
			}
		}
		else { sprintf(segType[n], "F"); }
	}
	FOR_ABP(n) {
		pArr = &P2A(chAbp,n,0,nChAb);
		cnt = 1;
		if (pArr[0] < 0) { cnt++; }
		if (pArr[1] < 0) { cnt++; }
		if (motSA.gTgl != 0) {
			if (pArr[0] < 0 && pArr[1] < 0 && pArr[3] < 0 && pArr[4] < 0) {
				cnt = 4;
			}
		}
		switch(pArr[2]) {
			case 0: sprintf(segType[n + nActVmd], "CR%d", cnt);	break;
			case 1: sprintf(segType[n + nActVmd], "BU%d", cnt);	break;
			case 2: sprintf(segType[n + nActVmd], "MO%d", cnt);	break;
		}
	}
	
	if (recConfVmd.gTglInfo > 0 && currTimeStep > 1) {
		for(m = 0; m < 2; m++) {
			nRec = (m == 0) ? nAct : nAbp;
			pntRec = (m == 0) ? recActG : recAbpG;
			for(n = 0; n < nRec; n++) {
				temp = TrimDblVal(P2A(pntRec,n,2,nD), recConfVmd.minF, 
						recConfVmd.maxF);
				mass[n + ((m == 1) ? nAct : 0)] 
						= 1. + 20. * (temp  - recConfVmd.minF) 
						/ (recConfVmd.maxF - recConfVmd.minF);
			}
		}
	}
	else {
		for(n = 0; n < nAct; n++) {
			mass[n] = 12.;
		}
		for(n = 0; n < nAbp; n++) {
			mass[n + nAct] = 14.;
		}
	}
	for(n = 0; n < nAct + nAbp; n++) { occu[n] = 1.; }
	if (rho.tgl != 0) {
		for(n = 0; n < nAct; n++) {
			if (actRho[n] != -1) { occu[n] = 2.; }
		}
		for(n = 0; n < nAbp; n++) {
			if (abpRho[n] != -1) { occu[nAct + n] = 2.; }
		}
	}

	// Write PDB file
	sprintf(fn, "%s.pdb", fnIn);
	fOut = fopen(GenFileName(fn), ((recConfVmd.mode == 0) ? "w" : "a"));
	if (recConfVmd.mode != 0) {
		fprintf(fOut, "------------------------------------ %6d ------------"
				"------------------------\n", 
				(int)(currTimeStep / recConfVmd.prd));
	}
	fprintf (fOut, "REMARK\n");

	FOR_NDIM(k) {	
		domCen[k] = 0.5 * mul[k] * dimDom[k];
	}
	// actin
	for(n = 0; n < nActVmd; n++) {
		ind = (int)(n / ((mode != 0) ? nChAcX : 1)) % nAct;
	    VV3SUB(&P2(rActVmd,n,0), domCen);
		fprintf(fOut, "ATOM %6d  C   LYS     1    %8.3f%8.3f%8.3f  %4.2f  "
				"0.00      %s\n", n + 1, P2(rActVmd,n,0) * amp, 
				P2(rActVmd,n,1) * amp, P2(rActVmd,n,2) * amp, 
				occu[ind], segType[n]);
	}
	// ABP
	for(n = 0; n < nAbpVmd; n++) {
		VV3SUB(&P2(rAbpVmd,n,0), domCen);
		fprintf(fOut, "ATOM %6d  N   LYS     1    %8.3f%8.3f%8.3f  %4.2f  "
				"0.00      %s\n", n + nActVmd + 1, P2(rAbpVmd,n,0) * amp, 
				P2(rAbpVmd,n,1) * amp, P2(rAbpVmd,n,2) * amp, 
				occu[(n % nAbp) + nAct], segType[(n % nAbp) + nActVmd]);
	}
	// If a coloring method is used, a reference particle is necessary.
	if (recConfVmd.gTglInfo > 0) {
		fprintf(fOut, "ATOM %6d  O   LYS     1    %8.3f%8.3f%8.3f  %4.2f  "
				"0.00      %s\n", nActVmd + nAbpVmd + 1, 0., 0., 0., 1., "REF");
	}
	// Visualization of boundaries using cylinders
	if (recConfVmd.gTglBnd != 0) {
		for(n = 0; n < 8; n++) {
			V3IND_ASSIGN_CONST_INT(n, 2, ind2);
			V3SET_ALL(sft, 0.);
			FOR_NDIM(k) {
				rPos[k] = rGrid[k][(ind2[k] == 0) ? 0 : nGrid[k] - 1] 
						- domCen[k] + sft[k] 
						+ ((ind2[k] == 1) ? dimDom[k] * (mul[k] - 1) : 0.); 
			}
			fprintf(fOut, "ATOM %6d  H   LYS     1    %8.3f%8.3f%8.3f  %4.2f"
					"  0.00      %s\n", nAbpVmd + nActVmd + n + facInfoConfVmd 
					+ 1, rPos[0], rPos[1], rPos[2], 1., "BND");
		}
	}

	fprintf(fOut, "ENDMDL\n");
	fclose(fOut);

	// Write PSF file
	sprintf(fn, "%s.psf", fnIn);
	fOut = fopen(GenFileName(fn), ((recConfVmd.mode == 0) ? "w" : "a"));
	if (recConfVmd.mode != 0) {
		fprintf(fOut, "------------------------------------ %6d ------------"
				"------------------------\n", 
				(int)(currTimeStep / recConfVmd.prd));
	}
	fprintf(fOut, "PSF CMAP\n\n");
	fprintf(fOut, "%8d !NTITLE\n\n", 7);
	fprintf(fOut, "%8d !NATOM\n", nActVmd + nAbpVmd 
			+ facInfoConfVmd + facBndConfVmd + nBeadVmd);
	// actin
	for (n = 0; n < nActVmd; n++) {
		ind = (int)(n / ((mode != 0) ? nChAcX : 1)) % nAct;
		fprintf(fOut, "%8d %-4s 1    LYS  C      11%11.5f%14.4f%11d\n", 
				n + 1, segType[n], 0., mass[ind], 0);
	}
	// ABP
	for (n = 0; n < nAbpVmd; n++) {
		fprintf(fOut, "%8d %-4s 1    LYS  N      38%11.5f%14.4f%11d\n", 
				n + nActVmd + 1, segType[(n % nAbp) + nActVmd], 
				0.,	mass[(n % nAbp) + nAct], 0);
	}
	// A reference particle for coloring method
	if (recConfVmd.gTglInfo > 0) {
		fprintf(fOut, "%8d %-4s 1    LYS  O      38%11.5f%14.4f%11d\n", 
				nAbpVmd + nActVmd + 1, "REF", 0., 16., 0);
	}
	// Visualization of boundaries
	if (recConfVmd.gTglBnd != 0) {
		for(n = 0; n < 8; n++) {
			fprintf(fOut, "%8d %-4s 1    LYS  H      38%11.5f%14.4f%11d\n", 
					nAbpVmd + nActVmd + n + 1 + facInfoConfVmd, 
					"BND", 0., 1., 0);
		}
	}
	chain.c = 0;
	// Chain information from actin (chActVmd)
	nCh = (mode != 0) ? nChAcY + 2 : nChAc;
	for (n = 0; n < nActVmd; n++) {
		CONT(!(P2A(chActVmd,n,0,nCh) > -1 && P2A(chActVmd,n,1,nCh) < 0));
		curr = n;
		while(P2A(chActVmd,curr,0,nCh) > -1) {
			V3SUB(dr, &P2(rActVmd,curr,0), 
					&P2(rActVmd,P2A(chActVmd,curr,0,nCh),0));
			FOR_NDIM(k) {
				BREAK(pbc[k] != 0 && fabs(dr[k]) > dimDomH[k])
			}
			if (k == NDIM) { 
				P2A(chain.l,chain.c,0,2) = curr + 1;
				P2A(chain.l,chain.c,1,2) = 
						P2A(chActVmd,curr,0,nCh) + 1;
				chain.c++; 
			}
			curr = P2A(chActVmd,curr,0,nCh);	
		}
	}
	// Chain information from ABP (chAbpVmd)
	for (n = 0; n < nAbpVmd; n++) {
		pArr = &P2A(chAbpVmd,n,0,nChAb);
		for(m = 0; m < nChAb - 1; m++) {
			CONT(m == 2);
			CONT(pArr[m] < 0);
			CONT(!(ISMTF(pArr[2])) && m > 2);
			if (m < 2) { 
				V3SUB(dr, &P2(rAbpVmd,n,0), &P2(rActVmd,pArr[m],0));
			}
			else {
				V3SUB(dr, &P2(rAbpVmd,n,0),	&P2(rAbpVmd,pArr[m],0));
			}
			FOR_NDIM(k) {
				BREAK(fabs(dr[k]) > dimDomH[k])
			}
			if (k == NDIM) { 
				P2A(chain.l,chain.c,0,2) = n + nActVmd + 1;
				P2A(chain.l,chain.c,1,2) = pArr[m] + 1 
						+ ((m < 2) ? 0 : nActVmd);
				(chain.c)++; 
			}
		}
    }
	// Chain information is also required for cylinders showing boundaries.
	if (recConfVmd.gTglBnd != 0) {
		ind = nActVmd + nAbpVmd + 1 + facInfoConfVmd;
		for(n = 0; n < 12; n++) {
			P2A(chain.l,chain.c,0,2) = P2A(bndCh,n,0,2) + ind;
			P2A(chain.l,chain.c,1,2) = P2A(bndCh,n,1,2) + ind;
			(chain.c)++; 
		}	
	}
	fprintf(fOut, "\n%8d !NBOND: bonds\n", chain.c);
	for(n = 0; n < chain.c; n++) {
		fprintf(fOut, "%8d%8d", P2A(chain.l,n,0,2), P2A(chain.l,n,1,2));
		CONT((n + 1) % 4 != 0);
		fprintf(fOut, "\n"); 
	}
	if (chain.c % 4 != 0) {
		fprintf(fOut, "\n"); 
	}	
   	fclose(fOut);

	free(rActVmd);
	free(rAbpVmd);
	free(chActVmd);
	free(chAbpVmd);
	free(chain.l);
	free(occu);
	free(mass);
	free(rAct);
	free(rAbp);
	free(chAct);
	free(chAbp);
	free(iFila);
if (rho.tgl != 0) {
	free(actRho);
	free(abpRho);
}
	for (n = 0; n < nActVmd + nAbp; n++) { free(segType[n]); }
	free(segType);
	if (recConfVmd.gTglInfo > 0) {
		free(recActG);
		free(recAbpG);
	}
  }
  free(nActAll);
  free(nAbpAll);
  free(chActAlt);
}

/*-------------------------- Recording network structure ---------------------*/

/*------------------- Dynamic behaviors of actin and ABPs --------------------*/

void RecordAbpBindEvent(int abpInd, int side, int mode) {
  int n, locAbpInd, actInd, locActInd, sumCnt, nextActInd, actSide;
  int chk[26], *cntAll;
  double *recvArr, rPos[NDIM], dr[NDIM];
  FILE *fOut;  

  if (mode < 2) {
	locAbpInd = iAbp[abpInd];
	V5SET(&P2A(bindLog.l,bindLog.c,0,26), mode, (double)currTimeStep, 
			(double)abpInd, K_ABP(locAbpInd), (double)side);
	V3COPY(&P2A(bindLog.l,bindLog.c,5,26), &P2(abp.r,locAbpInd,0));
	for(n = 0; n < 2; n++) {
		actInd = P2A(abp.ch,locAbpInd,(n == 0 ? side : 1 - side),nChAb);
		if (actInd > -1) {
			locActInd = iAct[actInd];
			nextActInd = P2A(act.ch,locActInd,0,nChAc);
			actSide = FindAbpActinChain(locActInd, abpInd, 0);
			CalcPosOnActSegSide(&P2(act.r,locActInd,0), 
					&P2(act.r,iAct[nextActInd],0), rPos, actSide);
			CalcUnitVec(dr, &P2(act.r,iAct[nextActInd],0), 
					&P2(act.r,locActInd,0));
			V2SET(&P2A(bindLog.l,bindLog.c,8 + 9 * n,26), (double)actInd, 
					(double)((int)((actSide - 2) / nChAcY)) / (double)nChAcX);
			V3COPY(&P2A(bindLog.l,bindLog.c,10 + 9 * n,26), rPos);
			V3COPY(&P2A(bindLog.l,bindLog.c,13 + 9 * n,26), dr);
			P2A(bindLog.l,bindLog.c,16 + 9 * n,26) = (double)act.iF[locActInd];
		}
		else {
			SetAllValue1dArrayDouble(&P2A(bindLog.l,bindLog.c,8 + 9 * n,26), 
					9, 0.);
		}
	}
	(bindLog.c)++;
  }
  else {
	MALLOC(cntAll,int,nCpu); 
	MPI_Gather(&bindLog.c, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		sumCnt = SumArrInt(cntAll, nCpu);
		MALLOC(recvArr,double,sumCnt * 26);
	}
	Gather2dArrayDoubleWoIndWoSort(bindLog.c, cntAll, bindLog.l, recvArr, 26);
	if (rank == 0) {
		qsort(recvArr, sumCnt, 26 * sizeof(double), CompDbl);
		fOut = fopen(GenFileName("AbpBind"), "a");
		memset(chk, 0, sizeof(int) * 26);
		chk[0] = chk[2] = chk[3] = chk[8] = chk[16] = chk[17] = chk[25] = 1;
		Fprintf2dArrayIntDouble(fOut, recvArr, sumCnt, 26, chk);
		fclose(fOut);
		free(recvArr);
	}
	free(cntAll);
	bindLog.c = 0;
  }
}

void RecordAbpUnbindEvent(int abpInd, int side, int mode) {
  int n, locAbpInd, actInd, locActInd, sumCnt, nextActInd, actSide;
  int chk[29], *cntAll;
  double *recvArr, rPos[NDIM], dr[NDIM];
  FILE *fOut;  

  double len, len2, dr2[NDIM], drAxis[NDIM], fac, rPos2[NDIM], fInst[2];
  double fi[NDIM], f, f2;
  int kind, k;

  if (mode == 0) {
	locAbpInd = iAbp[abpInd];
	V4SET(&P2A(unbLog.l,unbLog.c,0,29), 0, (double)currTimeStep, (double)abpInd,
			(double)side);
	V3COPY(&P2A(unbLog.l,unbLog.c,4,29), &P2(abp.r,locAbpInd,0));
	for(n = 0; n < 2; n++) {
		actInd = P2A(abp.ch,locAbpInd,(n == 0 ? side : 1 - side),nChAb);
		if (actInd > -1) {
			locActInd = iAct[actInd];
			nextActInd = P2A(act.ch,locActInd,0,nChAc);
			actSide = FindAbpActinChain(locActInd, abpInd, 0);
			CalcPosOnActSegSide(&P2(act.r,locActInd,0), 
					&P2(act.r,iAct[nextActInd],0), rPos, actSide);
			CalcUnitVec(dr, &P2(act.r,iAct[nextActInd],0), 
					&P2(act.r,locActInd,0));
			V2SET(&P2A(unbLog.l,unbLog.c,7 + 11 * n,29), (double)actInd, 
					(double)((int)((actSide - 2) / nChAcY)) / (double)nChAcX);
			V3COPY(&P2A(unbLog.l,unbLog.c,9 + 11 * n,29), rPos);
			V3COPY(&P2A(unbLog.l,unbLog.c,12 + 11 * n,29), dr);
			V2COPY(&P2A(unbLog.l,unbLog.c,15 + 11 * n,29), 
					&P2A(recInstSprFabp,locAbpInd,side * 2,4));
			P2A(unbLog.l,unbLog.c,17 + 11 * n,29) = (double)act.iF[locActInd];
		}
		else {
			SetAllValue1dArrayDouble(&P2A(unbLog.l,unbLog.c,7 + 11 * n,29), 
					11, 0.);
		}
	}
	(unbLog.c)++;
  }
  else if (mode == 1) {
	locAbpInd = iAbp[abpInd];
	V4SET(&P2A(unbLog.l,unbLog.c,0,29), 1, (double)currTimeStep, (double)abpInd,
			(double)side);
	V3COPY(&P2A(unbLog.l,unbLog.c,4,29), &P2(abp.r,locAbpInd,0));
	for(n = 0; n < 2; n++) {
		actInd = P2A(abp.ch,locAbpInd,(n == 0 ? side : 1 - side),nChAb);
		if (actInd > -1) {
			locActInd = iAct[actInd];
			nextActInd = P2A(act.ch,locActInd,0,nChAc);
			kind = K_ABP(locAbpInd);

			// Find the location of ABP on the actin segment.
			actSide = FindAbpActinChain(locActInd, abpInd, 0);
			CalcPosOnActSegSide(&P2(act.r,locActInd,0), 
					&P2(act.r,iAct[nextActInd],0), rPos, actSide);
			len = CalcVecDist(dr, rPos, &P2(abp.r,locAbpInd,0), 0);
			if (ISMTF(kind)) {
				CalcVec(drAxis, &P2(act.r,locActInd,0), 
						&P2(act.r,iAct[nextActInd],0));
				fac = V3DOT(dr, drAxis) / V3LEN_SQ(drAxis);
				VS3SUB(rPos2, rPos, drAxis, fac);
				ApplyBoundCondVector(rPos2, -1, 0);
				len = CalcVecDist(dr, rPos2, &P2(abp.r,locAbpInd,0), 0);
				len2 = CalcVecDist(dr2, rPos, rPos2, 0);
			}	
			f = SPRING(abpF.spr[kind].stf, len, abpF.spr[kind].eq);
			if (ISMTF(kind)) {
				f2 = SPRING(abpF.spr[kind].stf2, len2, 0.);
			}

			f /= len;
			if (ISMTF(kind) && len2 > 0) { f2 /= len2; }
			FOR_NDIM(k) {
				fi[k] = f * dr[k];
				if (ISMTF(kind)) { fi[k] += f2 * dr2[k]; }
			}
			fInst[0] = REVSIGN(f * len);
			if (ISMTF(kind)) {
				fInst[1] = f2 * len2 * ((fac < 0) ? -1. : 1.);
			}
			else {
				if (f < 0) {
					CalcUnitVec(drAxis, &P2(act.r,locActInd,0),
							&P2(act.r,iAct[nextActInd],0));
					fInst[1] = V3DOT(fi, drAxis);
				}
				else {
					fInst[1] = 0;
				}
			}
			CalcUnitVec(drAxis, &P2(act.r,iAct[nextActInd],0), 
					&P2(act.r,locActInd,0));
			V2SET(&P2A(unbLog.l,unbLog.c,7 + 11 * n,29), (double)actInd, 
					(double)((int)((actSide - 2) / nChAcY)) / (double)nChAcX);
			V3COPY(&P2A(unbLog.l,unbLog.c,9 + 11 * n,29), rPos);
			V3COPY(&P2A(unbLog.l,unbLog.c,12 + 11 * n,29), drAxis);
			V2COPY(&P2A(unbLog.l,unbLog.c,15 + 11 * n,29), fInst);
			P2A(unbLog.l,unbLog.c,17 + 11 * n,29) = (double)act.iF[locActInd];
		}
		else {
			SetAllValue1dArrayDouble(&P2A(unbLog.l,unbLog.c,7 + 11 * n,29), 
					11, 0.);
		}
	}
  }
  else {
	MALLOC(cntAll,int,nCpu); 
	MPI_Gather(&unbLog.c, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		sumCnt = SumArrInt(cntAll, nCpu);
		MALLOC(recvArr,double,sumCnt * 29);
	}
	Gather2dArrayDoubleWoIndWoSort(unbLog.c, cntAll, unbLog.l, recvArr, 29);
	if (rank == 0) {
		qsort(recvArr, sumCnt, 29 * sizeof(double), CompDbl);
		fOut = fopen(GenFileName("AbpUnb"), "a");
		memset(chk, 0, sizeof(int) * 29);
		chk[0] = chk[2] = chk[3] = chk[7] = chk[17] = chk[18] = chk[28] = 1;
		Fprintf2dArrayIntDouble(fOut, recvArr, sumCnt, 29, chk);
		fclose(fOut);
		free(recvArr);
	}
	free(cntAll);
	unbLog.c = 0;
  }
}

// Record the turnover of ABPs: to find the displacement from the location 
// of unbinding to the location of following binding and whether unbound ABPs
// bind to the same filament or not.
// [Columns of abpTurn]
// 0: the moment of unbinding 
// 1, 2, 3: previous position, 4: force applied at the moment of unbinding
// 5: previous shear strain, 6: previous filament number 
// mode = 0: called at the moment of unbinding, 1: called at that of binding
void RecordAbpTurnover(int abpInd, int actInd, int side, int mode) {
  int locAbpInd, sameFila, *cntAll, sumCnt, chk[14];
  double dr[NDIM], len, *recvArr;
  FILE *fOut;

  locAbpInd = iAbp[abpInd];
  // At the moment of unbinding
  if (mode == 0) {
	P2A(abpTurn,locAbpInd,0,7) = (double)currTimeStep;
	V3COPY(&P2A(abpTurn,locAbpInd,1,7), &P2(abp.r,locAbpInd,0));
	V3SET(&P2A(abpTurn,locAbpInd,4,7), 0., 
			P2A(recInstSprFabp,locAbpInd,side * 2,4),
			(double)act.iF[iAct[actInd]]);
  }
  // At the moment of binding
  else if (mode == 1) {
	V3SUB(dr, &P2(abp.r,locAbpInd,0), &P2A(abpTurn,locAbpInd,1,7));
    ApplyBoundCondVecDiff(dr);
	len = V3LEN(dr);
	sameFila = (act.iF[iAct[actInd]] 
			== (int)P2A(abpTurn,locAbpInd,6,7)) ? 1 : 0;
	V6SET(&P2A(toLog.l,toLog.c,0,14), P2A(abpTurn,locAbpInd,0,7), 
			(double)currTimeStep, (double)abpInd, 
			F_S2PN(P2A(abpTurn,locAbpInd,5,7)), L_S2NM(len), (double)sameFila);
	V3COPY(&P2A(toLog.l,toLog.c,6,14), &P2(abp.r,locAbpInd,0));
	V4COPY(&P2A(toLog.l,toLog.c,9,14), &P2A(abpTurn,locAbpInd,1,7));
	P2A(toLog.l,toLog.c,13,14) = 0.;
	(toLog.c)++;
	SetAllValue1dArrayDouble(&P2A(abpTurn,locAbpInd,0,7), 7, -1.);
  }
  else {
	MALLOC(cntAll,int,nCpu); 
	MPI_Gather(&toLog.c, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		sumCnt = SumArrInt(cntAll, nCpu);
		MALLOC(recvArr,double,sumCnt*14);
	}
	Gather2dArrayDoubleWoIndWoSort(toLog.c, cntAll, toLog.l, recvArr, 14);
	if (rank == 0) {
		qsort(recvArr, sumCnt, 14 * sizeof(double), CompDbl);
		fOut = fopen(GenFileName("AbpTurnover"), "a");
		memset(chk, 0, sizeof(int) * 14);
		chk[0] = chk[1] = chk[2] = chk[5] = 1;
		Fprintf2dArrayIntDouble(fOut, recvArr, sumCnt, 14, chk);
		fclose(fOut);
		free(recvArr);
	}
	free(cntAll);
	toLog.c = 0;
  }
}

/*------------------- Dynamic behaviors of actin and ABPs --------------------*/

/*--------------------- Structural properties of network ---------------------*/

// Record the distribution of actin filaments
void RecordFilamentLength(int period) {
  int n, len, curr, cntFila, maxLen, *filaLen, *nActAll, *pArr;
  double avgFilaLen, stdFilaLen, len2;
  FILE *fOut;

  maxLen = -1;
  MALLOC(nActAll,int,nCpu);
  if (rank == 0) { MALLOC(chAct,int,nAct*nChAc); }
  // Gather information about chains.
  MPI_Gather(&nActMe, 1, MPI_INT, nActAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  Gather2dArrayInt(nActMe, nActAll, nChAc, act.id, act.ch, chAct); 

  if (rank == 0) {
	MALLOC(filaLen,int,nAct+1);
	memset(filaLen, 0, sizeof(int) * (nAct + 1) );
	cntFila = 0; 
	avgFilaLen = 0.; 
	stdFilaLen = 0.;

	FOR_ACT(n) {
		pArr = &P2A(chAct,n,0,nChAc);
		if (pArr[0] < 0 && pArr[1] < 0) { 
			filaLen[0]++; 
		} 
		else if (pArr[0] > -1 && pArr[1] < 0) { 
			len = 1; 
			curr = pArr[0];
			while (P2A(chAct,curr,0,nChAc) > -1) {
				curr = P2A(chAct,curr,0,nChAc);
				len++;
			}
			filaLen[len]++;
			if (maxLen < len) { maxLen = len; }
		}
	}
	fOut = fopen(GenFileName("ActFilaLen"), "a");
	// Record the number of lines recorded.
    Fprintf1dArrayIntWFil(fOut, &maxLen, 1, 0, 3);
	// Calculate average length and standard deviation.
	for(n = 1; n < maxLen + 1; n++) {
		len2 = L_S2UM(n);
	    fprintf(fOut, "%d\t%g\t%d\n", n - 1, len2, filaLen[n]);
	    cntFila += filaLen[n];
		avgFilaLen += (double)(len2 * filaLen[n]);
		stdFilaLen += (double)(SQR(len2) * filaLen[n]);
	}
	if (cntFila == 0) { 
		avgFilaLen = 0.; 
		stdFilaLen = 0.;
	}
	else { 
		avgFilaLen /= (double)cntFila; 
		stdFilaLen = sqrt(stdFilaLen / (double)cntFila - SQR(avgFilaLen));
	}
	fprintf(fOut, "%d\t%g\t%g\n", maxLen - 1, avgFilaLen, stdFilaLen);
	fclose(fOut);
	free(filaLen);
	free(chAct);
  }
  free(nActAll);
}

/*--------------------- Structural properties of network ---------------------*/

/*----------------------------- Energy calculation ---------------------------*/

// mode = 0: measure energy of entire network
// mode = 1: measure energy of selected portion of the network
//			 to the other 
void RecordMechEnergy(int mode) {
  int n, k, locActInd, *pArr, side, kind;
  double len, lenEq, ang, fac;
  double E[4], Esum[4], *Eall, dr[2][NDIM + 1];
  double rPos[NDIM], rPos2[NDIM], *rPnt[2];
  FILE *fOut;
 
  if (rank == 0) { MALLOC(Eall,double,nCpu*4); }
  V4SET_ALL(E, 0.);

  // Bending of actin filament
  FOR_ACTME(n) {
	pArr = &P2A(act.ch,n,0,nChAc);
    CONT(!(pArr[0] > -1 && pArr[1] > -1));
	CalcVec(dr[0], &P2(act.r,n,0), &P2(act.r,iAct[pArr[0]],0));
	CalcVec(dr[1], &P2(act.r,iAct[pArr[1]],0), &P2(act.r,n,0));
	ang = V3ANG(dr[0], dr[1]);
	E[0] += 0.5 * actF.bend.stf * SQR(ang);
  }

  // Extension of actin filament (double-counting)
  FOR_ACTME(n) {
	pArr = &P2A(act.ch,n,0,nChAc);
	CONT(pArr[0] < 0);
	len = CalcDist(&P2(act.r,n,0), &P2(act.r,iAct[pArr[0]],0), 0);
	E[1] += 0.5 * actF.spr.stf * SQR(len - actF.spr.eq);
  }

  // Bending of ABP (actin-ABP-actin)
  FOR_ABPME(n) {
	pArr = &P2A(abp.ch,n,0,nChAb);
	CONT(!(pArr[0] > -1 && pArr[1] > -1));
	kind = pArr[2];
	CONT(ISMTF(kind));
	for(k = 0; k < 2; k++) {
		dr[k][NDIM] = CalcVecDistActinAbp(dr[k], pArr[k], abp.id[n], 0);
	}
	V3REVSIGN(dr[1]);
	ang = V3ANG(dr[0], dr[1]);
	E[2] += 0.5 * abpF.bend[kind].stf * SQR(ang - abpF.bend[kind].eq);
  }
  // Bending of motor-assembled structures
  if (motSA.gTgl != 0) {
	for(n = 0; n < nAbpMe; n++) {
		pArr = &P2A(abp.ch,n,0,nChAb);
		CONT(pArr[2] != 2);
	    CONT(!(pArr[3] > -1 && pArr[4] > -1));
		CalcVec(dr[0], &P2(abp.r,n,0), &P2(abp.r,iAbp[pArr[3]],0));
		CalcVec(dr[1], &P2(abp.r,iAbp[pArr[4]],0), &P2(abp.r,n,0));
		ang = V3ANG(dr[0], dr[1]);
		E[2] += 0.5 * motSA.bend.stf * SQR(ang);
	}
  }

  // Extension of ABP
  FOR_ABPME(n){
	CONT(ISABPM(n));
	pArr = &P2A(abp.ch,n,0,nChAb);
	kind = pArr[2];
	for(k = 0; k < 2; k++) {
		CONT(pArr[k] < 0);
		locActInd = iAct[pArr[k]];
		side = FindAbpActinChain(locActInd, abp.id[n], 0);
		rPnt[0] = &P2(act.r,locActInd,0);
		rPnt[1] = &P2(act.r,iAct[P2A(act.ch,locActInd,0,nChAc)],0);
		CalcPosOnActSegSide(rPnt[0], rPnt[1], rPos, side);
  		CalcVec(dr[0], rPos, &P2(abp.r,n,0));
		if (ISMTF(kind)) {
			CalcVec(dr[1], rPnt[0], rPnt[1]);
			fac = V3DOT(dr[0], dr[1]) / V3LEN_SQ(dr[1]);
			VS3SUB(rPos2, rPos, dr[1], fac);
			ApplyBoundCondVector(rPos2, -1, 0);
			CalcVec(dr[0], rPos2, &P2(abp.r,n,0));
			CalcVec(dr[1], rPos, rPos2);			
		}
		len = V3LEN(dr[0]);
		E[3] += 0.5 * abpF.spr[kind].stf 
				* SQR(len - abpF.spr[kind].eq);
		if (ISMTF(kind)) {	
			len = V3LEN(dr[1]);
			E[3] += 0.5 * abpF.spr[kind].stf2 * SQR(len);
		}
	}
  }
  
  // Extension of motor-assembled structures
  if (motSA.gTgl != 0) {
	for(n = 0; n < nAbpMe; n++) {
		pArr = &P2A(abp.ch,n,0,nChAb);
		CONT(pArr[2] != 2);
		CONT(pArr[3] < 0);
		len = CalcDist(&P2(abp.r,n,0), &P2(abp.r,iAbp[pArr[3]],0), 0);
		lenEq = (abp.mId[n] == 0 && abp.mId[iAbp[pArr[3]]] == 0) 
				? motSA.cenDist : motSA.spr.eq;
		E[3] += 0.5 * motSA.spr.stf * SQR(len - lenEq);
	}
  }

  MPI_Gather(E, 4, MPI_DOUBLE, Eall, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    V4SET_ALL(Esum, 0.); 
    for(n = 0; n < nCpu; n++) {
        for(k = 0; k < 4; k++) {
            Esum[k] += Eall[4 * n + k];
        }
    }
	if (mode == 0) {
	    fOut = fopen(GenFileName("MechEall"), "a");
	}
	else if (mode == 1) {
	    fOut = fopen(GenFileName("MechEperc"), "a");
	}
	else {
	    fOut = fopen(GenFileName("MechEsupp"), "a");
	}
    fprintf(fOut, "%lld\t", currTimeStep);
	Fprintf1dArrayDouble(fOut, Esum, 4, 0);
    fclose(fOut);
    free(Eall);
  }
}

/*----------------------------- Energy calculation ---------------------------*/

/*--------------------------- Accumulated information ------------------------*/

// Accumulate forces acting on actin and ABP. This information is not recorded 
// to file, but used for coloring methods via VMD.
void RecordAccuLengthForces(void) {
  int n, k, cntCh;
  int actInd, nextActInd, nextAbpInd;
  double lenAccu, len;
  // actin
  FOR_ACTME(n) {
	CONT(ISACTM(n));
	VV3ADD(&P2A(recAct.allF,n,0,NDIM + 1), &P2(act.f,n,0));
	P2A(recAct.allF,n,NDIM,NDIM + 1) += V3LEN(&P2(act.f,n,0));
	nextActInd = P2A(act.ch,n,0,nChAc);
	if (nextActInd > -1) {
		len = CalcDist(&P2(act.r,n,0), &P2(act.r,iAct[nextActInd],0), 0);
		recAct.len[n] += len;
		recAct.cnt[n] += 1;
	}
	else {
		recAct.len[n] = 0.;
		recAct.cnt[n] = 0;
	}
  } 
  // ABP
  FOR_ABPME(n) {
	CONT(ISABPM(n));
	VV3ADD(&P2A(recAbp.allF,n,0,NDIM + 1), &P2(abp.f,n,0));
	P2A(recAbp.allF,n,NDIM,NDIM + 1) += V3LEN(&P2(abp.f,n,0));
	cntCh = 0;
	lenAccu = 0.;
	for(k = 0; k < 2; k++) {
		actInd = P2A(abp.ch,n,k,nChAb);
		CONT(actInd < 0);
		len = CalcDistActinAbp(actInd, abp.id[n], 0);
		lenAccu += len / abpF.spr[K_ABP(n)].eq;
		cntCh++;
	}
	if (cntCh > 0) {
		P2A(recAbp.len,n,0,recAbp.nL) += lenAccu / (double)cntCh;
		recAbp.cnt[n] += 1;
	}
	else {
		P2A(recAbp.len,n,0,recAbp.nL) = 0.;
		if (!(ISMTF(K_ABP(n)))) {
			recAbp.cnt[n] = 0;
		}
	}
	CONT(!(ISMTF(K_ABP(n))));
	for(k = 0; k < 2; k++) {
		nextAbpInd = P2A(abp.ch,n,k + 3,nChAb);
		if (nextAbpInd > -1) {
			len = CalcDist(&P2(abp.r,n,0), &P2(abp.r,iAbp[nextAbpInd],0), 0);
			P2A(recAbp.len,n,1,recAbp.nL) += len;
			if (cntCh == 0) { 
				recAbp.cnt[n] += 1;
			}
		}
		else {
			P2A(recAbp.len,n,1,recAbp.nL) = 0.;
			if (cntCh == 0) {
				recAbp.cnt[n] = 0;
			}
		}
	}
  }
}

void ResetAccuLengthForces(void) {
  int n;

  FOR_ACTC(n) {
	recAct.sprF[n] = 0.;
	recAct.bendF[n] = 0.;
	V4SET_ALL(&P2A(recAct.allF,n,0,NDIM + 1), 0.);
	recAct.len[n] = 0.;
	recAct.cnt[n] = 0.;
  }
  FOR_ABPC(n) {
	recAbp.sprF[n] = 0.;
	recAbp.bendF[n] = 0.;
	V4SET_ALL(&P2A(recAbp.allF,n,0,NDIM + 1), 0.);
	SetAllValue1dArrayDouble(&P2A(recAbp.len,n,0,recAbp.nL), recAbp.nL, 0.);
	recAbp.cnt[n] = 0.;
  }
}

// Calculate and record the longitudinal forces acting on ACP or motor. This
// information is very useful to analyze the activity of motors since 
// the longitudinal forces are known to affect the rate of walking events.
void RecordAccuLongSpringForces(int period) {
  int n, k, *nAbpAll;
  double *f, *fAbp;
  FILE *fOut;

  // Accmulate the calculated longitudinal forces.
  FOR_ABPME(n) {
	CONT(ISABPM(n));
    for (k = 0; k < 2; k++) {
		CONT(P2A(abp.ch,n,k,nChAb) < 0);
        // find force
        f = &P2A(recInstSprFabp,n,k * 2,4);
		P2A(recLongSprFabp,n * 2 + k,0,2) += f[0];
		P2A(recLongSprFabp,n * 2 + k,1,2) += f[1];
	}
  }
  if (currTimeStep % period == 0) {
	MALLOC(nAbpAll,int,nCpu);
	if (rank == 0) { MALLOC(fAbp,double,nAbp * 2 * 2); }
	MPI_Gather(&nAbpMe, 1, MPI_INT, nAbpAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
	Gather2dArrayDouble(nAbpMe, nAbpAll, 4, abp.id, recLongSprFabp, fAbp);	
    for(n = 0; n < nAbpC * 4; n++) {
        recLongSprFabp[n] = 0.;
    }

	if (rank == 0) {
		for(n = 0; n < nAbp * 4; n++) {
			fAbp[n] = F_S2PN(fAbp[n]) / (double)period;
		}
		fOut = fopen(GenFileName("AbpLongFor"), "a");
		if (currTimeStep / period == 1) {
    		Fprintf1dArrayIntWFil(fOut, &nAbp, 1, 0, 5);
		}
		Fprintf2dArrayDouble(fOut, fAbp, nAbp, NDIM + 1, 0, 0);
		fclose(fOut);
		free(fAbp);
	}
	free(nAbpAll);
  }
}

/*--------------------------- Accumulated information ------------------------*/

/*------------------ Information for filaments and segments  -----------------*/

void RecordIndvSegFilaInformation(int period) {
  int n, k, l, curr, ind, iPercSum, iSuppSum, cnt, cntFor, CS, side;
  int abpInd, actInd, prevAbpInd, prevActInd, sft[NDIM], cntXlink[3];
  int *iFila, *nActAll, *nAbpAll, *cntAct, *cntAbp, *iMot;
  int *binAct, *binAbp, *binActAll, *binAbpAll, *pArr, *pArr2;
  double len, lenCtr, lenEE, lenSeg, ang, angSum[NDIM], crsAng[2];
  double tensFactSum, tensFabpSum, bendFactSum, bendFabpSum;
  double dr[NDIM], dr2[NDIM], drCrs[NDIM], drSum[NDIM];
  double viscFsum[NDIM], rPos[2][NDIM], ratio[2];
  double *fAct, *fAbp, *fBendAct, *fBendAbp;
  double *lAct, *lAbp, *tensF, *bendF, *viscF;
  FILE *fOut;
  int *actRho, *abpRho;

  MALLOC(nActAll,int,nCpu);
  MALLOC(nAbpAll,int,nCpu);
  MALLOC(binAct,int,nActC);
  MALLOC(binAbp,int,nAbpC);
  MPI_Gather(&nActMe, 1, MPI_INT, nActAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&nAbpMe, 1, MPI_INT, nAbpAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	MALLOC(fAct,double,nAct * (NDIM + 1));
	MALLOC(fAbp,double,nAbp * (NDIM + 1));
	MALLOC(fBendAct,double,nAct);
	MALLOC(fBendAbp,double,nAbp);
	MALLOC(rAct,double,nAct * NDIM);
	MALLOC(rAbp,double,nAbp * NDIM);
	MALLOC(lAct,double,nAct);
	MALLOC(lAbp,double,nAbp * recAbp.nL);
	MALLOC(chAct,int,nAct*nChAc);
	MALLOC(chAbp,int,nAbp*nChAb);
	MALLOC(cntAct,int,nAct);
	MALLOC(cntAbp,int,nAbp);
	MALLOC(tensF,double,nAct + nAbp);
	MALLOC(bendF,double,nAct + nAbp);
	MALLOC(viscF,double,(nAct + nAbp) * NDIM);
	MALLOC(fixAct,int,nAct);
	MALLOC(iFila,int,nAct);
	MALLOC(binActAll,int,nAct);
	MALLOC(binAbpAll,int,nAbp);
	if (motSA.gTgl != 0) {
		MALLOC(iMot,int,nAbp);
	}
	if (rho.tgl != 0) {
		MALLOC(actRho,int,nAct);
		MALLOC(abpRho,int,nAbp);
	}
  }
  FOR_ACTME(n) { 
	binAct[n] = BinaryPackActinChainArray(act.id[n]);
  }
  FOR_ABPME(n) { 
	binAbp[n] = BinaryPackAbpChainArray(abp.id[n]);
  }
  GatherActChainPosition(nActAll, chAct, rAct);
  GatherAbpChainPosition(nAbpAll, chAbp, rAbp);
  Gather2dArrayDouble(nActMe, nActAll, NDIM + 1, act.id, recAct.allF, fAct);
  Gather2dArrayDouble(nAbpMe, nAbpAll, NDIM + 1, abp.id, recAbp.allF, fAbp);
  Gather2dArrayDouble(nActMe, nActAll, 1, act.id, recAct.bendF, fBendAct);
  Gather2dArrayDouble(nAbpMe, nAbpAll, 1, abp.id, recAbp.bendF, fBendAbp);
  Gather2dArrayDouble(nActMe, nActAll, 1, act.id, recAct.len, lAct);
  Gather2dArrayDouble(nAbpMe, nAbpAll, recAbp.nL, abp.id, recAbp.len, lAbp);
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, recAct.cnt, cntAct);
  Gather2dArrayInt(nAbpMe, nAbpAll, 1, abp.id, recAbp.cnt, cntAbp);
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, act.iF, iFila);
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, act.fix, fixAct);
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, binAct, binActAll);
  Gather2dArrayInt(nAbpMe, nAbpAll, 1, abp.id, binAbp, binAbpAll);
  if (rho.tgl != 0) {
	Gather2dArrayInt(nActMe, nActAll, 1, act.id, act.rho, actRho);  
	Gather2dArrayInt(nAbpMe, nAbpAll, 1, abp.id, abp.rho, abpRho);  
  }

  if (rank == 0) {
	if (motSA.gTgl != 0) {
		memset(iMot, -1, sizeof(int) * nAbp);
		cnt = 0;
		FOR_ABP(n) {
			pArr = &P2A(chAbp,n,0,nChAb);
			CONT(!(pArr[3] > -1 && pArr[4] < 0));
			curr = n;
			while(P2A(chAbp,curr,3,nChAb) > -1) {
				iMot[curr] = cnt;
				curr = P2A(chAbp,curr,3,nChAb);
			}
			iMot[curr] = cnt;
			cnt++;
		}
	}
	fOut = fopen(GenFileName("InfoIndv"), "a");
	if (currTimeStep / period == 1) {
		fprintf(fOut, "%d\t%d\t%d\t", nAct + nAbp, nAct, nAbp);
		Fprintf1dFillerInt(fOut, 0, 17, 0);
	}
	FOR_ACT(n) {
		pArr = &P2A(chAct,n,0,nChAc);
		if ((pArr[0] < 0 && pArr[1] < 0) || cntAct[n] == 0) {
			fprintf(fOut, "%d\t", n);
			Fprintf1dFillerInt(fOut, 0, 11, 1);
			Fprintf1dFillerInt(fOut, -1, 8, 0);
			tensF[n] =	0.;
			bendF[n] =	0.;
			V3SET_ALL(&P2(viscF,n,0), 0.);
			continue;
		}
		cnt = 0;
		V3SET_ALL(drSum, 0.);
		for(k = 0; k < 2; k++) {
			CONT(pArr[k] < 0);
			CalcUnitVec(dr, &P2(rAct,n,0), &P2(rAct,pArr[k],0));
			if (k == 1) { V3REVSIGN(dr); }
			VV3ADD(drSum, dr);
			cnt++;
		}
		NormVec(drSum);
		tensF[n] = F_S2PN(actF.spr.stf * (lAct[n] / (double)cntAct[n] 
				- actF.spr.eq));
		bendF[n] = F_S2PN(fBendAct[n] / (double)cntAct[n]);
		VS3COPY(&P2(viscF,n,0), &P2A(fAct,n,0,NDIM + 1), INV(cntAct[n]));
		FOR_NDIM(k) { P2(viscF,n,k) = F_S2PN(P2(viscF,n,k)); }
		fprintf(fOut, "%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t"
				"%g\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", n, P2(rAct,n,0), 
				P2(rAct,n,1), P2(rAct,n,2), drSum[0], drSum[1], drSum[2], 
				tensF[n], bendF[n], P2(viscF,n,0), P2(viscF,n,1), P2(viscF,n,2),
				pArr[0], pArr[1], 0, 0, fixAct[n], iFila[n], 
				binActAll[n], actRho[n]);
	}
	FOR_ABP(n) {
		pArr = &P2A(chAbp,n,0,nChAb);
		if ((ISMTF(pArr[2]) && pArr[0] < 0 && pArr[1] < 0 && pArr[3] < 0 
				&& pArr[4] < 0) || (!(ISMTF(pArr[2])) && pArr[0] < 0 
				&& pArr[1] < 0) || cntAbp[n] == 0) {
			fprintf(fOut, "%d\t", n + nAct);
			Fprintf1dFillerInt(fOut, 0, 11, 1);
			Fprintf1dFillerInt(fOut, -1, 8, 0);
			tensF[n + nAct] = 0.;
			bendF[n + nAct] = 0.;
			V3SET_ALL(&P2(viscF,n + nAct,0), 0.);
			continue;
		}
		cnt = 0;
		V3SET_ALL(drSum, 0.);
		ind = (ISMTF(pArr[2])) ? 3 : 0;
		for(k = 0; k < 2; k++) {
			CONT(pArr[k + ind] < 0);
			if (ind == 0) {
				CalcUnitVecActinAbp(dr, pArr[k], n, 1);
			}
			else {
				CalcUnitVec(dr, &P2(rAbp,pArr[k + ind],0), &P2(rAbp,n,0)); 
			}
			if (k == 1) { V3REVSIGN(dr); }
			VV3ADD(drSum, dr);
			cnt++;
		}
		NormVec(drSum);
		ind = n + nAct;
		tensF[ind] = F_S2PN(abpF.spr[pArr[2]].stf * abpF.spr[pArr[2]].eq 
				* (P2A(lAbp,n,0,recAbp.nL) / (double)cntAbp[n] - 1.));
		bendF[ind] = F_S2PN(fBendAbp[n] / (double)cntAbp[n]);
		VS3COPY(&P2(viscF,ind,0), &P2A(fAbp,n,0,NDIM + 1), INV(cntAbp[n]));
		FOR_NDIM(k) { P2(viscF,ind,k) = F_S2PN(P2(viscF,ind,k)); }
		for(k = 0; k < 2; k++) {
			if (pArr[k] > -1) {
				side = FindAbpActinChain(pArr[k], n, 1);
				ratio[k] = (side > -1) ? (double)((side - 2) / nChAcY) 
						/ (double)nChAcX : -1.;
			}
			else {
				ratio[k] = -1.;
			}
		}
		fprintf(fOut, "%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t"
				"%g\t%d\t%d\t%d\t%d\t%g\t%g\t%d\t%d\n", ind, P2(rAbp,n,0), 
				P2(rAbp,n,1), P2(rAbp,n,2), drSum[0], drSum[1], drSum[2], 
				tensF[ind], bendF[ind], P2(viscF,ind,0), P2(viscF,ind,1), 
				P2(viscF,ind,2), pArr[0], pArr[1], 0, 0, 
				ratio[0], ratio[1],	binAbpAll[n], abpRho[n]);
	}
	fclose(fOut);

	// Record information in unit of filaments
	fOut = fopen(GenFileName("InfoFila"), "a");
	// Count and print the number of actin filaments
	cnt = 0;
	FOR_ACT(n) {
		CONT(!(P2A(chAct,n,0,nChAc) > -1 && P2A(chAct,n,1,nChAc) < 0));
		cnt++;
	}
    Fprintf1dArrayIntWFil(fOut, &cnt, 1, 0, 26);
	// Record the information
	cnt = 0;
	FOR_ACT(n) {
		CONT(!(P2A(chAct,n,0,nChAc) > -1 && P2A(chAct,n,1,nChAc) < 0));
		tensFactSum = 0.;
		tensFabpSum = 0.;
		bendFactSum = 0.;
		bendFabpSum = 0.;
		lenSeg = 0.;
		lenCtr = 0.;
		iPercSum = 0;
		iSuppSum = 0;
		V3SET_ALL(angSum, 0.);
		V3SET_ALL(cntXlink, 0);
		V3SET_ALL(drSum, 0.);
		V3SET_ALL(viscFsum, 0.);
		V3SET_ALL(sft, 0);
		curr = n;
		while(P2A(chAct,curr,0,nChAc) > -1) {
			pArr = &P2A(chAct,curr,0,nChAc);
			len = CalcVecDist(dr, &P2(rAct,pArr[0],0), &P2(rAct,curr,0), 0);
			lenCtr += len;
			CheckCrossBound(sft, &P2(rAct,curr,0), &P2(rAct,pArr[0],0));
			if (pArr[1] > -1) {
				// Calculate bending angles
				CalcVec(dr2, &P2(rAct,curr,0), &P2(rAct,pArr[1],0));
				V3CROSS(drCrs, dr, dr2);
				NormVec(drCrs);
				ang = V3ANG(dr, dr2);
				VVS3ADD(angSum, drCrs, ang);
			}
			VVS3ADD(drSum, dr, INV(len));
			tensFactSum += tensF[curr];
			bendFactSum += bendF[curr];
			for(k = 2; k < nChAc; k++) {
				abpInd = pArr[k];
				CONT(abpInd < 0);
				pArr2 = &P2A(chAbp,abpInd,0,nChAb);
				CONT(!(ISMTF(pArr2[2])) && pArr2[1] < 0);
				tensFabpSum += tensF[nAct + abpInd];
				bendFabpSum += bendF[nAct + abpInd];
				iSuppSum += 0;
				cntXlink[(int)(pArr2[2] / 2)]++; 
				CONT(ISMTF(pArr2[2]));
				actInd = pArr2[(pArr2[0] == curr) ? 1 : 0];
				CalcVec(dr2, &P2(rAct,P2A(chAct,actInd,0,nChAc),0), 
						&P2(rAct,actInd,0));
				if (V3DOT(dr, dr2) > 0.) { cntXlink[2]++; }
			}
			VV3ADD(viscFsum, &P2(viscF,curr,0));
			iPercSum += 0.;
			iSuppSum += 0;
			lenSeg += 1.;
			curr = pArr[0];
		}
		// Normalize the orientation vector.
		NormVec(drSum);
		// Measure the end-to-end distance of the filament.
		V3SUB(dr, &P2(rAct,curr,0), &P2(rAct,n,0));
		VSV3ADD(dr, dr, sft, dimDom);
		lenEE = V3LEN(dr);
		if (lenSeg > 0.) {
			tensFactSum /= lenSeg;
			bendFactSum /= lenSeg;
			V3SCALE(viscFsum, INV(lenSeg));
		}
		if (cntXlink[0] + cntXlink[1] > 0) {
			tensFabpSum /= (double)(cntXlink[0] + cntXlink[1]);
			bendFabpSum /= (double)(cntXlink[0] + cntXlink[1]);
		}
		fprintf(fOut, "%d\t%d\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g"
				"\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%d\t%d\t%d\t%d\n", cnt, 
				iFila[n], n, curr, drSum[0], drSum[1], drSum[2], tensFactSum, 
				tensFabpSum, bendFactSum, bendFabpSum, L_S2UM(lenSeg), 
				L_S2UM(lenCtr), L_S2UM(lenEE), angSum[0], angSum[1], angSum[2],
				viscFsum[0], viscFsum[1], viscFsum[2], cntXlink[0], cntXlink[1],
				cntXlink[2], iPercSum, 0, iSuppSum);
		cnt++;
	}
	fclose(fOut);

	// Record information in unit of filaments
	fOut = fopen(GenFileName("InfoFilaSeg"), "a");
	// Count and print the number of actin filaments
	cnt = 0;
	FOR_ACT(n) {
		CONT(!(P2A(chAct,n,0,nChAc) > -1 && P2A(chAct,n,1,nChAc) < 0));
		CS = 0;
		curr = n;
		while(P2A(chAct,curr,0,nChAc) > -1) {
			for(k = 0; k < nChAcX; k++) {
				for(l = 0; l < nChAcY; l++) {
					abpInd = P2A(chAct,curr,2 + k * nChAcY + l,nChAc);
					CONT(abpInd < 0);
					CONT(!(ISMTF(P2A(chAbp,abpInd,2,nChAb)))
							&& P2A(chAbp,abpInd,1,nChAb) < 0);
					if (CS == 1) { 
						cnt++;
					}
					CS = 1;
					prevAbpInd = abpInd;
				}
			}
			curr = P2A(chAct,curr,0,nChAc);
		}
	}
    Fprintf1dArrayIntWFil(fOut, &cnt, 1, 0, 22);
	// Record the information
	cnt = 0;
	FOR_ACT(n) {
		CONT(!(P2A(chAct,n,0,nChAc) > -1 && P2A(chAct,n,1,nChAc) < 0));
		CS = 0;
		curr = n;
		while(P2A(chAct,curr,0,nChAc) > -1) {
			pArr = &P2A(chAct,curr,0,nChAc);
			len = CalcVecDist(dr, &P2(rAct,pArr[0],0), &P2(rAct,curr,0), 0);
			if (CS == 1) {
				VVS3ADD(drSum, dr, INV(len));
				tensFactSum += tensF[curr];
				bendFactSum += bendF[curr];
				VV3ADD(viscFsum, &P2(viscF,curr,0));
				cntFor++;
				iPercSum += 0;
				iSuppSum += 0;
			}
			for(k = 0; k < nChAcX; k++) {
				for(l = 0; l < nChAcY; l++) {
					abpInd = pArr[2 + k * nChAcY + l];
					CONT(abpInd < 0);
					pArr2 = &P2A(chAbp,abpInd,0,nChAb);
					CONT(!(ISMTF(pArr2[2])) && pArr2[1] < 0);
					CalcPosOnActSeg(&P2(rAct,curr,0), &P2(rAct,pArr[0],0), 
							rPos[1], (double)k / (double)nChAcX);
					crsAng[1] = -1.;
					if (!(ISMTF(pArr2[2]))) {
						actInd = pArr2[(pArr2[0] == curr) ? 1 : 0];
						CalcVec(dr2, &P2(rAct,P2A(chAct,actInd,0,nChAc),0), 
								&P2(rAct,actInd,0));
						crsAng[1] = V3ANG(dr, dr2);
					}
					if (CS == 1) { 
						if (curr != prevActInd) {
							CheckCrossBound(sft, &P2(rAct,curr,0), rPos[1]);
						}
						V3SUB(dr2, rPos[1], rPos[0]);
						VSV3ADD(dr2, dr2, sft, dimDom);
						lenEE = V3LEN(dr2);
						if (cntFor > 0) {
							tensFactSum /= (double)cntFor;
							bendFactSum /= (double)cntFor;
							V3SCALE(viscFsum, INV((double)cntFor));
						}
						NormVec(drSum);
						fprintf(fOut, "%d\t%d\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t"
								"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%g\t"
								"%g\n", cnt, iFila[n], prevActInd, curr, 
								drSum[0], drSum[1], drSum[2], tensFactSum, 
								bendFactSum, L_S2UM(lenSeg), L_S2UM(lenCtr), 
								L_S2UM(lenEE), angSum[0], angSum[1], angSum[2],
								viscFsum[0], viscFsum[1], viscFsum[2],
								prevAbpInd * 10 + P2A(chAbp,prevAbpInd,2,nChAb),
								abpInd * 10 + pArr2[2], crsAng[0], crsAng[1]);
						cnt++;
					}
					lenCtr = 0.;
					lenSeg = 0.;
					tensFactSum = tensF[curr];
					bendFactSum = bendF[curr];
					V3COPY(viscFsum, &P2(viscF,curr,0));
					cntFor = 1;
					V3SET_ALL(angSum, 0.);
					VS3COPY(drSum, dr, INV(len));
					iPercSum = 0;
					iSuppSum = 0;
					V3SET_ALL(sft, 0);
					prevAbpInd = abpInd;
					prevActInd = curr;
					V3COPY(rPos[0], rPos[1]);
					crsAng[0] = crsAng[1];
					CS = 1;
				}
				if (CS == 1) {
					lenCtr += len / (double)nChAcX;
					lenSeg += 1. / (double)nChAcX;
				}
			}
			if (CS == 1) { 
				CheckCrossBound(sft, &P2(rAct,curr,0), &P2(rAct,pArr[0],0));
				if (pArr[1] > -1) {
					// Calculate bending angles
					CalcVec(dr2, &P2(rAct,curr,0), &P2(rAct,pArr[1],0));
					V3CROSS(drCrs, dr, dr2);
					NormVec(drCrs);
					ang = V3ANG(dr, dr2);
					VVS3ADD(angSum, drCrs, ang);
				}
			}
			curr = pArr[0];
		}
	}
	fclose(fOut);

  }

  if (rank == 0) {
	free(rAct);
	free(rAbp);
	free(fAct);
	free(fAbp);
	free(fBendAct);
	free(fBendAbp);
	free(lAct);
	free(lAbp);
	free(chAct);
	free(chAbp);
	free(cntAct);
	free(cntAbp);
	free(tensF);
	free(bendF);
	free(viscF);
	free(fixAct);
	free(iFila);
	free(binActAll);
	free(binAbpAll);
	if (motSA.gTgl != 0) {
		free(iMot);
	}
	if (rho.tgl != 0) {
		free(actRho);
		free(abpRho);
	}
  }
  free(nActAll);
  free(nAbpAll);
  free(binAct);
  free(binAbp);
}

/*------------------ Information for filaments and segments  -----------------*/

