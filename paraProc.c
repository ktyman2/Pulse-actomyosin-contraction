// ##################################################
// #   update.c - finally revised on Apr 2022       #
// #   coded by Taeyoon Kim                         #
// #   Copyright (C) 2005 - 2022, Taeyoon Kim,      #
// #   All rights reserved.                         #
// ##################################################
// This file contains functions which update things.

/*----------------------- Related to handling long chains --------------------*/

// Subroutine for UpdateLongChainNormalSubroutine()
// This finds the ranks where a particle is possibly located at a next time 
// step based on its location. The ranks found here are stored in pRank.
void UpdateLongChainNormalSubSubroutine1(double *rMol, ListInt *pRank, 
		double limit) {
  int k, n, oftInd[][2] = {{-1, 0}, {0}, {0, 1}}, cntOftInd[] = {2, 1, 2};
  int oft[NDIM], oft2[NDIM], offset[NDIM], idx[NDIM];
  double r[NDIM], bound;

  pRank->c = 0;
  V3COPY(r, rMol);
  V3SET_ALL(offset, 1);
  V3SET_ALL(oft, 1);

  // oft[] indicates a current rank, and oft2 represents possible ranks.
  FOR_NDIM(k) {
	if (iCell[k] > 0 && iCell[k] < nCell[k] - 1) {
		if (r[k] < P2A(bnd.r,0,k,NDIM)) { oft[k] = 0; }
		else if (r[k] > P2A(bnd.r,1,k,NDIM)) { oft[k] = 2; }
		for(n = 0; n < 2; n++) {
			CONT(!(fabs(r[k] - P2A(bnd.r,n,k,NDIM)) < limit));
			if (r[k] < P2A(bnd.r,n,k,NDIM)) { offset[k] = 2; }
			else { offset[k] = 0; }
		}
	}
	else if ((iCell[k] == 0 || iCell[k] == nCell[k] - 1) && nCell[k] > 1) {
		if (iCell[k] == 0) {
			if (r[k] >= rGrid[k][1] && r[k] < rGrid[k][2]) { oft[k] = 2; }
			else if (r[k] >= rGrid[k][2]) { oft[k] = 0; }
		}
		else {
			if (r[k] >= rGrid[k][nGrid[k] - 3] && r[k] < rGrid[k][nGrid[k] - 2])
			{ oft[k] = 0; }
			else if (r[k] < rGrid[k][nGrid[k] - 3]) { oft[k] = 2; }
		}
        for(n = 0; n < 2; n++) {
			if (n == 0) {
				bound = (iCell[k] == 0) ? rGrid[k][nGrid[k] - 1]
						: rGrid[k][nGrid[k] - 2];
			}
			else {
				bound = (iCell[k] == 0) ? rGrid[k][1] : rGrid[k][0];
			}
			if (r[k] > rGrid[k][nGrid[k] - 1] - limit) { offset[k] = 2; }
			else if (r[k] < rGrid[k][0] + limit) { offset[k] = 0; }
			else if (fabs(r[k] - bound) < limit) {
				if (r[k] < bound) { offset[k] = 2; }
				else { offset[k] = 0; }
			}
		}
	}
  }
  for(idx[0] = 0; idx[0] < cntOftInd[offset[0]]; idx[0]++) {
	for(idx[1] = 0; idx[1] < cntOftInd[offset[1]]; idx[1]++) {
		for(idx[2] = 0; idx[2] < cntOftInd[offset[2]]; idx[2]++) {
			FOR_NDIM(k) { 
				oft2[k] = oft[k] + oftInd[offset[k]][idx[k]]; 
				if (oft2[k] == 3) { oft2[k] = 1; }
				else if (oft2[k] == -1) { oft2[k] = 1; }
			}
			V3IND_BACK_CONST_INT(pRank->l[pRank->c], oft2, NDIM);
			(pRank->c)++;
		}
	}
  }
}

// Subroutine for UpdateLongChainNormalSubroutine()
void UpdateLongChainNormalSubSubroutine2(int rankMol, int rankL, int ind1, 
		int ind2, int ind3) {
  int CS;

  if (rankMol != rank) {
	// "longChIntMsg.l" has the list of actins which a 
	// current subdomain has to transfer to other subdomains 
	// where ABP can belong at next step
	if (rankL == 13) {
		CS = Find2ElementArray(longChIntMsg.l,
				longChIntMsg.c, ind1, ind2, 0, 2);
		if (CS < 0) {
			V2SET(&P2A(longChIntMsg.l,longChIntMsg.c,0,2), ind1, ind2);
			(longChIntMsg.c)++;
		}
	}
	// "longChExtMsg.l" has the list of actins which the 
	// other subdomains have to transfer to subdomains where 
	// ABP can belong at next step
	else {
		if (ind2 != ind3) {
			CS = Find3ElementArray(longChExtMsg.l,
					longChExtMsg.c, ind1, ind2, ind3, 0, 3);
			if (CS < 0) {
				V3SET(&P2A(longChExtMsg.l,longChExtMsg.c,
						0,3), ind1, ind2, ind3);
				(longChExtMsg.c)++;
			}
		}
	}
  }
}

// Subroutine for UpdateLongChainNormal()
void UpdateLongChainNormalSubroutine(int ind1, int ind2, double len) {
  int n, k, CS, CS2; 
  int abpRankInd, actRankInd, abpRankMol, actRankMol, ind;
  double r[NDIM];
  ListInt abpRank, actRank;

  MALLOC(abpRank.l, int, nCpu);
  MALLOC(actRank.l, int, nCpu);

  CS = 0; 
  for(n = 0; n < 2; n++) { 
	ind = (n == 0) ? ind1 : ind2;
	if (ind < nAct) { V3COPY(r, &P2(act.r,iAct[ind],0)); }
	else if (ind >= nAct && ind < nAct + nAbp) 
	{ V3COPY(r, &P2(abp.r,iAbp[ind - nAct],0)); }
	// either of two elements should be placed near boundary
	FOR_NDIM(k) {
		CONT(r[k] >= P2A(bnd.r,0,k,NDIM) + P2A(edge,0,k,NDIM) 
				&& r[k] < P2A(bnd.r,1,k,NDIM) - P2A(edge,1,k,NDIM));
		CS++;	
		break; 
	}
	BREAK(CS > 0);
  }

  // If the current chain length is long enough, they should be in longCh.l
  if (CS > 0 && len >= neiEdge * 0.9) {
	CS2 = InsertLongChain(ind1, ind2, len * 1.1);
	// If it is in longCh.l already, the chain length is just updated
	if (CS2 > -1) { longChDist[CS2] = len * 1.1; }
  }
  // If chain length is not long enough, or if both of them exist at outside
  // of synchronized regions, they are removed from longCh.l
  else { DeleteLongChain(ind1, ind2); }

  // If this is an ABP-actin chain, another actin located in a barbed direction
  // should be transferred for calculation of bending in actin-ABP-actin 
  // or right angles between ABP arms and axes of actin filaments.
  if (CS > 0 && ind1 >= nAct && ind1 < nAct + nAbp && ind2 < nAct) {
	// Find rank to which ABP or actin belongs 
	abpRankMol = CalcRankMolecule(&P2(abp.r,iAbp[ind1 - nAct],0));
	abpRank.c = 0;
	// Find all possible ranks to which ABP or actin can belong at next step
	UpdateLongChainNormalSubSubroutine1(&P2(abp.r,iAbp[ind1 - nAct],0), 	
			&abpRank, maxDisp);
	ind = P2A(act.ch,iAct[ind2],0,nChAc);
	actRank.c = 0;
	if (iAct[ind] < 0) {
		actRankMol = -1;
		// If the current subdomain doesn't have information of actin,
		// it scans all possible subdomains using other actin. 
		UpdateLongChainNormalSubSubroutine1(&P2(act.r,iAct[ind2],0), 
				&actRank, maxActCh);
	}
	else {
		actRankMol = CalcRankMolecule(&P2(act.r,iAct[ind],0)); 
		// Find all possible ranks to which actin can belong at next step
		UpdateLongChainNormalSubSubroutine1(&P2(act.r,iAct[ind],0),
				&actRank, maxDisp);
	}
	for(n = 0; n < actRank.c; n++) {
		actRankInd = P2A(adjRank,actRank.l[n],0,2);
		// If the subdomain does not exist, skip this procedure
		CONT(actRank.l[n] != 13 && actRankInd < 0);
		for(k = 0; k < abpRank.c; k++) {
			abpRankInd = P2A(adjRank,abpRank.l[k],0,2);
			// If the subdomain does not exist, skip this procedure
			CONT(abpRank.l[k] != 13 && abpRankInd < 0);
			// If ABP is not in the current subdomain
			UpdateLongChainNormalSubSubroutine2(abpRankMol, actRank.l[n],
					ind, abpRankInd, actRankInd);
			// If actin belongs to other subdomains
			UpdateLongChainNormalSubSubroutine2(actRankMol, abpRank.l[k],
					ind1, actRankInd, abpRankInd);
		}
	}
  }
  free(actRank.l);
  free(abpRank.l);
}

// normal method.
void UpdateLongChainNormal(void) {
  int n, k;
  int abpInd, actInd, locActInd, locAbpInd;
  int *pArr, *chkL;
  double len;

  chkL = allIntL;
  memset(chkL, -1, sizeof(int) * nActMe * 2);

  CheckArraySize(&longChExtMsg, &longChExtMsg.siz, 3, 0);
  CheckArraySize(&longChIntMsg, &longChIntMsg.siz, 2, 0);

  // Initialize messgages
  longChExtMsg.c = 0;
  longChIntMsg.c = 0;

  FOR_ABPME(n) {
	CONT(ISABPM(n));
	pArr = &P2A(abp.ch,n,0,nChAb);
	// If motors don't walk, inactive motors cannot have long chains, so
	// they don't need to be considered here
	for (k = 0; k < 2; k++) {
		actInd = pArr[k];
		CONT(actInd < 0);
		locActInd = iAct[actInd];
		CONT(locActInd < 0);
		len = CalcDist(&P2(abp.r,n,0), &P2(act.r,locActInd,0), 0);
		UpdateLongChainNormalSubroutine(abp.id[n] + nAct, actInd, len);
	}
	if (ISMTF(K_ABP(n)) && motSA.cenDist * 1.5 > neiEdge) {
		if (abp.mId[n] == 0) {
			for(k = 0; k < 2; k++) {
				abpInd = pArr[k + 3];
				CONT(abpInd < 0);
				locAbpInd = iAbp[abpInd];
				CONT(locAbpInd < 0);
				CONT(abp.mId[locAbpInd] != 0);
				len = CalcDist(&P2(abp.r,locAbpInd,0), &P2(abp.r,n,0), 0);
				UpdateLongChainNormalSubroutine(abp.id[n] + nAct, 
						abpInd + nAct, len);
			}	
		}
	}
  }

  FOR_ACTME(n) {
	CONT(ISACTM(n));
	// Actins usually have very high extensional stiffness. However, If a chain
	// between actins can be extended a lot due to motors or external forces. 
	// The following part takes care of such long actin chains.
	if (maxActCh > neiEdge * 1.25) {
		for(k = 0; k < 2; k++) {
			CONT(P2A(chkL,n,k,2) == 1);
			actInd = P2A(act.ch,n,k,nChAc);
			CONT(actInd < 0);
			locActInd = iAct[actInd];
			CONT(locActInd < 0);
			len = CalcDist(&P2(act.r,locActInd,0), &P2(act.r,n,0), 0);
			UpdateLongChainNormalSubroutine(actInd, act.id[n], len);
			if (locActInd < nActMe) {
				P2A(chkL,locActInd,1 - k,2) = 1;
			}
		}
	}

	// This is for a case with actin in a current subdomain and 
	// ABP in other subdomain
	for(k = 0; k < nChAc - 2; k++) {
		abpInd = P2A(act.ch,n,k + 2,nChAc);
		CONT(abpInd < 0);
		locAbpInd = iAbp[abpInd];
		CONT(locAbpInd < nAbpMe);
		len = CalcDist(&P2(abp.r,locAbpInd,0), &P2(act.r,n,0), 0);
		UpdateLongChainNormalSubroutine(abpInd + nAct, act.id[n], len);
	}
  }
}
	
// Delete an element in "longCh.l"
void DeleteLongChain(int ind1, int ind2) {
  int n, k, *pArr;

  for(n = 0; n < longCh.c; n++) {
	pArr = &P2A(longCh.l,n,0,2);
	CONT(!((pArr[0] == ind1 && pArr[1] == ind2) 
			|| (pArr[0] == ind2 && pArr[1] == ind1)));
	DeleteElementArrayByIndex(longCh.l, &longCh.c, n, 2);
	for(k = n; k < longCh.c; k++) {
		longChDist[k] = longChDist[k + 1];
	}
	break;			
  }
}

// Insert an element in "longCh.l"
int InsertLongChain(int ind1, int ind2, double len) {
  int CS, ind[2];

  if (ind1 > ind2) { V2SET(ind, ind1, ind2); }
  else { V2SET(ind, ind2, ind1);  }
  CS = Find2ElementArray(longCh.l, longCh.c, ind[0], ind[1], 0, 2);
  if (CS < 0) {
	V2SET(&P2A(longCh.l,longCh.c,0,2), ind[0], ind[1]);
	longChDist[longCh.c] = len;
	(longCh.c)++;
  }
  return CS;
}

/*----------------------- Related to handling long chains --------------------*/

/*--------------------- Move particles between subdomains --------------------*/

// Subroutine for MovePartilces()
void MoveParticlesSubroutine1(ListInt *longChMv, double *longChDistMv, 
		int ind) {
  int n, k, oppInd; 

  for(n = longCh.c - 1; n >= 0; n--) {
	for(k = 0; k < 2; k++) {
		BREAK(P2A(longCh.l,n,k,2) == ind);
	}
	CONT(k == 2);
	// copy information from longCh.l to longChMv.l
	V2COPY(&P2A(longChMv->l,longChMv->c,0,2), &P2A(longCh.l,n,0,2));
	longChDistMv[longChMv->c] = longChDist[n];
	(longChMv->c)++;
	// delete them in longCh.l
	oppInd = P2A(longCh.l,n,1 - k,2);
	if (oppInd < nAct) { CONT(iAct[oppInd] > -1); }
	else if (oppInd >= nAct && oppInd < nAct + nAbp) 
	{ CONT(iAbp[oppInd - nAct] > -1); }
	DeleteElementArrayByIndex(longCh.l, &longCh.c, n, 2);
	for(k = n; k < longCh.c; k++) {
		longChDist[k] = longChDist[k + 1];
	}
  }
}

// mode = 0: actin, 1: ABP
void MoveParticlesSubroutine2(int locInd, int mode) {
  int m, n;
  
  if (mode == 0) {
	for(n = locInd; n < nActMe - 1; n++) {
	    act.id[n] = act.id[n + 1];
	    V3COPY(&P2(act.r,n,0), &P2(act.r,n + 1,0));
	    V3COPY(&P2(act.f,n,0), &P2(act.f,n + 1,0));
	    Copy1dArrayInt(&P2A(act.ch,n,0,nChAc), 
				&P2A(act.ch,n + 1,0,nChAc), nChAc);
	    V3COPY(&P2(act.rPrev,n,0), &P2(act.rPrev,n + 1,0));
	    act.fix[n] = act.fix[n + 1];
		if (rho.tgl != 0) { act.rho[n] = act.rho[n + 1]; }
	    act.iF[n] = act.iF[n + 1];
  		if (confVmdInfo.tgl != 0) {
	        recAct.len[n] = recAct.len[n + 1];
	        recAct.sprF[n] = recAct.sprF[n + 1];
	        recAct.bendF[n] = recAct.bendF[n + 1];
	        V4COPY(&P2A(recAct.allF,n,0,NDIM + 1), 
					&P2A(recAct.allF,n + 1,0,NDIM + 1));
	        recAct.cnt[n] = recAct.cnt[n + 1];
	    }
	    iAct[act.id[n]] = n;
	}
	nActMe--;
  }
  else if (mode == 1) {
	for(m = locInd; m < nAbpMe - 1; m++) {
		abp.id[m] = abp.id[m + 1];
		V3COPY(&P2(abp.r,m,0), &P2(abp.r,m + 1,0));
		V3COPY(&P2(abp.f,m,0), &P2(abp.f,m + 1,0));
		V4COPY(&P2A(recInstSprFabp,m,0,4), &P2A(recInstSprFabp,m + 1,0,4));
		Copy1dArrayInt(&P2A(abp.ch,m,0,nChAb), 
				&P2A(abp.ch,m + 1,0,nChAb), nChAb);
		if (rho.tgl != 0) { 
			abp.rho[m] = abp.rho[m + 1]; 
			abp.kind[m] = abp.kind[m + 1]; 
		}
		if (motSA.gTgl != 0) {
			abp.mId[m] = abp.mId[m + 1];
		}
  		if (confVmdInfo.tgl != 0) {
			Copy1dArrayDouble(&P2A(recAbp.len,m,0,recAbp.nL), 
					&P2A(recAbp.len,m + 1,0,recAbp.nL), recAbp.nL);
			recAbp.sprF[m] = recAbp.sprF[m + 1];
			recAbp.bendF[m] = recAbp.bendF[m + 1];
	        V4COPY(&P2A(recAbp.allF,m,0,NDIM + 1), 
					&P2A(recAbp.allF,m + 1,0,NDIM + 1));
			recAbp.cnt[m] = recAbp.cnt[m + 1];
		}
		if (recLongF.tgl != 0) {
			V4COPY(&P2A(recLongSprFabp,m,0,4), &P2A(recLongSprFabp,m + 1,0,4));
		}
		if (tglRecAbpTurn != 0) {
			Copy1dArrayDouble(&P2A(abpTurn,m,0,7), &P2A(abpTurn,m + 1,0,7), 7);
		}
		iAbp[abp.id[m]] = m;
	}
	nAbpMe--;
  }
}
// Transfer particles which move out of the current subdomain to adjacent
// subdomains.
void MoveParticles(void) {
  int m, n, k, l, cnt, cntRecvMsg, cntElem, CS, *pArr, kind, tempInt;
  int ind[3], ind2, actInd, abpInd, indAdjCellMv; 
  int *chkL, oftMv[NDIM];
  int nActMvSub, nAbpMvSub, nChAbMv;
  int begin, end, rep, posi, tag = 0;
  double r[NDIM], *longChDistMv;
  ListInt longChMv, actCh;

  nChAbMv = (motSA.gTgl != 0) ? nChAb : nChAb - 2;
  MALLOC(chkL,int,cntAdjRank[0]);
  longChMv.l = allIntL;
  longChDistMv = allDblL;
  if (modeActCh != 0) {
	MALLOC(actCh.l,int,2 * nChAc - 1);
  }

  FOR_ACTCP(n) { iAct[act.id[n + nActMe]] = -1; }
  FOR_ABPCP(n) { iAbp[abp.id[n + nAbpMe]] = -1; }

  nActCp = 0;
  nAbpCp = 0;

  insNeiPar.c = 0; 

  rep = 1;
  for(m = 0; m < rep; m++) {
	begin = 0;   
	end = NDIM; 
	if (cntAdjRank[m] > 0) {
		// Move particles
		memset(cntMvPar, 0, sizeof(int) * cntAdjRank[m] * 3);
		for(n = 0; n < cntAdjRank[m]; n++) { 
			mvPar[n].c = 0; 
		}
		// Find particles to be transferred
		for(n = 0; n < nActMe + nAbpMe; n++) {
			kind = SetKind(n, nActMe, nAbpMe);
			if (kind == 0) { 
				CONT(ISACTM(n));
				V3COPY(r, &P2(act.r,n,0)); 
			}
			else if (kind == 1) { 
				CONT(ISABPIM(n - nActMe));
				V3COPY(r, &P2(abp.r,n - nActMe,0)); 
			}
			CS = 1;
			V3SET_ALL(oftMv, 1);

			for(k = begin; k < end; k++) {
				// If there is no periodic boundary condition, particles 
				// crossing the outer boundary of the domain are ignored.
				if (pbc[k] == 0 && ((iCell[k] == 0 && r[k] < rGrid[k][0])
						|| (iCell[k] == nCell[k] - 1 
						&& r[k] >= rGrid[k][nGrid[k] - 1]))) { 
					CS = 0;   
					continue;   
				}
				// Treat particles crossing the outer boundary.
				else if (nCell[k] > 1) {
					// For middle cells
					if (iCell[k] > 0 && iCell[k] < nCell[k] - 1) {
						if (r[k] < P2A(bnd.r,0,k,NDIM)) { oftMv[k] -= 1; }
						else if (r[k] >= P2A(bnd.r,1,k,NDIM))
						{ oftMv[k] += 1; }
					}
					// For cells located at edges
					else if (iCell[k] == 0 && r[k] >= P2A(bnd.r,1,k,NDIM)) {
						if (r[k] >= rGrid[k][2]) { oftMv[k] -= 1; }
						else { oftMv[k] += 1; }
					}
					else if (iCell[k] == nCell[k] - 1 
							&& r[k] < P2A(bnd.r,0,k,NDIM)) {
						if (r[k] >= rGrid[k][nGrid[k] - 3]) { oftMv[k] -= 1; }
						else { oftMv[k] += 1; }
					}
				}
			}
			CS = 0;
			// Normal method
			V3IND_BACK_CONST_INT(indAdjCellMv, oftMv, NDIM);
		    ind2 = P2A(adjRank,indAdjCellMv,1,2);
			// "indAdjCellMv == 13" corrsponds to the current cell. (1,1,1)
			if (indAdjCellMv != 13 && ind2 > -1) { CS = 1; }
			// Add particles in list
			if (CS == 1) {
	            InsertElement1dArrayWoChk(mvPar[ind2].l, &mvPar[ind2].c, n);
				P2A(cntMvPar,ind2,kind,3)++;
			}
		}
	}
	// Send chosen particles
	// The information of the particles is packed into a message
	for(n = 0; n < cntAdjRank[m]; n++) {
		posi = 0;
		V3COPY(ind, &P2A(cntMvPar,n,0,3));
		longChMv.c = 0;

		MPI_PACK_INT(&P2A(cntMvPar,n,0,3), 3, n);
		// actin
		for(k = 0; k < ind[0]; k++) {
			actInd = mvPar[n].l[k];
			MPI_PACK_INT(&act.id[actInd], 1, n);
			MPI_PACK_DBL(&P2(act.r,actInd,0), NDIM, n);
			MPI_PACK_DBL(&P2(act.f,actInd,0), NDIM, n);
			if (modeActCh == 0) { 
				MPI_PACK_INT(&P2A(act.ch,actInd,0,nChAc), nChAc, n);
			}
			else {		
				V2COPY(&actCh.l[1], &P2A(act.ch,actInd,0,nChAc));
				actCh.c = 2;
				for(l = 2; l < nChAc; l++) {
					if (P2A(act.ch,actInd,l,nChAc) > -1) {
						V2SET(&actCh.l[actCh.c + 1], l, 
								P2A(act.ch,actInd,l,nChAc));
						actCh.c += 2;
					}
				}
				actCh.l[0] = actCh.c;
				MPI_PACK_INT(actCh.l, actCh.c + 1, n);
			}
			MPI_PACK_DBL(&P2(act.rPrev,actInd,0), NDIM, n);
			MPI_PACK_INT(&act.fix[actInd], 1, n);
			if (rho.tgl != 0) { MPI_PACK_INT(&act.rho[actInd], 1, n); }

			MPI_PACK_INT(&act.iF[actInd], 1, n);
  			if (confVmdInfo.tgl != 0) {
				MPI_PACK_DBL(&recAct.len[actInd], 1, n);
				MPI_PACK_DBL(&recAct.sprF[actInd], 1, n);
				MPI_PACK_DBL(&recAct.bendF[actInd], 1, n);
				MPI_PACK_DBL(&P2A(recAct.allF,actInd,0,NDIM + 1), NDIM + 1, n);
				MPI_PACK_INT(&recAct.cnt[actInd], 1, n);
			}
			if (tglActMoDyn != 0) {
				CS = FindElementArray(noActDyn.l, noActDyn.c, 
						act.id[actInd], 0, 2);
				MPI_PACK_INT(&CS, 1, n);
				if (CS > -1) {
					MPI_PACK_INT(&P2A(noActDyn.l,CS,0,2), 2, n);
					DeleteElementArrayByIndex(noActDyn.l, &noActDyn.c, CS, 2);
				}
			}
			if (tglAbpAcInaDyn != 0 || acpMoBind.tgl != 0 || motMoBind.tgl != 0 
					|| (actDis.tgl != 0 && actDis.facKWA > 0.)) { 
				CS = FindElementArray(noAbpDyn.l, noAbpDyn.c, 
						act.id[actInd], 1, 3);
				MPI_PACK_INT(&CS, 1, n);
				if (CS > -1) {
					MPI_PACK_INT(&P2A(noAbpDyn.l,CS,0,3), 3, n);
		        }
			}
			MoveParticlesSubroutine1(&longChMv, longChDistMv, 
					act.id[actInd]);
			iAct[act.id[actInd]] = -1;
			act.id[actInd] = -1;
		}
		// ABP
		for(k = 0; k < ind[1]; k++) {
			abpInd = mvPar[n].l[k + ind[0]] - nActMe;
			pArr = &P2A(abp.ch,abpInd,0,nChAb);
			MPI_PACK_INT(&abp.id[abpInd], 1, n);
			MPI_PACK_DBL(&P2(abp.r,abpInd,0), NDIM, n);
			MPI_PACK_DBL(&P2(abp.f,abpInd,0), NDIM, n);
			MPI_PACK_DBL(&P2A(recInstSprFabp,abpInd,0,4), 4, n);
			MPI_PACK_INT(pArr, nChAbMv, n);
			if (rho.tgl != 0) { 
				MPI_PACK_INT(&abp.rho[abpInd], 1, n); 
				MPI_PACK_INT(&abp.kind[abpInd], 1, n); 
			}
			if (motSA.gTgl != 0) {
				MPI_PACK_INT(&abp.mId[abpInd], 1, n);
			}
  			if (confVmdInfo.tgl != 0) {
				MPI_PACK_DBL(&P2A(recAbp.len,abpInd,0,recAbp.nL), recAbp.nL, n);
				MPI_PACK_DBL(&recAbp.sprF[abpInd], 1, n);
				MPI_PACK_DBL(&recAbp.bendF[abpInd], 1, n);
				MPI_PACK_DBL(&P2A(recAbp.allF,abpInd,0,NDIM + 1), NDIM + 1, n);
				MPI_PACK_INT(&recAbp.cnt[abpInd], 1, n);
			}
			if (recLongF.tgl != 0) {
				MPI_PACK_DBL(&P2A(recLongSprFabp,abpInd,0,4), 4, n);
			}
			MoveParticlesSubroutine1(&longChMv, longChDistMv,
					abp.id[abpInd] + nAct);
			if (pArr[0] > -1 && pArr[1] > -1) {
				if (tglNeiAbpDC != 0) {
					DeleteElementInNeighborList(abp.id[abpInd], 1);
				}
			}
			else {
				if (tglNeiAbpSC != 0) {
					DeleteElementInNeighborList(abp.id[abpInd], 1);
				}
				if (K_ABP(abpInd) == 2) { 
					if (pArr[0] > -1) { nMotInaMe--; }
					else { nMotMme--; }
				}
				else {
					if (pArr[0] > -1) { nAcpInaMe--; }
					else { nAcpMme--; }
				}
			}
			if ((tglAcpAcInaDyn != 0 && K_ABP(abpInd) != 2) 
					|| (tglMotAcInaDyn != 0 && K_ABP(abpInd) == 2)
					|| (actDis.tgl != 0 && actDis.facKWA > 0.)) {
				CS = FindElementArray(noAbpDyn.l, noAbpDyn.c, 
						abp.id[abpInd], 0, 3);
				MPI_PACK_INT(&CS, 1, n);
				if (CS > -1) {
					MPI_PACK_INT(&P2A(noAbpDyn.l,CS,0,3), 3, n);
		        }
			}
			if (tglRecAbpTurn != 0) {
				MPI_PACK_DBL(&P2A(abpTurn,abpInd,0,7), 7, n);
			}
			iAbp[abp.id[abpInd]] = -1;
			abp.id[abpInd] = -1;
		}
		// Additional information to be packed
		MPI_PACK_INT(&longChMv.c, 1, n);
		MPI_PACK_INT(longChMv.l, longChMv.c * 2, n);
		MPI_PACK_DBL(longChDistMv, longChMv.c, n);
		// If it's the subdomain where actin in longChExtMsg.l can exist
		// at next step, the message is transferred. 
		longChMv.c = 0;
		for(k = 0; k < longChExtMsg.c; k++) {
			if (iRank[n] == P2A(longChExtMsg.l,k,2,3)) { 
				V2COPY(&P2A(longChMv.l,longChMv.c,0,2),
						&P2A(longChExtMsg.l,k,0,3));
				(longChMv.c)++;
			}
		}
		MPI_PACK_INT(&longChMv.c, 1, n);
		MPI_PACK_INT(longChMv.l, longChMv.c * 2, n);
		// Send and receive the messages to adjacent subdomains
		MPI_Isend(bufSendMsg[n], posi, MPI_PACKED, iRank[n], 
				tag, MPI_COMM_WORLD, &sReq[n]);
		MPI_Irecv(bufRecvMsg[n], sizeBufMsg, MPI_PACKED, 
				iRank[n], tag, MPI_COMM_WORLD, &rReq[n]);
	}
	if (cntAdjRank[m] > 0) {
		// Delete moved particles
		// actin
		for(n = nActMe - 1; n >= 0; n--) {  
			CONT(!(act.id[n] < 0));
			MoveParticlesSubroutine2(n, 0);
		}
		// ABP
		for(n = nAbpMe - 1; n >= 0; n--) {
			CONT(!(abp.id[n] < 0));
			MoveParticlesSubroutine2(n, 1);
		}

		// Receive information of particles transferred from adjacent domains
		cnt = 0;
		cntRecvMsg = 0;
		memset(mpiTestRecvFlag, 0, sizeof(int) * cntAdjRank[m]);
		memset(mpiTestSendFlag, 0, sizeof(int) * cntAdjRank[m]);
		memset(chkL, 0, sizeof(int) * cntAdjRank[m]);
		while(cntRecvMsg < cntAdjRank[m]) {
			if (mpiTestSendFlag[cnt] == 0) {
				MPI_Test(&sReq[cnt], &mpiTestSendFlag[cnt], &status);
			}
			if (mpiTestRecvFlag[cnt] == 0 && mpiTestSendFlag[cnt] != 0 
					&& chkL[cnt] == 0) {
				MPI_Test(&rReq[cnt], &mpiTestRecvFlag[cnt], &status);
				if (mpiTestRecvFlag[cnt] != 0) { 
					posi = 0;
			  		MPI_UNPACK_INT(&nActMvSub, 1, cnt);
			  		MPI_UNPACK_INT(&nAbpMvSub, 1, cnt);
			  		MPI_UNPACK_INT(&tempInt, 1, cnt);
					// actin
					for(n = 0; n < nActMvSub; n++) {
						actInd = nActMe + n;
				  		pArr = &P2A(act.ch,actInd,0,nChAc);
				  		MPI_UNPACK_INT(&act.id[actInd], 1, cnt);
				  		MPI_UNPACK_DBL(&P2(act.r,actInd,0), NDIM, cnt);
				  		MPI_UNPACK_DBL(&P2(act.f,actInd,0), NDIM, cnt);
						if (modeActCh == 0) {
				  			MPI_UNPACK_INT(pArr, nChAc, cnt);
						}
						else {
							memset(&pArr[2], -1, sizeof(int) * (nChAc - 2));
				  			MPI_UNPACK_INT(&actCh.c, 1, cnt);
				  			MPI_UNPACK_INT(actCh.l, actCh.c, cnt);
							V2COPY(pArr, actCh.l);
							for(k = 2; k < actCh.c; k += 2) {
								pArr[actCh.l[k]] =  actCh.l[k + 1];
							}
						}
						iAct[act.id[actInd]] = actInd;
				  		MPI_UNPACK_DBL(&P2(act.rPrev,actInd,0), NDIM, cnt);
					  	MPI_UNPACK_INT(&act.fix[actInd], 1, cnt);
						if (rho.tgl != 0) {
						  	MPI_UNPACK_INT(&act.rho[actInd], 1, cnt);
						}
					  	MPI_UNPACK_INT(&act.iF[actInd], 1, cnt);
			  			if (confVmdInfo.tgl != 0) {
					  		MPI_UNPACK_DBL(&recAct.len[actInd], 1, cnt);
					  		MPI_UNPACK_DBL(&recAct.sprF[actInd], 1, cnt);
					  		MPI_UNPACK_DBL(&recAct.bendF[actInd], 1, cnt);
					  		MPI_UNPACK_DBL(&P2A(recAct.allF,actInd,0,NDIM + 1), 
									NDIM + 1, cnt);
					  		MPI_UNPACK_INT(&recAct.cnt[actInd], 1, cnt);
						}
						if (tglActMoDyn != 0) {
					  		MPI_UNPACK_INT(&CS, 1, cnt);
							if (CS > -1) {
					  			MPI_UNPACK_INT(&P2A(noActDyn.l, 
										noActDyn.c,0,2), 2, cnt);
								(noActDyn.c)++;
							}
						}
						if (tglAbpAcInaDyn != 0 || acpMoBind.tgl != 0 
								|| motMoBind.tgl != 0 
								|| (actDis.tgl != 0 && actDis.facKWA > 0.)) { 
					  		MPI_UNPACK_INT(&CS, 1, cnt);
							if (CS > -1) {
					  			MPI_UNPACK_INT(&P2A(noAbpDyn.l, 
									noAbpDyn.c,0,3), 3, cnt);
								(noAbpDyn.c)++;
							}
						}
					}
					// ABP
					for(n = 0; n < nAbpMvSub; n++) {
						abpInd = nAbpMe + n;
						pArr = &P2A(abp.ch,abpInd,0,nChAb);
				  		MPI_UNPACK_INT(&abp.id[abpInd], 1, cnt);
				  		MPI_UNPACK_DBL(&P2(abp.r,abpInd,0), NDIM, cnt);
				  		MPI_UNPACK_DBL(&P2(abp.f,abpInd,0), NDIM, cnt);
				  		MPI_UNPACK_DBL(&P2A(recInstSprFabp,abpInd,0,4), 4, cnt);
				  		MPI_UNPACK_INT(pArr, nChAbMv, cnt);
						if (rho.tgl != 0) {
						  	MPI_UNPACK_INT(&abp.rho[abpInd], 1, cnt);
						  	MPI_UNPACK_INT(&abp.kind[abpInd], 1, cnt);
						}
						if (motSA.gTgl != 0) {
				  			MPI_UNPACK_INT(&abp.mId[abpInd], 1, cnt);
						}
						else {
							V2SET_ALL(&pArr[3], -1);
						}
			  			if (confVmdInfo.tgl != 0) {
					  		MPI_UNPACK_DBL(&P2A(recAbp.len,abpInd,0,recAbp.nL), 
									recAbp.nL, cnt);
					  		MPI_UNPACK_DBL(&recAbp.sprF[abpInd], 1, cnt);
					  		MPI_UNPACK_DBL(&recAbp.bendF[abpInd], 1, cnt);
					  		MPI_UNPACK_DBL(&P2A(recAbp.allF,abpInd,0,NDIM + 1), 
									NDIM + 1, cnt);
					  		MPI_UNPACK_INT(&recAbp.cnt[abpInd], 1, cnt);
						}
						if (recLongF.tgl != 0) {
				  			MPI_UNPACK_DBL(&P2A(recLongSprFabp,
									abpInd,0,4), 4, cnt);
						}
						iAbp[abp.id[abpInd]] = abpInd;
		
						if (pArr[0] > -1 && pArr[1] > -1) {
							if (tglNeiAbpDC != 0) {
								InsertElement1dArrayWoChk(insNeiPar.l, 
										&insNeiPar.c, nAct + abp.id[abpInd]);
							}
						}
						else {
							if (tglNeiAbpSC != 0) {
								InsertElement1dArrayWoChk(insNeiPar.l, 
									&insNeiPar.c, nAct + abp.id[abpInd]);
							}
							if (pArr[2] == 2) { 
								if (pArr[0] > -1) { nMotInaMe++; }
								else { nMotMme++; }
							}
							else {
								if (pArr[0] > -1) { nAcpInaMe++; }
								else { nAcpMme++; }
							}
						}
						if ((tglAcpAcInaDyn != 0 && pArr[2] != 2) 
								|| (tglMotAcInaDyn != 0 && pArr[2] == 2)
								|| (actDis.tgl != 0 && actDis.facKWA > 0.)) {
					  		MPI_UNPACK_INT(&CS, 1, cnt);
							if (CS > -1) {
					  			MPI_UNPACK_INT(&P2A(noAbpDyn.l, 
									noAbpDyn.c,0,3), 3, cnt);
								(noAbpDyn.c)++;
							}
						}
						if (tglRecAbpTurn != 0) {
							MPI_UNPACK_DBL(&P2A(abpTurn,abpInd,0,7), 7, cnt);
						}
					}
					// Additional information to be received
					MPI_UNPACK_INT(&cntElem, 1, cnt);
					MPI_UNPACK_INT(longChMv.l, cntElem * 2, cnt);
					MPI_UNPACK_DBL(longChDistMv, cntElem, cnt);
					for(n = 0; n < cntElem; n++) {
						pArr = &P2A(longChMv.l,n,0,2);
						InsertLongChain(pArr[0], pArr[1], longChDistMv[n]);
					}
					// The transferred longChExtMsg is stored in 
					// longChIntMsg.l.
					MPI_UNPACK_INT(&cntElem, 1, cnt);
					MPI_UNPACK_INT(&P2A(longChIntMsg.l,longChIntMsg.c
							,0,2) , cntElem * 2, cnt);
					longChIntMsg.c += cntElem;					
					nActMe += nActMvSub;
					nAbpMe += nAbpMvSub;
					cntRecvMsg++;
					chkL[cnt] = 1;
				}

			}
			cnt++;
			if (cnt == cntAdjRank[m]) { cnt = 0; }
		}
	}
  } 
  free(chkL);
  if (modeActCh != 0) {
	free(actCh.l);
  }
}

/*--------------------- Move particles between subdomains --------------------*/

/*--------------------- Copy particles between subdomains --------------------*/

// Subroutine for CopyPartilcesSubroutine()
void CopyParticlesSubSubroutine(int id, int *idx, int ind) {
  int k, CS;

  // If actin, two actins are transferred
  if (id < nActMe + nActCp) {
	for(k = 0; k < 2; k++) {
		CONT(!(idx[k] > -1));
		CS = InsertElement1dArrayWChk(cpPar[ind].l, &cpPar[ind].c, idx[k]);
		if (CS == -1) { P2A(cntCpPar,ind,0,3)++; }
	}
  }
  // If ABP, only one actin is transferred
  else if (id >= nActMe + nActCp
		&& id < nActMe + nActCp + nAbpMe + nAbpCp) {
	for(k = 0; k < 2; k++) {
		CONT(!(idx[k] > -1));
		CS = FindElementArray(cpPar[ind].l, 
				P2A(cntCpPar,ind,0,3), idx[k], 0, 1);
		CONT(!(CS == -1));
		// Squeeze it in the list
		InsertElementArrayByIndex(cpPar[ind].l,
				&cpPar[ind].c, &idx[k], P2A(cntCpPar,ind,0,3), 1);
		P2A(cntCpPar,ind,0,3)++;
	}
  }
}

// Subroutine for CopyParticles() with normal method
// mode: 0 - normal, 1 - long
void CopyParticlesNormalSubroutine(int id, int *oftCp, int *idx, int mode) { 
  int m, n, k, l, ind, kind, cntCpuList, *cpuList;
  int *cpuListTemp, CS, indAdjCellCp;
  int oftInd[][2] = {{0}, {1}, {2}, {0, 2}}, cntOftInd[4] = {1, 1, 1, 2};
  int offCpuList[][7] = OFFSET_CPU_LIST, offCpuLen[] = OFFSET_CPU_LEN;

  if (mode == 0 && oftCp[0] != 3 && oftCp[1] != 3 && oftCp[2] != 3) {
	indAdjCellCp = ((oftCp[0] * NDIM) + oftCp[1]) * NDIM + oftCp[2];
	if (indAdjCellCp != 13) {
		cntCpuList = offCpuLen[indAdjCellCp];
		cpuList = offCpuList[indAdjCellCp];
	}
	else { return; }
  }
  else {
	MALLOC(cpuList,int,27);
	MALLOC(cpuListTemp,int,27);
	memset(cpuListTemp, -1, sizeof(int) * 27);
	for(m = 0; m < cntOftInd[oftCp[0]]; m++) {
		for(n = 0; n < cntOftInd[oftCp[1]]; n++) {
			for(k = 0; k < cntOftInd[oftCp[2]]; k++) {
	  			indAdjCellCp = ((oftInd[oftCp[0]][m] * NDIM) 
						+ oftInd[oftCp[1]][n]) * NDIM + oftInd[oftCp[2]][k];
				if (indAdjCellCp != 13) {
					for(l = 0; l < offCpuLen[indAdjCellCp]; l++) {
						cpuListTemp[offCpuList[indAdjCellCp][l]] = 1;
					}
				}
			}
		}
	}
	cntCpuList = 0;
	for(n = 0; n < 27; n++) {
		CONT(!(cpuListTemp[n] == 1));
		cpuList[cntCpuList] = n;
		cntCpuList++;
	}
	free(cpuListTemp);
  }
  for(n = 0; n < cntCpuList; n++) {
	ind = P2A(adjRank,cpuList[n],1,2);
	CONT(!(ind > -1));
	CS = InsertElement1dArrayWChk(cpPar[ind].l, &cpPar[ind].c, id);
	if (CS == -1) {
  		kind = SetKind(id, nActMe + nActCp, nAbpMe + nAbpCp);
		P2A(cntCpPar,ind,kind,3)++; 
	}
	CopyParticlesSubSubroutine(id, idx, ind);
  }
  if (!(mode == 0 && oftCp[0] != 3 && oftCp[1] != 3 && oftCp[2] != 3)) {
	free(cpuList); 
  }
}

void CopyParticlesSubroutine(int begin, int end, int *oftCp,
	double *r, double dist) {
  int k;

  V3SET_ALL(oftCp, 1);
  for(k = begin; k < end; k++) {
	if (!(pbc[k] == 0 && iCell[k] == 0) && r[k] >= P2A(bnd.r,0,k,NDIM)
			&& r[k] < P2A(bnd.r,0,k,NDIM) + dist) {
		oftCp[k]--;
	}
	if (!(pbc[k] == 0 && iCell[k] == nCell[k] - 1)
			&& r[k] >= P2A(bnd.r,1,k,NDIM) - dist
			&& r[k] < P2A(bnd.r,1,k,NDIM)) {
		oftCp[k] = (oftCp[k] == 1) ? 2 : 3;
	}
  }
}

// Send information of particles located near boundaries to adjacent boundaries
// for synchronization.
void CopyParticles(void) {
  int m, n, k, l, cnt, cntRecvMsg, CS, oftCp[NDIM], kind;
  int ind[3], ind2, idx[5], *chkL, *pArr, tempInt;
  int actInd, absInd, oppInd;
  int begin, end, rep,  posi, maxPosi, tag = 0, nChAbCp;
  int nActCpSub, nAbpCpSub, nActCpPre, nAbpCpPre;
  double r[NDIM], dist;
  ListInt actCh;

  maxPosi = NEG_LARGE_VALUE; 
  nChAbCp = (motSA.gTgl != 0) ? nChAb : nChAb - 2;
  if (modeActCh != 0) {
	MALLOC(actCh.l,int,2 * nChAc - 1);
  }
  MALLOC(chkL,int,cntAdjRank[0]);
  rep = 1;

  for(m = 0; m < rep; m++) {
	begin = 0;	
	end = NDIM; 
	if (cntAdjRank[m] > 0) {
		memset(cntCpPar, 0, sizeof(int) * cntAdjRank[m] * 3); 
		for(n = 0; n < cntAdjRank[m]; n++) { 
			cpPar[n].c = 0; 
		}
		// Find particles to be synchronized
		for(n = 0; n < nActMe + nActCp + nAbpMe + nAbpCp; n++) {
			kind = SetKind(n, nActMe + nActCp, nAbpMe + nAbpCp);
			if (n < nActMe) { 
				CONT(ISACTM(n));
			}
			else if (n >= nActMe + nActCp && n < nActMe + nActCp + nAbpMe) {
				ind2 = n - nActMe - nActCp;
				CONT(ISABPIM(ind2));
			}
			if (kind == 0) { 
				V3COPY(r, &P2(act.r,n,0)); 
				absInd = act.id[n]; 
			}
			else if (kind == 1) {
				V3COPY(r, &P2(abp.r,n - nActMe - nActCp,0)); 
				absInd = abp.id[n - nActMe - nActCp] + nAct;
			}
			CopyParticlesSubroutine(begin, end, oftCp, r, neiEdge);
			V2SET_ALL(idx, -1);
			// actin
			if (kind == 0) {
				pArr = &P2A(act.ch,n,0,nChAc);
				// If the actin is not the end of filaments,
				if (pArr[0] > -1) {
					// If no ABP is attached on the actin,
					// another actin connected to the actin is transferred.
					// This is for calculation of bending of actin filament.
					if (pArr[1] > -1) {
						for(k = 0; k < 2; k++) {
							CONT(!(iAct[pArr[k]] > -1 
									&& iAct[pArr[1 - k]] < 0));
							idx[k] = iAct[pArr[k]];
							break;
						}
					}
					// If ABP is attached, but if the ABP doesn't exist
					// in the subdomain, both actins connected to 
					// the actin are transferred.
					// This is for calculation of bending of right angle.
					if (idx[0] == -1) {
						for(k = 2; k < nChAc; k++) {
							CONT(pArr[k] < 0);
							CONT(iAbp[pArr[k]] > -1);
							idx[0] = iAct[pArr[0]];
							break;
						}
					}
				}
			}
			// ABP
			else if (kind == 1) {
				ind2 = n - nActMe - nActCp;
				if (abpF.bend[K_ABP(ind2)].facStf > 0. 
						&& !(ISMTF(K_ABP(ind2)))) {
					pArr = &P2A(abp.ch,ind2,0,nChAb);
					// If the ABP is connected to two act.filaments
					if (pArr[0] > -1 && pArr[1] > -1) {
						// If one actin belongs to a current domain,
						// and if the other actin belongs to the other domain,
						// the actin in the current domain should be 
						// transferred. This is for calculating bending between
						// two arms of ABPs.
						for(k = 0; k < 2; k++) {
							CONT(!(iAct[pArr[k]] > -1
									&& iAct[pArr[1 - k]] < 0));
							idx[0] = iAct[pArr[k]];
							idx[1] = iAct[P2A(act.ch,iAct[pArr[k]],0,nChAc)];
							break;
						}
					}
				}
			}
			// Very long chains are handled here..
			CopyParticlesNormalSubroutine(n, oftCp, idx, 0); 
			// Check longCh.l and transfer to more adjacent subdomains
			// if necessary
			for(k = 0; k < longCh.c; k++) {
				pArr = &P2A(longCh.l,k,0,2);
				CS = -1;
				for(l = 0; l < 2; l++) {
					CONT(absInd != pArr[l]);
					oppInd = pArr[1 - l];
					if (oppInd < nAct) {
						if (iAct[oppInd] < 0) { CS = 1; }
					}
					else if (oppInd >= nAct && oppInd < nAct + nAbp) {
						if (iAbp[oppInd - nAct] < 0) { CS = 1; }
					}
					break;
				}
				if (CS == 1) {
					dist = longChDist[k] + 1.5;
					CopyParticlesSubroutine(0, NDIM, oftCp, r, dist);
					CopyParticlesNormalSubroutine(n, oftCp, idx, 1); 
				}
			}
			// If actin, check longChIntMsg.l and transfer actin 
			// if it's in the list.
			for(k = 0; k < longChIntMsg.c; k++) {
				pArr = &P2A(longChIntMsg.l,k,0,2);
				CONT(pArr[0] != absInd);
				// Find the destination
				ind2 = -1;
				for(l = 0; l < 27; l++) {
					CONT(!(pArr[1] == P2A(adjRank,l,0,2)));
					ind2 = P2A(adjRank,l,1,2);
					break;
				}
				CONT(!(ind2 > -1));
				CS = InsertElement1dArrayWChk(cpPar[ind2].l, 
						&cpPar[ind2].c, n);
				CONT(!(CS == -1));
				P2A(cntCpPar,ind2,kind,3)++; 
			}
		}
	}
	// Pack information of the chosen particles into a message
	for(n = 0; n < cntAdjRank[m]; n++) {
		posi = 0;
		V3COPY(ind, &P2A(cntCpPar,n,0,3));
		MPI_PACK_INT(&P2A(cntCpPar,n,0,3), 3, n);
		if (tglActFormDyn != 0) {
			MPI_PACK_INT(&cntNucAss, 1, n);
		}
		// actin
		for(k = 0; k < ind[0]; k++) {
			MPI_PACK_INT(&act.id[cpPar[n].l[k]], 1, n);
		}
		for(k = 0; k < ind[0]; k++) {
			MPI_PACK_DBL(&P2(act.r,cpPar[n].l[k],0), NDIM, n);
		}
		// Send the whole or filled chain information, depending on which
		// is a more efficient way..
		if (modeActCh == 0) { 
			for(k = 0; k < ind[0]; k++) {
				MPI_PACK_INT(&P2A(act.ch,cpPar[n].l[k],0,nChAc), nChAc, n);
			}
		}
		else {		
			for(k = 0; k < ind[0]; k++) {
				actInd = cpPar[n].l[k];
				V2COPY(&actCh.l[1], &P2A(act.ch,actInd,0,nChAc));
				actCh.c = 2;
				for(l = 2; l < nChAc; l++) {
					CONT(!(P2A(act.ch,actInd,l,nChAc) > -1));
					V2SET(&actCh.l[actCh.c + 1], l, 
							P2A(act.ch,actInd,l,nChAc));
					actCh.c += 2;
				}
				actCh.l[0] = actCh.c;
				MPI_PACK_INT(actCh.l, actCh.c + 1, n);
			}
		}
		for(k = 0; k < ind[0]; k++) {
			MPI_PACK_INT(&act.iF[cpPar[n].l[k]], 1, n);
		}
		// ABP
		for(k = 0; k < ind[1]; k++) {
			MPI_PACK_INT(&abp.id[cpPar[n].l[k + ind[0]] 
					- nActMe - nActCp], 1, n);
		}
		for(k = 0; k < ind[1]; k++) {
			MPI_PACK_DBL(&P2(abp.r,cpPar[n].l[k + ind[0]] 
					- nActMe - nActCp,0), NDIM, n);
		}
		for(k = 0; k < ind[1]; k++) {
			MPI_PACK_INT(&P2A(abp.ch,cpPar[n].l[k + ind[0]] 
					- nActMe - nActCp,0,nChAb), nChAbCp, n);
		}
		if (motSA.gTgl != 0) {
			for(k = 0; k < ind[1]; k++) {
				MPI_PACK_INT(&abp.mId[cpPar[n].l[k + ind[0]] 
						- nActMe - nActCp], 1, n);
			}
		}
		// Send and receive the messages to adjacent subdomains
		MPI_Isend(bufSendMsg[n], posi, MPI_PACKED, iRank[n], 
				tag, MPI_COMM_WORLD, &sReq[n]);
		MPI_Irecv(bufRecvMsg[n], sizeBufMsg, MPI_PACKED, 
				iRank[n], tag, MPI_COMM_WORLD, &rReq[n]);
		if (posi > maxPosi) { maxPosi = posi; }
	}
	if (cntAdjRank[m] > 0) {
		// Receive messages and unpack the transferred information.
		cnt = 0;
		cntRecvMsg = 0;
		memset(mpiTestRecvFlag, 0, sizeof(int) * cntAdjRank[m]);
		memset(mpiTestSendFlag, 0, sizeof(int) * cntAdjRank[m]);
		memset(chkL, 0, sizeof(int) * cntAdjRank[m]);

		nActCpPre = nActCp;
		nAbpCpPre = nAbpCp;
		while(cntRecvMsg < cntAdjRank[m]) {
			if (mpiTestSendFlag[cnt] == 0) {
				MPI_Test(&sReq[cnt], &mpiTestSendFlag[cnt], &status);
			}
			if (mpiTestRecvFlag[cnt] == 0 && mpiTestSendFlag[cnt] != 0 
					&& chkL[cnt] == 0) {
				MPI_Test(&rReq[cnt], &mpiTestRecvFlag[cnt], &status);
				if (mpiTestRecvFlag[cnt] != 0) { 
					posi = 0;
			  		MPI_UNPACK_INT(&nActCpSub, 1, cnt); 
			  		MPI_UNPACK_INT(&nAbpCpSub, 1, cnt);
			  		MPI_UNPACK_INT(&tempInt, 1, cnt);
					if (tglActFormDyn != 0) { 
						MPI_UNPACK_INT(&cntNucAss, 1, cnt);
					}
					// actin
			  		MPI_UNPACK_INT(&act.id[nActMe + nActCp], 
							nActCpSub, cnt);
			  		MPI_UNPACK_DBL(&P2(act.r,nActMe + nActCp,0), 
							nActCpSub * NDIM, cnt);
					if (modeActCh == 0) {
			  			MPI_UNPACK_INT(&P2A(act.ch,nActMe + nActCp,0,nChAc), 
								nActCpSub * nChAc, cnt);
					}
					else {
						memset(&P2A(act.ch,nActMe + nActCp,0,nChAc), -1, 
								sizeof(int) * nActCpSub * nChAc);
						for(k = 0; k < nActCpSub; k++) {
							pArr = &P2A(act.ch,nActMe + nActCp + k,0,nChAc);
				  			MPI_UNPACK_INT(&actCh.c, 1, cnt);
				  			MPI_UNPACK_INT(actCh.l, actCh.c, cnt);
							V2COPY(pArr, actCh.l);
							for(l = 2; l < actCh.c; l += 2) {
								pArr[actCh.l[l]] = actCh.l[l + 1];
							}
						}
					}
			  		MPI_UNPACK_INT(&act.iF[nActMe + nActCp], nActCpSub, cnt);
					// ABP
			  		MPI_UNPACK_INT(&abp.id[nAbpMe + nAbpCp], nAbpCpSub, cnt);
			  		MPI_UNPACK_DBL(&P2(abp.r,nAbpMe + nAbpCp,0), 
							nAbpCpSub * NDIM, cnt);
					if (motSA.gTgl == 0) {
						for(k = 0; k < nAbpCpSub; k++) {
				  			MPI_UNPACK_INT(&P2A(abp.ch,nAbpMe + nAbpCp + k,0,
								nChAb),	nChAb - 2, cnt);
							V2SET_ALL(&P2A(abp.ch,nAbpMe + nAbpCp + k,3,nChAb), 
									-1);
						}
					}
					else {
			  			MPI_UNPACK_INT(&P2A(abp.ch,nAbpMe + nAbpCp,0,nChAb), 
								nAbpCpSub * nChAb, cnt);
		  				MPI_UNPACK_INT(&abp.mId[nAbpMe + nAbpCp], nAbpCpSub, 
								cnt);
					}
					if (tglActFormDyn != 0) {
						if (cntNucAss > 0) {
							for(k = 0; k < nActCpSub; k++) {
								InsertElement1dArrayWoChk(insNeiPar.l, 
									&insNeiPar.c, act.id[nActMe + nActCp + k]);
							}
						}
					}
				
					nActCp += nActCpSub;
					nAbpCp += nAbpCpSub;
					cntRecvMsg++;
					chkL[cnt] = 1;
				}
			}
			cnt++;
			if (cnt == cntAdjRank[m]) { cnt = 0; }
		}
		for(n = nActCpPre; n < nActCp; n++) {
			iAct[act.id[n + nActMe]] = n + nActMe;
		}
		for(n = nAbpCpPre; n < nAbpCp; n++) {
			iAbp[abp.id[n + nAbpMe]] = n + nAbpMe;
		}
	}
  }
  free(chkL);
  if (modeActCh != 0) { free(actCh.l); }
  
  if (maxPosi > sizeBufMsg / 2) {
	end = cntAdjRank[0];
	sizeBufMsg *= 2;
	for(n = 0; n < end; n++) {
		free(bufSendMsg[n]);
		free(bufRecvMsg[n]);
		MALLOC(bufSendMsg[n],char,sizeBufMsg);
		MALLOC(bufRecvMsg[n],char,sizeBufMsg);
	}	
  }
  for(n = 0; n < insNeiPar.c; n++) {
	ind[0] = insNeiPar.l[n];
	if (ind[0] < nAct) {
		ind2 = iAct[ind[0]];
		pArr = &P2A(act.ch,ind2,0,nChAc);
		for(k = 0; k < 2; k++) {
			CONT(pArr[k] < 0);
			CONT(iAct[pArr[k]] < 0);
			CS = FindElementArray(act.cyl.l, act.cyl.c, ind[0], k, 2);
			CONT(CS > -1);
			InsertElementInNeighborList(ind[0], pArr[k], k);
		}
	}
	else {
		ind[0] -= nAct;
		InsertElementInNeighborList(ind[0], -1, 2);
	}	
  }
}

/*--------------------- Copy particles between subdomains --------------------*/

/*-------------------- Process conflicts and dynamic events ------------------*/

void UpdateAbpUnbRebLists(int abpInd, int actInd, int side, int mode) {
  V3SET(&P2A(sendAbpDyn.l,sendAbpDyn.c,0,3), abpInd, actInd, side);
  (sendAbpDyn.c)++;
  V3SET(&P2A(noAbpDyn.l,noAbpDyn.c,0,3), abpInd, actInd, currTimeStep);
  (noAbpDyn.c)++;
}

// Subroutine for UpdateActinAbpDynamicsEventsSubroutine()
void UpdateActinAbpDynamicsEventsSubSubroutine(ListInt *all, int abpInd, 
		int actInd, int side, ListInt *confAct, int *CS) {
  int locActInd, actInd2, locActInd2, locAbpInd, abpInd2, locAbpInd2, loc[2];
  int abpRankMol, abpRankMol2, actRankMol, abpSide, *pArr, *pArr2, CS2;

  locActInd = iAct[actInd];
  locAbpInd = iAbp[abpInd];
  abpInd2 = P2A(act.ch,locActInd,side,nChAc);
  actRankMol = CalcRankMolecule(&P2(act.r,locActInd,0));
  // It is possible that multiple ABPs belonging to different subdomains 
  // try to bind to the same actin. Check the possibility here.
  // If so, only ABP in the same subdomain as that of actin can bind to actin.
  if (abpInd2 > -1 && abpInd2 != abpInd) {
	CS2 = Find2ElementArray(confAct->l, confAct->c, locActInd, side, 0, 2);
	if (CS2 == -1) {
		V2SET(&P2A(confAct->l,confAct->c,0,2), locActInd, side);
		(confAct->c)++;
	}
	//P2A(confAct,locActInd,side - 2,2) = 1;
	locAbpInd2 = iAbp[abpInd2];
	if (locAbpInd2 > -1) {
		pArr2 = &P2A(abp.ch,locAbpInd2,0,nChAb); 
		abpRankMol2 = CalcRankMolecule(&P2(abp.r,locAbpInd2,0));
		// If they belong to different subdomains
		if (abpRankMol2 != actRankMol) {
			CS2 = 1;
			if (pArr2[0] == actInd) { abpSide = 0; }
			else if (pArr2[1] == actInd) { abpSide = 1; }
			else { 
				CS2 = -1;
			}
			if (CS2 > -1) {
				// Sever a link between actin and the pre-existing ABP
				UpdateActinDisassemblySubroutine2(abpInd2, actInd);
				// Check whether the pre-existing ABP originates from walking
				CS2 = FindElementArray(all->l, all->c, abpInd2, 0, 3);
			}
			if (CS2 > -1 && CS2 < all->c - 1) {
				pArr = &P2A(all->l,CS2,0,3);
				// If it originates from walking, the pre-existing ABP should
				// return to the previous actin.
				loc[0] = (int)((P2A(pArr,0,2,3) - 2) / nChAcY);
				loc[1] = (int)((P2A(pArr,1,2,3) - 2) / nChAcY);
				if (abpInd2 == P2A(pArr,1,0,3) 
						&& ((P2A(pArr,0,1,3) != P2A(pArr,1,1,3) 
						&& loc[0] == nChAcX - 1 && loc[1] == 0) 
						|| (P2A(pArr,0,1,3) == P2A(pArr,1,1,3) 
						&& loc[1] - loc[0] == 1))) { 
					// actInd2 is the previous actin to return.
					actInd2 = P2A(pArr,0,1,3);
					locActInd2 = iAct[actInd2];
					if (locActInd2 > -1) {
						P2A(act.ch,locActInd2,P2A(pArr,0,2,3),nChAc) = abpInd2;
					}
					if (abpSide == 0) {
						pArr2[1] = pArr2[0];
					}
					pArr2[abpSide] = actInd2;
					if (locAbpInd2 < nAbpMe) {
						if (pArr2[0] > -1 && pArr2[1] > -1) {
							if (tglNeiAbpSC != 0 && tglNeiAbpDC == 0) {
								DeleteElementInNeighborList(abpInd2, 1);
							}
						}
						(motWalk.cntMe)--;
						if (pArr2[1] < 0) { (motInaUnb.cntMe)--; } 
						else { (motUnb.cntMe)--; }
					}
					V3SET(&P2A(noAbpDyn.l,noAbpDyn.c,0,3), 
							abpInd2, actInd2, currTimeStep);
					(noAbpDyn.c)++;
					if (actInd != actInd2
							&& ((locActInd2 > -1 && locActInd2 < nActMe) 
							|| (locAbpInd2 > -1 && locAbpInd2 < nAbpMe))) {
						InsertLongChain(abpInd2 + nAct, actInd2, 
								minDimDomC * 0.9);
					}
				}
			}
			// if (abpRankMol2 != actRankMol)
			// and if information of the current ABP is available
			if (locAbpInd > -1) {
				abpRankMol = CalcRankMolecule(&P2(abp.r,locAbpInd,0));
				*CS = (abpRankMol == actRankMol) ? 1 : 0;
			}
			// if information of the current ABP is not available
			else { *CS = 1; }
		}
		// if (abpRankMol2 == actRankMol)
		else { *CS = 0; }
	}
  }
  // If there was a conflict on the binding spot, but if it is currently
  // avalable
  else if (abpInd2 < 0) {
	CS2 = Find2ElementArray(confAct->l, confAct->c, locActInd, side, 0, 2);
	if (CS2 > -1) {
		if (locAbpInd > -1) {
			abpRankMol = CalcRankMolecule(&P2(abp.r,locAbpInd,0));
			*CS = (actRankMol == abpRankMol) ? 1 : 0;
		}
		else { *CS = 1; }
	}
	else { *CS = 1; }
  }
}

// Subroutine for UpdateActinAbpDynamicsEvents()
void UpdateActinAbpDynamicsEventsSubroutine(ListInt *all, int mode) {
  int n, k, side, walk, begin, CS, CS2, *pArr, *pArr2, loc[2];
  int abpInd, actInd, actInd2, locAbpInd, locActInd, nextActInd, locNextActInd;
  ListInt confAct;

  confAct.l = allIntL;
  confAct.c = 0;
  begin = sendAbpDyn.c + ((mode == 0) ? sendActDyn.c : 0);
  for(n = begin; n < all->c; n++) {
	// ABP unbinding/binding/walking
	if (P2A(all->l,n,2,3) >= 0) {
		abpInd = P2A(all->l,n,0,3);
		actInd = P2A(all->l,n,1,3);
		locAbpInd = iAbp[abpInd];
		locActInd = iAct[actInd];
		// Check whether it is a walking action or not.
		// If the first element is the same in two consecutive rows,
		// it indicates the walking.
		walk = 0;
		CS2 = 1;
		if (n < (all->c) - 1 && motWalk.tgl != 0) {
			if (abpInd == P2A(all->l,n + 1,0,3) && P2A(all->l,n + 1,2,3) >= 0) {
				loc[0] = (int)((P2A(all->l,n,2,3) - 2) / nChAcY);
				loc[1] = (int)((P2A(all->l,n + 1,2,3) - 2) / nChAcY);
				if ((actInd != P2A(all->l,n + 1,1,3) && loc[0] == nChAcX - 1
						&& loc[1] == 0) || (actInd == P2A(all->l,n + 1,1,3)
						&& loc[1] - loc[0] == 1)) {
					nextActInd = P2A(all->l,n + 1,1,3);
					walk = 1;
					if (iAct[nextActInd] > -1) {
						side = P2A(all->l,n + 1,2,3);
						// Check whether a conflict exists on the binding spot 
						// to walk.
						UpdateActinAbpDynamicsEventsSubSubroutine(all, abpInd,
								nextActInd, side, &confAct, &CS2);
					}
				}
			}
		}
		if (locAbpInd >= nAbpMe || locAbpInd < 0) {
			CS = 1;
			if (locActInd > -1 && !(CS2 == 0 && walk == 1)) {
				pArr = &P2A(act.ch,locActInd,0,nChAc);
				side = P2A(all->l,n,2,3);
				UpdateActinAbpDynamicsEventsSubSubroutine(all, abpInd, 
						actInd, side, &confAct, &CS);
				// Modify act.ch
				if (CS == 1) {
					if (locActInd < nActMe) {
						if (pArr[side] == abpInd) {
							DeleteLongChain(abpInd + nAct, actInd);
						}
						else {
							InsertLongChain(abpInd + nAct, actInd, 
									minDimDomC * 0.9);
						}
					}
					pArr[side] = (pArr[side] == abpInd) ? -1 : abpInd;
				}
			}	
			if (locAbpInd >= nAbpMe) {
				// Modify abp.ch
				pArr2 = &P2A(abp.ch,locAbpInd,0,nChAb);
				// If it is a successful walking motion
				if (CS2 == 1 && walk == 1) {
					pArr2[(pArr2[0] == actInd) ? 0 : 1] = nextActInd;
				}
				// If it is a simple binding or unbinding
				else if (CS == 1 && walk == 0) {
					if (pArr2[0] == actInd) {
				       	pArr2[0] = pArr2[1];
				       	pArr2[1] = -1;	
					}
					else if (pArr2[1] == actInd) {
						pArr2[1] = -1; 
					}
					else { 
						pArr2[(pArr2[0] < 0) ? 0 : 1] = actInd; 
					}
				}
			}
			V3SET(&P2A(noAbpDyn.l,noAbpDyn.c,0,3), abpInd, actInd, 
					currTimeStep);
			(noAbpDyn.c)++;
		}
		// From actin dissembly or severing
		// (Among ABP dynamics, unbinding of ABPs by actin disassembly and
		// severing is only a thing which can occur at the side of actin.
		if (locAbpInd > -1 && locAbpInd < nAbpMe && locActInd > -1) {
			pArr2 = &P2A(abp.ch,locAbpInd,0,nChAb);
			// By chance, ABP might have been unbound already by its unbinding 
			// or walking at the same time step, so it needs to check the bond.
			if (pArr2[0] == actInd || pArr2[1] == actInd) {
				// Sever the chain
				UpdateActinDisassemblySubroutine2(abpInd, actInd);
			}
		}
		// If walking event
		if (walk == 1) {
			if (CS2 == 1) {
				locNextActInd = iAct[nextActInd];
				if (locNextActInd > -1) {
					side = P2A(all->l,n + 1,2,3);
					P2A(act.ch,locNextActInd,side,nChAc) = abpInd;
					if (locNextActInd < nActMe) {
						if (actInd != nextActInd) {
							InsertLongChain(abpInd + nAct, nextActInd, 
									minDimDomC * 0.9);
						}
					}
					V3SET(&P2A(noAbpDyn.l,noAbpDyn.c,0,3), 
							abpInd, nextActInd, -1 * currTimeStep);
					(noAbpDyn.c)++;
				}
			}	
			n++; 
		}
	}
	// Actin disassembly
	else if (P2A(all->l,n,2,3) == -1 || P2A(all->l,n,2,3) == -2) {
		actInd = P2A(all->l,n,0,3);
		actInd2 = P2A(all->l,n,1,3);
		side = -1 * P2A(all->l,n,2,3) - 1;
		locActInd = iAct[actInd];
		pArr = &P2A(act.ch,locActInd,0,nChAc);
		// This must be among copied particles
		if (locActInd > -1) { 
	        if (pArr[0] > -1 || pArr[1] > -1) {
				// Check whether or not ABPs are bound on it
				if (side == 1) {
					for(k = 2; k < nChAc; k++) {
						abpInd = pArr[k];
						CONT(!(abpInd > -1));
						UpdateActinDisassemblySubroutine2(abpInd, actInd);
					}
				}
		        V3SET_ALL(&P2(act.r,locActInd,0), 0.);
		        SetAllValue1dArrayInt(pArr, nChAc, -1);
				act.iF[locActInd] = -1;
			}
		  	V2SET(&P2A(noActDyn.l,noActDyn.c,0,2), actInd, currTimeStep);
		  	(noActDyn.c)++;
		}
		// This can be among copied particles or among current subdomain
		if (iAct[actInd2] > -1) { 
			pArr2 = &P2A(act.ch,iAct[actInd2],0,nChAc);
	        if (pArr2[0] > -1 || pArr2[1] > -1) {
				// Check whether or not ABPs are bound on it
				if (side == 0) {
					for(k = 2; k < nChAc; k++) {
						abpInd = pArr2[k];
						CONT(!(abpInd > -1));
						UpdateActinDisassemblySubroutine2(abpInd, actInd2);
					}
				}
				UpdateActinDisassemblySubroutine(actInd2, side);
			}
		}
		// Update the neighboring list
		// Find and delete an actin cylinder corresponding the disassembled one
		//DeleteElementInNeighborList(P2A(all->l,n,0,3), 0);
		DeleteElementInNeighborList(actInd, 0);
	}
  }
  CheckArraySize(&sendAbpDyn, &sendAbpDyn.siz, 3, 0);
  sendAbpDyn.c = 0;
  if (mode == 0) {
	CheckArraySize(&sendActDyn, &sendActDyn.siz, 3, 0);
  }
  sendActDyn.c = 0;
}

// This function let other subdomains know the happening of the dynamic 
// behaviors of ABPs (e.g. unbinding, binding, and walking). 
void UpdateActinAbpDynamicsEvents(int mode) {
  int n, sizeArr;
  ListInt all;

  sizeArr = sendAbpDyn.siz + ((mode == 0) ? sendActDyn.siz : 0.);
  sizeArr *= cntAdjRank[0];
  MALLOC(all.l,int,sizeArr);

  all.c = sendAbpDyn.c;
  for(n = 0; n < sendAbpDyn.c; n++) {
	V3COPY(&P2A(all.l,n,0,3), &P2A(sendAbpDyn.l,n,0,3));
  }
  if (mode == 0) {
	(all.c) += sendActDyn.c;
	for(n = 0; n < sendActDyn.c; n++) {
		V3COPY(&P2A(all.l,n + sendAbpDyn.c,0,3), &P2A(sendActDyn.l,n,0,3));
	}
  }
  CollectArrayIntFromAdjacentSubdomain(&all, 3);
  UpdateActinAbpDynamicsEventsSubroutine(&all, mode);
  free(all.l);
}

/*-------------------- Process conflicts and dynamic events ------------------*/

/*-------------------- Balancing counters between subdomains -----------------*/

// Fill gap by information of transferred elements
// mode = 0: actin, 1: ACP, 2: motor
void UpdateActinAbpMonomerListSubroutine2(int locInd, int ind, int mode) {
  if (mode == 0) {
	V3SET_ALL(&P2(act.r,locInd,0), 0.);
	SetAllValue1dArrayInt(&P2A(act.ch,locInd,0,nChAc), nChAc, -1);
	act.id[locInd] = ind;
	act.iF[locInd] = -1;
	act.fix[locInd] = -1;
	if (rho.tgl != 0) { act.rho[locInd] = -1; }
	iAct[ind] = locInd;
	if (confVmdInfo.tgl != 0) {
		recAct.len[locInd] = 0.;
		recAct.sprF[locInd] = 0.;
		recAct.bendF[locInd] = 0.;
		V4SET_ALL(&P2A(recAct.allF,locInd,0,NDIM + 1), 0.);
		recAct.cnt[locInd] = 0;
	}
  }
  else {
	V3SET_ALL(&P2(abp.r,locInd,0), 0.);
	SetAllValue1dArrayInt(&P2A(abp.ch,locInd,0,nChAb), nChAb, -1);
	P2A(abp.ch,locInd,2,nChAb) = mode - 1;
	if (rho.tgl != 0) { abp.rho[locInd] = -1; }
	if (motSA.gTgl != 0) {
		abp.mId[locInd] = -1;
	}
	abp.id[locInd] = ind;
	iAbp[ind] = locInd;

	if (confVmdInfo.tgl != 0) {
		SetAllValue1dArrayDouble(&P2A(recAbp.len,locInd,0,recAbp.nL), 
				recAbp.nL, 0.);
		recAbp.sprF[locInd] = 0.;
		recAbp.bendF[locInd] = 0.;
		V4SET_ALL(&P2A(recAbp.allF,locInd,0,NDIM + 1), 0.);
		recAbp.cnt[locInd] = 0;
	}
	if (recLongF.tgl != 0) {	
		V4SET_ALL(&P2A(recLongSprFabp,locInd,0,4), 0.);
	}
  }
}

// Shift all the information of copied elements from locInd2 to locInd1
// mode = 0: actin, 1: ABP
void UpdateActinAbpMonomerListSubroutine(int locInd1, int locInd2, int mode) {
  if (mode == 0) {
	V3COPY(&P2(act.r,locInd1,0), &P2(act.r,locInd2,0));
	Copy1dArrayInt(&P2A(act.ch,locInd1,0,nChAc), 
			&P2A(act.ch,locInd2,0,nChAc), nChAc);
	act.iF[locInd1] = act.iF[locInd2];
	act.id[locInd1] = act.id[locInd2];
	iAct[act.id[locInd2]] = locInd1;
  }
  else {
	V3COPY(&P2(abp.r,locInd1,0), &P2(abp.r,locInd2,0));
	Copy1dArrayInt(&P2A(abp.ch,locInd1,0,nChAb), 
			&P2A(abp.ch,locInd2,0,nChAb), nChAb);
	if (motSA.gTgl != 0) {
		abp.mId[locInd1] = abp.mId[locInd2];
	}
	abp.id[locInd1] = abp.id[locInd2];
	iAbp[abp.id[locInd2]] = locInd1;
  }
}

/*-------------------- Balancing counters between subdomains -----------------*/

/*--------------- Collect information from adjacent subdomains ---------------*/
void CollectArrayIntFromAdjacentSubdomain(ListInt *all, int col) {
  int m, n, posi, rep, tag = 0, idRank;
  int cnt, cntRecvMsg, cntElem, *chkList;

  MALLOC(chkList,int,cntAdjRank[0]);
  rep = 1;
  for(m = 0; m < rep; m++) {
	for(n = 0; n < cntAdjRank[m]; n++) {
		idRank = iRank[n];
		posi = 0;
		MPI_PACK_INT(&(all->c), 1, n);
		MPI_PACK_INT(all->l, (all->c) * col, n);
		MPI_Isend(bufSendMsg[n], posi, MPI_PACKED, idRank, 
				tag, MPI_COMM_WORLD, &sReq[n]);
		MPI_Irecv(bufRecvMsg[n], sizeBufMsg, MPI_PACKED, idRank,
				tag, MPI_COMM_WORLD, &rReq[n]);
	}
    if (cntAdjRank[m] > 0) {
        cnt = 0;
        cntRecvMsg = 0;
        memset(mpiTestRecvFlag, 0, sizeof(int) * cntAdjRank[m]);
        memset(mpiTestSendFlag, 0, sizeof(int) * cntAdjRank[m]);
        memset(chkList, 0, sizeof(int) * cntAdjRank[m]);
        while(cntRecvMsg < cntAdjRank[m]) {
            if (mpiTestSendFlag[cnt] == 0) {
                MPI_Test(&sReq[cnt], &mpiTestSendFlag[cnt], &status);
            }
            if (mpiTestRecvFlag[cnt] == 0 && mpiTestSendFlag[cnt] != 0
                    && chkList[cnt] == 0) {
                MPI_Test(&rReq[cnt], &mpiTestRecvFlag[cnt], &status);
                if (mpiTestRecvFlag[cnt] != 0) {
                    posi = 0;
					MPI_UNPACK_INT(&cntElem, 1, cnt)
					MPI_UNPACK_INT(&P2A(all->l,all->c,0,col), 
							cntElem * col, cnt);
					(all->c) += cntElem;
                    cntRecvMsg++;
                    chkList[cnt] = 1;
				}
			}
            cnt++;
            if (cnt == cntAdjRank[m]) { cnt = 0; }
		}
	}
  }
  free(chkList);
}

void CollectArrayDblFromSubdomainList(double *sendData, double *recvData, 
		int size, ListInt *subdList, int gotData) {
  int n, posi, tag = 0, *flag, cnt, cntRecvMsg, cntMsg, sizBufMsg;
  char *bufMsg, **bufMsg2;
  MPI_Request *req;

  // Gather the local sums to main CPU
  sizBufMsg = sizeof(double) * size;
  if (rank != subdList->l[0] && gotData != 0) {
    MALLOC(bufMsg,char,sizBufMsg);
	posi = 0;
	MPI_Pack(sendData, size, MPI_DOUBLE, bufMsg, sizBufMsg,
			&posi, MPI_COMM_WORLD);
	MPI_Send(bufMsg, posi, MPI_PACKED, subdList->l[0], tag, MPI_COMM_WORLD);
    free(bufMsg);
  }
  else if (rank == subdList->l[0]) {
	cntMsg = (subdList->c) - 1;
	for(n = 0; n < size; n++) {
		recvData[n] = sendData[n];
	}
    MALLOC(flag,int,cntMsg);
    MALLOC(req,MPI_Request,cntMsg);
    MALLOC2(bufMsg2,char,cntMsg);
    for(n = 0; n < cntMsg; n++) {
        MALLOC(bufMsg2[n],char,sizBufMsg);
        MPI_Irecv(bufMsg2[n], sizBufMsg, MPI_PACKED, subdList->l[n + 1], tag,
                MPI_COMM_WORLD, &req[n]);
    }
    cnt = 0;
    cntRecvMsg = 0;
    memset(flag, 0, sizeof(int) * cntMsg);
    while(cntRecvMsg < cntMsg) {
        if (flag[cnt] == 0) {
            MPI_Test(&req[cnt], &flag[cnt], &status);
            if (flag[cnt] != 0) {
                posi = 0;
                MPI_Unpack(bufMsg2[cnt], sizBufMsg, &posi, 
						&recvData[size * (cntRecvMsg + 1)], 
						size, MPI_DOUBLE, MPI_COMM_WORLD);
                cntRecvMsg++;
            }
        }
        cnt++;
        if (cnt == cntMsg) { cnt = 0; }
    }
    for(n = 0; n < cntMsg; n++) {
		free(bufMsg2[n]);
	}
	free(bufMsg2);
	free(req);
	free(flag);
  }
}

/*--------------- Collect information from adjacent subdomains ---------------*/
