// ##################################################
// #   boundary.c - finally revised on Apr 2022     #
// #   coded by Taeyoon Kim                         #
// #   Copyright (C) 2005 - 2022, Taeyoon Kim,      #
// #   All rights reserved.                         #
// ##################################################
// This file contains functions related to boundaries or domain

// Check whether the chain crosses a boundary or not.
void CheckCrossBound(int *sft, double *dr, double *dr2) {
  int k;
 
  FOR_NDIM(k) {
	CONT(pbc[k] != 1);
	CONT(!(fabs(dr2[k] - dr[k]) > 1.5 * dimDomH[k]));
	sft[k] += (dr2[k] > dr[k]) ? -1 : 1;
  }
}

// Apply periodic boundary condition (PBC) to a chain vector.
// If the absolute value of a component of the vector is greater than half of 
// domain width with PBC, it is added or subtracted by domain width.
void ApplyBoundCondVecDiff(double *dr) {
  int k;
  FOR_NDIM(k) {
	CONT(pbc[k] != 1);
	if (dr[k] >= dimDomH[k]) { dr[k] -= dimDom[k]; }
	else if (dr[k] <= REVSIGN(dimDomH[k])) { dr[k] += dimDom[k]; }
  }
}

// If mode is -1, neglect repulsive condition.
// Otherwise, if mode is 0, actin. If mode is 1-3, ABP.
void ApplyBoundCondVector(double *r, int mode, int mode2) {
  int k;
  double disp, dispAll, repF, drag, sft;

  FOR_NDIM(k) {
	disp = 0.;
	// If PBC
    if (pbc[k] == 1) {
		if (r[k] < disp + rGrid[k][0])
		{ r[k] += dimDom[k]; }
		else if (r[k] >= disp + rGrid[k][nGrid[k] - 1])
		{ r[k] -= dimDom[k]; }
	}
	// If no PBC and mode >= 0, repulsive condition is applied.
	else if (pbc[k] != 1 && mode >= 0) {
		sft = 0.;
		if (mode == 0) { drag = 1.; }
		else if (mode >= 1 && mode <= 3) { drag = abpF.drag[mode - 1].inv; }
		if ((r[k] < disp + rGrid[k][0] + sft) 
				|| (r[k] >= disp + rGrid[k][nGrid[k] - 1] - sft)) {
			if (r[k] < disp + rGrid[k][0] + sft) {
				dispAll = disp + rGrid[k][0] + sft;
			}
			else {
				dispAll = disp + rGrid[k][nGrid[k] - 1] - sft;
			}
			repF = bnd.stfRep * (r[k] - dispAll);
			r[k] -= repF * drag * dt;
		}
	}
  }
}

// Apply the periodic boundary condition on the positions of whole things
void ApplyBoundCondAll(void) {
  int n;

  FOR_ACTME(n) {
	CONT(ISACTM(n));
	ApplyBoundCondVector(&P2(act.r,n,0), 0, (act.fix[n] > -1 ? 0 : 1));
  }
 
  FOR_ABPME(n) {
	CONT(ISABPIM(n));
	ApplyBoundCondVector(&P2(abp.r,n,0), K_ABP(n) + 1, 1);
  }
}

// Regardless of shear deformation, offset positions of all particles into 
// retangular-solid domain
// If mode = 0, adjust the position only in shearing direction
// If mode > 0, apply PBC in other directions.
void ConvertRectDomainVector(double *r, int mode) {
  int k;

  if (mode == 0) {
  }
  else {
	FOR_NDIM(k) {
		if (pbc[k] != 0) {
			if (r[k] >= rGrid[k][nGrid[k] - 1]) { r[k] -= dimDom[k]; }
			else if (r[k] < rGrid[k][0]) { r[k] += dimDom[k]; }
		}
	}
  }
}

int CheckParticleInDomain(double *r) {
  int CS, k;
  double disp;
  CS = 1;
  FOR_NDIM(k) {
	CONT(!(pbc[k] == 0));
	disp = 0.;
	CONT(!(r[k] < disp + rGrid[k][0] || r[k] >= disp + rGrid[k][nGrid[k] - 1]));
	CS = 0;
	break;
  }
  return CS;
}

