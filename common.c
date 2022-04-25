// ##################################################
// #   common.c - finally revised on Apr 2022       #
// #   coded by Taeyoon Kim                         #
// #   Copyright (C) 2005 - 2022, Taeyoon Kim,      #
// #   All rights reserved.                         #
// ##################################################
// This file mainly contains the definitions of types and variables.

typedef struct {
  int **l, c, siz;
} ListInt2;
typedef struct {
  int *l, c, siz, cOri, cM, cMOri;
} ListInt;
typedef struct {
  double *l;
  int c;
} ListDbl;
typedef struct {
  int *cntU, *cntB, *cntW;
} AbpDyn;
// Duration
typedef struct {
  long long dur;
  double durR, facK;
} Duration;
typedef struct {
  int tgl, prd, nReg[3], nActiv, nTotReg, *reg;
  int cntActFor, cntActAss, cntMotNuc, cntMotAss, nFor, nForMe;
  int dlyMot, durMot, dlyAct, durAct;
  double prdR, facAct[2], facM, *widReg; 
  double dlyMotR, durMotR, dlyActR, durActR;
  double RFor, frac, fracM;
} Rho;
typedef struct {
  int tgl, gTgl,cntMe;
  double k, facP;
} Nucle;
typedef struct {
  int tgl, gTgl, cntMe;
  double k[2], p[2], pWA[2], facP[2], facKWA;
} AssDis;
// Dynamic behaviors of ABPs
typedef struct {
  int tgl, gTgl, cntMe;
  double facP, k, facK;
} InaUnbMbind;
typedef struct {
  int tgl, gTgl, gTglPa, gTglCrsAng, cntMe;
  double p, facK, por;
} Reb;
typedef struct {
  int tgl, gTgl, cntMe, maxF;
  double *p, facK0, facX, facK0c, facXc;
} Unb;
typedef struct {
  int tgl, gTgl, cntMe, maxF[2];
  double *p, facK0, facX;
} MotUnbWalk;
typedef struct {
  int *cnt, nL;
  double *len, *allF, *sprF, *bendF;
} RecLenForce;
typedef struct {
  int prd, tgl, gTgl, mode;
  double prdR;
} FuncCont;
typedef struct {
  int prd, tgl, gTgl, tgl2, gTglBnd, gTglInfo, mode;
  double prdR, minF, maxF;
} RecConfVmd;
typedef struct {
  double *r, *v, *f, *fBr, *maxFbr, *rPrev;
  int *ch, *id, *fix, *iF, *rho;
  ListInt cyl;
} MoleAct;
typedef struct {
  double *r, *v, *f, *fBr, *maxFbr;
  int *ch, *id, *mId, *rho, *kind;
} MoleAbp;
typedef struct {
  double n, inv;
} Double;
typedef struct {
  int gTgl;
  double stf, stf2, eq, eq2, lo, hi, facStf, drag, *val;
} Force;
typedef struct {
  int gotF;
  double stfRep, *r;
} Boundary;
typedef struct {
  double wid[3], base[3];
  int n[3], *l;
} Cell;
typedef struct {
  double *dia, *repDia, maxSiz;
  Double *len, *dtRepAct, *dtSprAct, *drag;
  Force *a90, *cr, *bend, *spr, rep;
} AbpF;
typedef struct {
  int gTgl;
  double *p, pf, k;
} MotTurnover;
typedef struct {
  Force bend, spr;
  MotTurnover to;
  int tgl, gTgl, *id;
  int cntNucMe, cntAssMe, cntTurnMe, nMotPerTF, nNucMe, nMotPerSide;
  int nNucMeOri, nNucMeMOri, nNucMeM;
  double cenDist, kAss, pAss;
} MotSelfAss;
typedef struct {
  int nHead;
  double k01, k10, k12, k21, k20;
} MotMechChem;
typedef struct {
  double dia, dragR;
  Force bend, spr, rep;
} ActF;

#include <mpi.h>
#include "header.h"
#include "boundary.c"
#include "calForce.c"
#include "dataMgr.c"
#include "error.c"
#include "gatPrint.c"
#include "init.c"
#include "paraProc.c"
#include "process.c"
#include "record.c"
#include "rng.c"
#include "tools.c"
#include "update.c"

// Variables related to time steps
long long currTimeStep;
double dt, dtReal;
time_t initTime;
Duration netForm, rheo;
// Given concentration
int nAbpGoal, nAbpDet[NK_ABP], nAbpMall[NK_ABP], nAbpGoalDet[NK_ABP];
double cAct, RAbp[NK_ABP];
// Number of particles
int nAct, nAbp, nMot; 
int nActMe, nActCp, nActMeCp, nActC;
int nAbpMe, nAbpCp, nAbpMeCp, nAbpC;
int nAcpInaMe, nMotInaMe, nAcpMme, nMotMme;
int nActGoal, nActMall, nActFilaMe;
int nActMin, nAbpMin; 
// Chain, position, forces, and lists of particles
ListInt actM, acpM, motM, iFilaP;
MoleAct act;
MoleAbp abp;
ActF actF;
AbpF abpF;
MotSelfAss motSA;
int *chAct, *chAbp, tglNeiAbpSC, tglNeiAbpDC;
int nChAc, nChAcX, nChAcY, nChAb, nActPerSeg;
double *rAct, *rAbp;
// Update neighboring list
ListInt neigh;
Cell cell;
double dispHeuSq, maxNeiLenSq;
// Arrays related to lists
ListInt sendAbpDyn, noAbpDyn, sendActDyn, noActDyn;
int *fixAct, *abpMotId;
// Domain (width, periodic boundary conditions, or repulsive force)
int pbc[NDIM], neiPbc[NDIM], dir2D;
double dimDom[NDIM], dimDomH[NDIM], minDimDomC;
Boundary bnd;
// For boundaries
double *rGridInit;

/*-------------------------- For parallel processing -------------------------*/
int nCpu;
int *iAct, *iAbp;
// Ranks of CPUs
int rank, *adjRank, *iRank, *cntAdjRank;
// For boundaries and indicies of subdomains
int nCell[NDIM], iCell[NDIM], nGrid[NDIM];
double edge[NDIM*2], **rGrid, neiEdge;
// Used in CopyParticles and MoveParticles
int *cntCpPar, *cntMvPar, modeActCh;
ListInt *cpPar, *mvPar, insNeiPar;
// Varilables related to messages
int sizeBufMsg, *mpiTestSendFlag, *mpiTestRecvFlag;
char **bufSendMsg, **bufRecvMsg;
MPI_Status status;
MPI_Request *sReq, *rReq;
// Variables related to longCh
ListInt longCh, longChExtMsg, longChIntMsg;
double *longChDist, maxDisp, maxActCh;
/*----------------------------------------------------------------------------*/

/*---------------- Dynamic behaviors of actin, ACP, and motor ----------------*/
// Dynamics of actin
AssDis actAss, actDis;
Nucle actNuc;
int tglActMoDyn, tglActFormDyn, gTglActDynNF;
int gTglActTherm, gTglAcpTherm, gTglMotTherm;
int durNoActDyn;
// Dynamics of ACPs and motors
Reb acpReb, motReb;
Unb acpUnb;
MotUnbWalk motUnb, motWalk;
InaUnbMbind acpInaUnb, motInaUnb, acpMoBind, motMoBind;
int tglAbpAcInaDyn, tglAbpInaMoDyn, tglAcpAcInaDyn, tglMotAcInaDyn;
int gTglAcpDynNF, gTglMotUnbRebNF, gTglMotWalkNF;
int gTglImpAcpM, gTglImpMotM;
int durNoAcpUnbReb, durNoMotUnbReb, durNoMotWalk;
ListDbl unbLog, bindLog, toLog;
MotMechChem motMC;
int cntNucAss;
// Rho activity
Rho rho;
/*----------------------------------------------------------------------------*/
 
/*---------------------------- For data recording ----------------------------*/
FuncCont recProg, recConf;
// For general records
char fnOut[80], dataFold[80];
// Record the accumulated chain lengths and forces
RecLenForce recAct, recAbp;
// Record configuration for VMD
RecConfVmd recConfVmd;
// Record the turnover of ACPs or motors 
int tglRecAbpTurn;
double *abpTurn;
// Record (longitudinal) forces of ACPs or motors
FuncCont recLongF;
double *recLongSprFabp, *recInstSprFabp;
// Record etc
FuncCont confVmdInfo, recFilaL, recE, recInfo;
FuncCont recAbpUnb, recAbpBind, recAbpTurn;
/*----------------------------------------------------------------------------*/

// Check errors
int stopSig;
double magUnstF;
// Misc.
int seed, *allIntL;
double *arrAcos, *allDblL;

