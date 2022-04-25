// ##################################################
// #   header.h - finally revised on Apr 2022       #
// #   coded by Taeyoon Kim                         #
// #   Copyright (C) 2005 - 2022, Taeyoon Kim,      #
// #   All rights reserved.                         #
// ##################################################
// This file contains the definitions of parameters, function prototypes
// variables, and functions.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>

/******************************************************************************/
/************************** Definition of Parameters **************************/
/******************************************************************************/

// General parameters and constants
#define NDIM 					3		
#define NK_ABP					3
#define KT						1.0			// Boltzmann Energy
#define KT_IN_J					4.142e-21	// Boltzmann Energy in J
#define VISCOSITY				0.8599e0	// Viscosity of water
#define PI						3.141593	
#define N_AVO					6.022e23
#define POS_LARGE_VALUE			1e9
#define NEG_LARGE_VALUE			-1e9
#define POS_SMALL_VALUE			1e-10		
// For neighboring list
#define DT_NL_UPDATE			0.4
#define DT_NL_UPDATE_BUF		0.1
// Forces
#define MAG_UNSTABLE_FORCE		10000.0e-12
// Record general data
#define OUTPUT_FILE				"Output"
#define DELETE_FILE				1		
// Extent of pre-assigned arrays 
#define DEG_ARR_ACOS			10000
// For parallelization (MoveParticles and CopyParticles)
#define OFFSET_LIST															\
	{{0,1,2,3}, {0,1,2,3,4}, {0,3,4}, {0,1,2,3}, {0,1,2,3,4}, {0,3,4}, 		\
	{0,1}, {0,1}, {0}, {0,1,2,3,5,6,7,8}, {0,1,2,3,4,5,6,7,8,9,10}, 		\
	{0,3,4,5,8,9,10}, {0,1,2,3,5,6,7,8,12,13},					 			\
	{0,1,2,3,4,5,6,7,8,9,10,11,12,13}, {0,3,4,5,8,9,10,11,12}, 		 		\
	{0,1,5,6,12,13}, {0,1,5,6,10,11,12,13}, {0,5,10,11,12}}
#define OFFSET_LEN	{4,5,3,4,5,3,2,2,1,8,11,7,10,14,9,6,8,5}
#define OFFSET_CPU_LIST														\
	{{0,1,3,4,9,10,12}, {1,4,10}, {1,2,4,5,10,11,14}, {3,4,12}, {4}, 		\
	{4,5,14}, {3,4,6,7,12,15,16}, {4,7,16}, {4,5,7,8,14,16,17}, {9,10,12},	\
	{10}, {10,11,14}, {12}, {-1}, {14}, {12,15,16}, {16}, {14,16,17},		\
	{9,10,12,18,19,21,22}, {10,19,22}, {10,11,14,19,20,22,23}, {12,21,22},	\
	{22}, {14,22,23}, {12,15,16,21,22,24,25}, {16,22,25}, 					\
	{14,16,17,22,23,25,26}}
#define OFFSET_CPU_LEN {7,3,7,3,1,3,7,3,7,3,1,3,1,1,1,3,1,3,7,3,7,3,1,3,7,3,7}

/*---------------- Dynamic behaviors of actin, ACP, and motor ----------------*/
// Unbinding/binding/binding of ACPs
#define MAX_ACP_UNB_FORCE   500.0e-12
#define K0_ACP_UNB			0.115
#define X_ACP_UNB			104.0e-12
#define K_ACP_BIND			1.0e2
// Binding of motors
#define K_MOT_BIND			1.0e8
// Turnover of motor filaments
#define MAX_MOT_TURN_FORCE   500.0e-12
// Lower and upper limits of distances allowing binding of ACPs or motors
#define DTL_ACP_REB			0.9	
#define DTH_ACP_REB			1.1	
#define DTL_MOT_REB			0.7037	
#define DTH_MOT_REB			1.2963	
// Time for which a dynamic behavior is prohibited after a previous one occurs
#define DUR_ACP_NO_UNB		2e-9    
#define DUR_MOT_NO_UNB		2e-9   
#define DUR_MOT_NO_WALK		2e-9    
#define DUR_ACT_NO_DYN		2e-9    
/*----------------------------------------------------------------------------*/

/*-- Mechanical stiffness and geometric parameters of actin, ACP, and motor --*/

// Mechanical stiffness for actin
#define STF_ACT_SPR			1.691e-1
#define STF_ACT_BEND		3.697e-26
#define STF_ACT_REP			1.691e-3
#define STF_ABP_REP			1.691e-3
// Geometric parameters for actin
#define ANGL_ACT_BEND		10.
#define DIA_CYL_ACT			7.0e-9
#define DIA_CYL_ACT_UM		(DIA_CYL_ACT * 1e6)
#define DIA_CYL_ACT_NM		(DIA_CYL_ACT * 1e9)
// Mechanical stiffness for ACP crosslinking filaments at right angle
#define STF_ACPC_SPR		4.2270e-4	
#define STF_ACPC_BEND		1.036e-18	 
#define STF_ACPC_90			4.142e-18	
#define STF_ACPC_TOR		4.142e-18
// Geometric parameters for ACP crosslinking filaments at right angle
#define L_ACPC_ARM			35.0e-9
#define ANG_ACPC_BEND		0
#define ANGL_ACPC_BEND		15.
#define ANGL_ACPC_90		10.
#define ANGL_ACPC_TOR		10.
#define DIA_CYL_ACPC		10.0e-9
// Mechanical stiffness for ACP bundling filaments in parallel
#define STF_ACPB_SPR		2.0e-3
#define STF_ACPB_BEND		1.036e-19
#define STF_ACPB_90			1.036e-19
#define STF_ACPB_TOR		4.142e-18
// Geometric parameters for ACP bundling filaments in parallel
#define L_ACPB_ARM			20.0e-9
#define ANG_ACPB_BEND		0.
#define ANGL_ACPB_BEND		15.
#define ANGL_ACPB_90		10.
#define ANGL_ACPB_TOR		10.
#define DIA_CYL_ACPB		10.0e-9
// Mechanical stiffness for motors
#define STF_MOT_SPR			1.0e-3
#define STF_MOT_SPR2		1.0e-3
#define STF_MOT_BEND		1.036e-18	 
#define STF_MOT_90			0
#define STF_MOT_TOR			4.142e-18
#define MOT_STALLF			5.04e-12
// Geometric parameters for motors
#define L_MOT_ARM			10.0e-9	
#define ANG_MOTBACK_TOR		180.
#define ANG_MOTBACK_BEND	0.
#define ANG_MOT_TOR			180.
#define ANG_MOT_BEND		0
#define ANGL_MOT_BEND		15.
#define ANGL_MOT_90			10.
#define ANGL_MOT_TOR		10.
#define DIA_CYL_MOT			10.0e-9
// Mechanical stiffness for motors
#define STF_MOTBACK_SPR		1.691e-2
#define STF_MOTBACK_BEND	5.072e-18
// Geometric parameters for motors
#define L_MOTBACK_DIST		42.0e-9	
#define L_MOTBACK_CEN_DIST	42.0e-9	
#define ANG_MOT_BEND		0
// Boundaries
#define STF_REP_BND			1e-3
// Degree of coarse-graining using cylindrical segments and 
// the corresponding length of one actin segments
#define L_SCALE_IN_M			(nActPerSeg * DIA_CYL_ACT)
#define L_SCALE_IN_UM			(nActPerSeg * DIA_CYL_ACT_UM)
#define L_SCALE_IN_NM			(nActPerSeg * DIA_CYL_ACT_NM)
/*----------------------------------------------------------------------------*/

/******************************************************************************/
/***************************** Definition of Tools ****************************/
/******************************************************************************/
// Alternative forms of arrays
#define P2(a,b,c)			a[(b)*NDIM+(c)]
#define P2A(a,b,c,d)		a[(b)*(d)+(c)]
// Frequently used expressions 
#define BREAK(a)			if (a) { break; }
#define CONT(a)             if (a) { continue; }
#define ISACTM(a)		    (P2A(act.ch,(a),0,nChAc) < 0  				\
		&& P2A(act.ch,(a),1,nChAc) < 0)
#define ISABPM(a)		    (P2A(abp.ch,(a),0,nChAb) < 0 				\
		&& P2A(abp.ch,(a),1,nChAb) < 0 && (K_ABP(a) != 2 				\
		|| motSA.gTgl == 0 || (P2A(abp.ch,(a),3,nChAb) < 0 				\
		&& P2A(abp.ch,(a),4,nChAb) < 0 && K_ABP(a) == 2 				\
		&& motSA.gTgl != 0)))
#define ISABPIM(a)		    (P2A(abp.ch,(a),0,nChAb) < 0 				\
		&& P2A(abp.ch,(a),1,nChAb) < 0 && ((K_ABP(a) != 2 				\
		&& gTglImpAcpM != 0) || (K_ABP(a) == 2 && gTglImpMotM != 0 		\
		&& motSA.gTgl == 0) || (P2A(abp.ch,(a),3,nChAb) < 0 			\
		&& P2A(abp.ch,(a),4,nChAb) < 0 && K_ABP(a) == 2 				\
		&& motSA.gTgl != 0)))
#define ISMTF(a)			(motSA.gTgl != 0 && (a) == 2)
#define K_ABP(a)			abp.ch[(a)*nChAb+2]
#define SPRING(a,b,c)		REVSIGN(a) * ((b) - (c))
#define CYL_DRAG(a,b)		(3.*PI*VISCOSITY*(a)*(3.+2.*(b)/(a))/5.)
// Frequently used for-loops
#define FOR_ACT(a)			for((a) = 0; (a) < nAct; (a)++)
#define FOR_ABP(a)			for((a) = 0; (a) < nAbp; (a)++)
#define FOR_NDIM(a)			for((a) = 0; (a) < NDIM; (a)++)
#define FOR_ACTME(a)		for((a) = 0; (a) < nActMe; (a)++)
#define FOR_ABPME(a)		for((a) = 0; (a) < nAbpMe; (a)++)
#define FOR_ACTCP(a)		for((a) = 0; (a) < nActCp; (a)++)
#define FOR_ABPCP(a)		for((a) = 0; (a) < nAbpCp; (a)++)
#define FOR_ACTMECP(a)		for((a) = 0; (a) < nActMe + nActCp; (a)++)
#define FOR_ABPMECP(a)		for((a) = 0; (a) < nAbpMe + nAbpCp; (a)++)
#define FOR_ACTC(a)			for((a) = 0; (a) < nActC; (a)++)
#define FOR_ABPC(a)			for((a) = 0; (a) < nAbpC; (a)++)
// Algebra
#define AVG2(a,b)			(0.5*((a)+(b)))
#define REVSIGN(a)			(-1.*(a))
#define INV(a)				(1./(a))
#define SQR(a)				((a)*(a))
#define CUBE(a)				((a)*(a)*(a))
// Alternative forms of memory allocation
#define MALLOC(a,b,c)		(a)=(b *)malloc(sizeof(b)*(c))
#define MALLOC2(a,b,c)		(a)=(b **)malloc(sizeof(b *)*(c))
#define MALLOC3(a,b,c)		(a)=(b ***)malloc(sizeof(b **)*(c))
#define MALLOC4(a,b,c)		(a)=(b ****)malloc(sizeof(b ***)*(c))
// Packing and unpacking data in meassages
#define MPI_PACK_DBL(a, b, c)		MPI_Pack((a), (b), MPI_DOUBLE,		\
	bufSendMsg[(c)], sizeBufMsg, &posi, MPI_COMM_WORLD);
#define MPI_PACK_INT(a, b, c)			MPI_Pack((a), (b), MPI_INT,		\
	bufSendMsg[(c)], sizeBufMsg, &posi, MPI_COMM_WORLD);
#define MPI_UNPACK_DBL(a, b, c)		MPI_Unpack(bufRecvMsg[(c)],			\
	sizeBufMsg,&posi,(a),(b), MPI_DOUBLE, MPI_COMM_WORLD);
#define MPI_UNPACK_INT(a, b, c)			MPI_Unpack(bufRecvMsg[(c)],		\
	sizeBufMsg,&posi,(a),(b), MPI_INT, MPI_COMM_WORLD);
// Convert between indices and numbers
#define V2IND_ASSIGN_INT(a,b,c,d)										\
	(b)=(a)/(d), (c)=(a)%(d)
#define V3IND_ASSIGN_DBL(a,b,c,d)										\
	(d)[0]=(double)(((a)/(b)[2])/(b)[1])*(c),							\
	(d)[1]=(double)(((a)/(b)[2])%(b)[1])*(c),							\
	(d)[2]=(double)((a)%(b)[2])*(c)							
#define V3IND_ASSIGN_INT(a,b,c,d)										\
	(d)[0]=(((a)/(b)[2])/(b)[1])*(c),									\
	(d)[1]=(((a)/(b)[2])%(b)[1])*(c),									\
	(d)[2]=((a)%(b)[2])*(c)							
#define V3IND_BACK_INT(a,b,c)											\
	(a)=((b)[0]*(c)[1]+(b)[1])*(c)[2]+(b)[2]
#define V3IND_ASSIGN_CONST_INT(a,b,c)									\
	(c)[0]=(((a)/b)/b),													\
	(c)[1]=(((a)/b)%b),													\
	(c)[2]=((a)%b)							
#define V3IND_BACK_CONST_INT(a,b,c)		(a)=((b)[0]*c+(b)[1])*c+(b)[2]
// Unit conversion between conventional units and simulation units
#define T_SEC2TS(a)			((int)((a) / dtReal))
#define T_S2SEC(a)			((a)*(CYL_DRAG(DIA_CYL_ACT, DIA_CYL_ACT		\
	* nActPerSeg)*SQR(L_SCALE_IN_M)/KT_IN_J))
#define T_SEC2S(a)			((a)/(CYL_DRAG(DIA_CYL_ACT, DIA_CYL_ACT 	\
	* nActPerSeg)*SQR(L_SCALE_IN_M)/KT_IN_J))
#define F_PN2S(a)			((a)/KT_IN_J*L_SCALE_IN_M/1.0e12)
#define F_S2PN(a)			((a)*KT_IN_J/L_SCALE_IN_M*1.0e12)
#define F_N2S(a)			((a)/KT_IN_J*L_SCALE_IN_M)
#define F_S2N(a)			((a)*KT_IN_J/L_SCALE_IN_M)
#define F_N2S(a)			((a)/KT_IN_J*L_SCALE_IN_M)
#define E_S2J(a)			((a)*KT_IN_J)
#define L_S2M(a)			((a)*L_SCALE_IN_M)
#define L_S2UM(a)			((a)*L_SCALE_IN_UM)
#define L_S2NM(a)			((a)*L_SCALE_IN_NM)
#define L_M2S(a)			((a)/L_SCALE_IN_M)
#define L_UM2S(a)			((a)/L_SCALE_IN_UM)
#define L_NM2S(a)			((a)/L_SCALE_IN_NM)
#define KS_S2NPM(a)			((a)*KT_IN_J/SQR(L_SCALE_IN_M))
#define KS_NPM2S(a)			((a)/KT_IN_J*SQR(L_SCALE_IN_M))
#define KB_S2NM(a)			((a)*KT_IN_J)
#define KB_NM2S(a)			((a)/KT_IN_J)
#define PA2S(a)				((a)/KT_IN_J*CUBE(L_SCALE_IN_M))
#define S2PA(a)				((a)*KT_IN_J/CUBE(L_SCALE_IN_M))
#define NPV2S(a)			((a)/KT_IN_J*pow(L_SCALE_IN_M, 4.))
#define S2NPV(a)			((a)*KT_IN_J/pow(L_SCALE_IN_M, 4.))
// Unit conversion of angles
#define DEG2RAD(a)			((a)*PI/180.)
#define RAD2DEG(a)			((a)/PI*180.)
// Probabilities and rates
#define K2P(a) 				(1.-exp(-1.*(a)*dtReal))
#define P2K(a) 				(log(1.-(a))/(-1.*dtReal))
#define FAC_P(a) 			(-1.*(a)/N_AVO/CUBE(L_SCALE_IN_M)*1.0e3*dtReal)
/*---------------------------- Vector calculation ----------------------------*/

// for 3-D vector
#define V3REVSIGN(a)		(a)[0]=-1.*(a)[0], (a)[1]=-1.*(a)[1],		\
	(a)[2]=-1.*(a)[2]
#define V3DOT(a,b)			((a)[0]*(b)[0]+(a)[1]*(b)[1]+(a)[2]*(b)[2])
#define V3CROSS(a,b,c)		(a)[0]=(b)[1]*(c)[2]-(b)[2]*(c)[1], 		\
	(a)[1]=(b)[2]*(c)[0]-(b)[0]*(c)[2], (a)[2]=(b)[0]*(c)[1]-(b)[1]*(c)[0]
#define V3ADD(a,b,c)													\
	(a)[0]=(b)[0]+(c)[0],(a)[1]=(b)[1]+(c)[1],(a)[2]=(b)[2]+(c)[2]
#define V3SUB(a,b,c)													\
	(a)[0]=(b)[0]-(c)[0],(a)[1]=(b)[1]-(c)[1],(a)[2]=(b)[2]-(c)[2]
#define V3LEN_SQ(a)			V3DOT(a,a)
#define V3LEN(a)			sqrt(V3DOT(a,a))
#define V3SET(a,b,c,d)		(a)[0]=(b), (a)[1]=(c), (a)[2]=(d)
#define V3SET_ALL(a,b)		(a)[0]=(b), (a)[1]=(b), (a)[2]=(b)
#define V3ZERO(a)			V3SET_ALL(a,0)
#define V3MUL(a,b,c)		(a)[0]=(b)[0]*(c)[0], (a)[1]=(b)[1]*(c)[1], \
	(a)[2]=(b)[2]*(c)[2]
#define V3DIV(a,b,c)		(a)[0]=(b)[0]/(c)[0], (a)[1]=(b)[1]/(c)[1], \
	(a)[2]=(b)[2]/(c)[2]
#define V3DIV_INT(a,b,c)	(a)[0]=(int)((b)[0]/(c)[0]), 				\
							(a)[1]=(int)((b)[1]/(c)[1]),				\
							(a)[2]=(int)((b)[2]/(c)[2])
#define V3SCALE(a,b)		(a)[0]*=(b), (a)[1]*=(b), (a)[2]*=(b);
#define VV3ADD(a,b)			V3ADD(a,a,b)
#define VV3SUB(a,b)			V3SUB(a,a,b)
#define VV3DIV(a,b)			V3DIV(a,a,b)
#define V3PROD(a)			((a)[0]*(a)[1]*(a)[2])
#define V3SUM(a)			((a)[0]+(a)[1]+(a)[2])
#define VS3ADD(a,b,c,d)													\
	(a)[0]=(b)[0]+(c)[0]*(d),(a)[1]=(b)[1]+(c)[1]*(d),					\
	(a)[2]=(b)[2]+(c)[2]*(d)
#define VSV3ADD(a,b,c,d)												\
	(a)[0]=(b)[0]+(c)[0]*(d)[0],(a)[1]=(b)[1]+(c)[1]*(d)[1],			\
	(a)[2]=(b)[2]+(c)[2]*(d)[2]
#define V3AVG(a,b,c)													\
	V3ADD(a,b,c), V3SCALE(a,0.5)
#define VVS3ADD(a,b,c)		VS3ADD(a,a,b,c)
#define VSS3ADD(a,b,c,d,e)												\
	(a)[0]=(b)[0]*(d)+(c)[0]*(e),(a)[1]=(b)[1]*(d)+(c)[1]*(e),			\
	(a)[2]=(b)[2]*(d)+(c)[2]*(e)
#define VS3SUB(a,b,c,d)													\
	(a)[0]=(b)[0]-(c)[0]*(d),(a)[1]=(b)[1]-(c)[1]*(d),					\
	(a)[2]=(b)[2]-(c)[2]*(d)
#define VVS3SUB(a,b,c)		VS3SUB(a,a,b,c)
#define VVSS3SUB(a,b)		(a)[0]-=(b), (a)[1]-=(b), (a)[2]-=(b)
#define VVSS3ADD(a,b)		(a)[0]+=(b), (a)[1]+=(b), (a)[2]+=(b)
#define V3COPY(a,b)			(a)[0]=(b)[0],(a)[1]=(b)[1],(a)[2]=(b)[2]
#define VS3COPY(a,b,c)													\
	(a)[0]=(b)[0]*(c),(a)[1]=(b)[1]*(c),(a)[2]=(b)[2]*(c)
#define V3COS(a,b)			(V3DOT((a),(b))/(V3LEN(a)*V3LEN(b)))
#define V3ANG(a,b)			Acos(V3COS(a,b))
// for 2-D vector
#define V2SET(a,b,c)		(a)[0]=(b), (a)[1]=(c)
#define V2SET_ALL(a,b)		(a)[0]=(b), (a)[1]=(b)
#define V2ZERO(a)			V2SET_ALL(a,0)
#define V2COPY(a,b)			(a)[0]=(b)[0],(a)[1]=(b)[1]
// for 4-D vector
#define V4SET(a,b,c,d,e)	(a)[0]=(b), (a)[1]=(c), (a)[2]=(d), (a)[3]=(e)
#define V4SET_ALL(a,b)		V4SET(a,b,b,b,b)
#define V4COPY(a,b)			(a)[0]=(b)[0],(a)[1]=(b)[1],(a)[2]=(b)[2],  \
	(a)[3]=(b)[3]
// for 5-D vector
#define V5SET(a,b,c,d,e,f)	(a)[0]=(b), (a)[1]=(c), (a)[2]=(d),			\
	(a)[3]=(e), (a)[4]=(f)
#define V5SET_ALL(a,b)		V4SET(a,b,b,b,b,b)
#define V5COPY(a,b)			(a)[0]=(b)[0],(a)[1]=(b)[1],(a)[2]=(b)[2],  \
	(a)[3]=(b)[3],(a)[4]=(b)[4]
#define V6SET(a,b,c,d,e,f,g)	(a)[0]=(b), (a)[1]=(c), (a)[2]=(d),			\
	(a)[3]=(e), (a)[4]=(f), (a)[5]=(g)
#define V6COPY(a,b)				(a)[0]=(b)[0],(a)[1]=(b)[1],(a)[2]=(b)[2],  \
	(a)[3]=(b)[3],(a)[4]=(b)[4], (a)[5]=(b)[5] 
#define V7COPY(a,b)				(a)[0]=(b)[0],(a)[1]=(b)[1],(a)[2]=(b)[2],  \
	(a)[3]=(b)[3],(a)[4]=(b)[4], (a)[5]=(b)[5], (a)[6]=(b)[6]

/******************************************************************************/
/***************************** Function Prototypes ****************************/
/******************************************************************************/
/*-------------------------------- process.c ---------------------------------*/
// Main process 
void MainProcess(void);
// Initialize and excute single steps & post-process after single steps
void SingleProcess(void), TaskAfterSingleProcess(void);
void InitSingleProcess(void);
void InitSingleProcessSubroutine(double *arr, int cnt) ;
// Switch toggles 
void InitToggleParameters(int mode);
// Procedures necessary without pre-assembled network data
void PrepareStateWoNetworkData(void);
/*---------------------------------- init.c ----------------------------------*/
// Initialize the number of cells based on the given number of cores
void InitCellNumber(void);
// Initialize process
void InitRun(void), InitOutput(void);
// Define arrays and assign values for them
void AllocPreArrays(void), AllocArrays(void), AssignArrayValues(void);
// Assign initial values of variables
void AssignInitValues(void), CheckPeriodToggleParameter(FuncCont *);
// Initialize record files and seeds for generating random numbers
void InitRecFiles(void), InitRandSeed(void);
// Check whether 'condition' and 'parallel' exist
void InitFileCheck(void);
// Initialize the unbinding and walking rates based on the number of heads
void InitMotorWalkUnbindRates(int);
/*-------------------------------- calForce.c --------------------------------*/
// Brownian forces (thermal fluctuation)
void CalcBrownForces(void);
void CalcBrownForcesSubroutine(double *, double, int, int, double *);
// Repulsive forces between cylindrical segments
void CalcRepulsiveForces(void);
void CalcRepulsiveForcesSubroutine(int *, int *, double [][NDIM],
        double *[4], double *, int);
double CalcRepulsiveForcesSubSubroutine(double [][NDIM], double *, 
		double *, double, int);
void CalcRepulsiveForcesSubSubSubroutine(double [][NDIM], int, int);
// Spring forces
void CalcSpringForces(void), CalcSpringForcesSubroutine(int, int, int, int *);
// Bending forces and necessary subroutines
void CalcFilaBendForces(void), CalcFilaBendForcesSubroutine(int, int, int);
void CalcAbpBendForces(void), CalcAbpBendForcesSubroutine1(int, int);
void CalcAbpBendForcesSubroutine2(int, int *);
// For motor thick filaments
void CalcMotorBackboneForces(void);
// Tools for force calculation
void CalcCosine(double*,double *,double *,double *,double *,double *,double *);
double CalcBendForceSubroutine(double *, double *, double, double, 
		double *, double *);
double CalcBendForce(double *, double *, double *, double *, 
		double *, double, double, double *, double *);
void AddSpringForce(double, double, double *, double *, double *);
/*---------------------------------- update.c --------------------------------*/
// Update locations of particles reflecting forces
void UpdateNewLocation(void);
// Update neighboring list 
void UpdatePrevLocationForNL(void), MeasureDisplacementForNL(void);
void UpdateNeighborListSubroutine(double [][NDIM], double *, int);
int UpdateNeighborListSubroutine2(int *, int *);
int UpdateNeighborListSubroutine3(int, int *, int *);
void UpdateNeighborList(void);
void DeleteActinSegmentInNeighborList(int *);
void DeleteElementInNeighborList(int, int);
void InsertElementInNeighborList(int, int, int);
// Dynamic behaviors of actins
void UpdateActinNucleation(void), UpdateActinAssembly(void);
void UpdateActinDisassembly(void), UpdateActinDisassemblySubroutine(int, int);
void UpdateActinDisassemblySubroutine2(int, int);
void UpdateNoActinDynamicsList(void);
int CheckActinAvailability(int, int);
// Dynamic behaviors (unbinding/binding/walking) of ACPs and motors
void UpdateAbpBinding(void), UpdateAbpBindingSubroutine(int *, int);
void UpdateActiveAbpUnbinding(void);
void UpdateActiveAbpUnbindingSubroutine(int, int);
void UpdateInactAbpUnbinding(void), UpdateInactAbpUnbindingSubroutine(int ,int);
void UpdateAbpMonomerBinding(void), UpdateMotorWalking(void);
void UpdateMotorAssembly(void), UpdateMotorTurnover(void);
void UpdateNoAbpUnbindList(void);
double UpdateMotorWalkingSubroutine(int, int);
double CalcUnbindingRate(int, int, int, int);
// Dynamics of rho activation
void UpdateRhoActivation(void);
int CheckRhoActiveRegion(double *);
// Updating chains and information
void UpdateChainList(void);
/*------------------------------- dataManage.c -------------------------------*/
// Load initial parameters
void LoadInitParameter(void);
void LoadInitParameterSubroutine(FILE *, const char *, int *);
void LoadInitParameterSubroutine2(char *, int, int *, int *);
void CheckAnswerYesNo(char *, const char *, int *);
// Eliminate inactive or free components
void AddFreeAbp(void), AddFreeActin(void);
void AddFreeActinAbpSubroutine(int, int *, int *);
// Pack information
int BinaryPackActinChainArray(int), BinaryPackAbpChainArray(int);
/*--------------------------------- record.c ---------------------------------*/
// Progress 
void RecordProgressSubroutine(int *, FILE *);
void RecordProgressSubroutine2(int *, int *, int, FILE *, int);
void RecordProgressSubroutine3(int *, int *, FILE *);
void RecordProgress(void);
// Record initial parameters
void RecordInitParameterSubroutine(const char *, int, int, FILE *); 
void RecordInitParameterSubroutine2(const char *, FuncCont *, FILE *);
void RecordInitParameter(void);
// Network configuration
void RecordConfig(char *, int), RecordConfigVmd(int);
int RecordConfigVmdSubroutine(double *, double *, int *, int *);
// Record information related to the unbinding event of ABPs
void RecordAbpUnbindEvent(int, int, int), RecordAbpBindEvent(int, int, int);
void RecordAbpTurnover(int, int, int, int);
// Morphological properties
void RecordFilamentLength(int);
// Energy 
void RecordMechEnergy(int);
// Accumulated information
void RecordAccuLengthForces(void), ResetAccuLengthForces(void);
void RecordAccuLongSpringForces(int);
// Information in unit of individual elements, filament segments, filaments
void RecordIndvSegFilaInformation(int);
/*--------------------------------- gatPrint.c -------------------------------*/
// Tools for recording
void Printf(const char *, ...), Printf0(const char *, ...);
void Fprintf1dFillerInt(FILE *, int, int, int);
void Fprintf1dArrayInt(FILE *, int *, int, int);
void Fprintf1dArrayIntWFil(FILE *, int *, int, int, int);
void Fprintf1dArrayIntWoRet(FILE *, int *, int); 
void Fprintf1dFillerDouble(FILE *, double, int, int);
void Fprintf1dArrayDouble(FILE *, double *, int, int);
void Fprintf1dArrayDoubleWFil(FILE *, double *, int, int, int); 
void Fprintf1dArrayDoubleWoRet(FILE *, double *, int); 
void Fprintf2dArrayIntDouble(FILE *, double *, int, int, int *);
void Fprintf2dArrayInt(FILE *, int *, int, int, int, int);
void Fprintf2dArrayIntWFil(FILE *, int *, int, int, int, int, int); 
void Fprintf2dArrayDouble(FILE *, double *, int, int, int, int);
void Fprintf2dArrayDoubleWFil(FILE *, double *, int, int, int, int, int); 
void RecordChainArrayInt(FILE *, const char *, int *, int, int, int); 
void Gather1dArrayIntWoInd(int, int *, int *, int *);
void RecordGather1dArrayIntWoInd(int, int *, int *, FILE *);
void RecordGather1dArrayIntWoIndWoCnt(int, int *, FILE *, const char *);
void Gather2dArrayInt(int, int *, int, int *, int *, int *);
void Gather2dArrayIntWoIndWoSort(int, int *, int *, int *, int);
void RecordGather2dArrayInt(int, int *, int, int *, int *, FILE *, int);
void RecordGather2dArrayIntWoCnt(int, int, int *, int *, FILE *, int);
void RecordChainArrayDouble(FILE *, const char *, double *, int, int, int);
void Gather1dArrayDoubleWoInd(int, int *, double *, double *);
void RecordGather1dArrayDoubleWoInd(int, int *, double *, FILE *);
void RecordGather1dArrayDoubleWoIndWoCnt(int, double *, FILE *, 
		const char *);
void Gather2dArrayDouble(int, int *, int, int *, double *, double *);
void Gather2dArrayDoubleWoIndWoSort(int, int *, double *, double *, int);
void RecordGather2dArrayDouble(int, int *, int, int *, double *, FILE *, int);
void RecordGather2dArrayDoubleWoCnt(int, int, int *, double *, FILE *, int);
void GatherActChainPosition(int *, int *, double *);
void GatherAbpChainPosition(int *, int *, double *);
/*---------------------------------- error.c ---------------------------------*/
// Print errors if any 
void RecordError(int);
void RecordErrorArrayInt(FILE *, const char *, int, int *, int);
void RecordErrorArrayDouble(FILE *, const char *, int, double *, int);
void RecordErrorSubroutine1(FILE *, int), RecordErrorSubroutine2(FILE *, int);
void RecordErrorSpringForce(int, int, double, double, int);
void RecordErrorRepForce(int *, int *, double (*)[NDIM], double, 
		double, double *, int);
void RecordErrorBendingForce(int, int, int, double, int);
void RecordErrorTotalForceSubroutine1(FILE *, int);
void RecordErrorTotalForceSubroutine2(FILE *, int);
void RecordErrorTotalForce(int, int), RecordErrorElement(int, int);
int CheckLargeForce(double, int);
void CheckLargeTotalForceAll(int), CheckNanForce(int);
/*----------------------------------- rng.c ----------------------------------*/
// Initialize the generation of rundom numbers
void init_genrand(unsigned long);
// Generate a Gaussian random number
void genrand_gauss(double *, double *);
// Generate a double (uniform) random number
double genrand_real3(void);
/*--------------------------------- boundary.c -------------------------------*/
// Handle boundary conditions
void CheckCrossBound(int *, double *, double *);
void ApplyBoundCondVector(double *, int, int);
void ApplyBoundCondVecDiff(double *), ApplyBoundCondAll(void);
void ConvertRectDomainVector(double *, int);
// Interactions between boundaries and others
int CheckParticleInDomain(double *);
/*--------------------------------- paraProc.c -------------------------------*/
// Handle the list of long chains between actins and ABPs
void UpdateLongChainNormal(void);
void UpdateLongChainNormalSubroutine(int, int, double);
void UpdateLongChainNormalSubSubroutine1(double *, ListInt *, double);
void DeleteLongChain(int, int);
int InsertLongChain(int, int, double);
// Move particles between subdomains
void MoveParticles(void), MoveParticlesSubroutine2(int, int);
void MoveParticlesSubroutine1(ListInt *, double *, int);
// Copy particles between subdomains
void CopyParticles(void), CopyParticlesSubSubroutine(int, int *, int);
void CopyParticlesNormalSubroutine(int, int *, int *, int);
// Process the messages regarding the dynamic events of actins and ABPs
// which will be transferred between subdomains
void UpdateActinAbpDynamicsEvents(int); 
void UpdateActinAbpDynamicsEventsSubroutine(ListInt *, int);
void UpdateActinAbpDynamicsEventsSubSubroutine(ListInt *, int, int, int, 
		ListInt *, int *);
void UpdateAbpUnbRebLists(int, int, int, int);
// Manage actin and ABP monomers
void UpdateActinAbpMonomerListSubroutine(int, int, int);
void UpdateActinAbpMonomerListSubroutine2(int, int, int);
// Collect arrays from adjacent subdomains
void CollectArrayIntFromAdjacentSubdomain(ListInt *, int);
void CollectArrayDblFromSubdomainList(double *, double *, int, ListInt *, int);

/*---------------------------------- tools.c ---------------------------------*/
// Find maximum and minimum values of values in an array
int FindMaxMin1dArrayInt(int *, int, int);
double FindMaxMin1dArrayDbl(double *, int, int);
// Set values for an array
void SetAllValue1dArrayInt(int *, int, int);
void SetAllValue1dArrayDouble(double *, int, double);
// Copy values in an array to other array
void Copy1dArrayInt(int *, int *, int);
void Copy1dArrayDouble(double *, double *, int);
// Find, insert, or delete elements in an array
int FindElementArray(int *, int, int, int, int);
int Find2ElementArray(int *, int, int, int, int, int);
int Find3ElementArray(int *, int, int, int, int, int, int);
void InsertElement1dArrayWoChk(int *, int *, int);
int InsertElement1dArrayWChk(int *, int *, int);
void InsertElementArrayByIndex(int *, int *, int *, int, int);
void DeleteElement1dArray(int *, int *, int);
void DeleteElementArrayByIndex(int *, int *, int, int);
int Delete2ElementArray(int *, int *, int, int, int, int);
// Sum or average values in array
double AvgArrDbl(double *, int), AvgArrInt(int *, int);
double SumArrDbl(double *, int);
int SumArrInt(int *, int);
// Check the size of array
void CheckArraySize(ListInt *, int *, int, int);
void Check2dArraySize(ListInt2 *, int *, int, int);
// Vector and geometry calculation (general)
double CalcSegPntDist(double [][NDIM], double *, double *, double *);
double CalcSegSegDist(double [][NDIM], double [][NDIM], double *, 
		double *, int);
void NormVec(double *), CalcVec(double *, double *, double *);
double CalcUnitVec(double *, double *, double *);
double CalcDist(double *, double *, int);
double CalcVecDist(double *, double *, double *, int);
// Vector calculation for Specific for elements
double CalcVecActinAbp(double *, int, int, int);
void CalcUnitVecActinAbp(double *, int, int, int);
double CalcDistActinAbp(int, int, int);
double CalcVecDistActinAbp(double *, int, int, int);
void CalcPosOnActSeg(double *, double *, double *, double);
void CalcPosOnActSegSide(double *, double *, double *, int);
void OffsetSegEndPnt(double *, double *, double *);
void CalcInactAbpPosition(double *, double *, double *, double *, int, int);
void CalcAbpArmEndPos(double *, int, int);
// Find information related to actin and ABP
int FindAbpActinChain(int, int, int), HowManyAbpActinChain(int, int);
// Generate random information
void GenRandDirecVec(double *), GenRandPosSubdom(double *);
int GenRandIntIndex(int);
// Calculate the rank or iCell of particles depending on location
void CalcIndMolecule(double *, int *);
int CalcRankMolecule(double *), CalcRankIndMolecule(double *, int *);
// Adjustment for rates of dynamic behaviors of actins and ABPs
double AdjustDynamicsRate(double);
// Matrix calculation
void MultiplyMatrixIntSerial(int *, int *, int *, int, int, int);
void MultiplySqMatrixIntParallel(int *, int *, int *, int);
// Misc.
char *GenFileName(const char *);
void SwapInt(int *, int *), SwapDbl(double *, double *);
int TrimIntVal(int, int, int), SetKind(int, int, int), SignDbl(double);
double TrimDblVal(double, double, double), Acos(double), Sqrt(double);
int CompInt(const void *, const void *);
int CompDbl(const void *, const void *);
char* IntToStr(int, int);

/******************************************************************************/
/****************************** Global Variables ******************************//******************************************************************************/
 
// Variables related to time steps
extern long long currTimeStep;
extern double dt, dtReal;
extern time_t initTime;
extern Duration netForm, rheo;
// Given concentration
extern int nAbpGoal, nAbpDet[NK_ABP], nAbpMall[NK_ABP], nAbpGoalDet[NK_ABP];
extern double cAct, RAbp[NK_ABP];
// Number of particles
extern int nAct, nAbp, nMot; 
extern int nActMe, nActCp, nActMeCp, nActC;
extern int nAbpMe, nAbpCp, nAbpMeCp, nAbpC;
extern int nAcpInaMe, nMotInaMe, nAcpMme, nMotMme;
extern int nActGoal, nActMall, nActFilaMe;
extern int nActMin, nAbpMin;
// Chain, position, forces, and lists of particles
extern ListInt actM, acpM, motM, iFilaP;
extern MoleAct act;
extern MoleAbp abp;
extern ActF actF;
extern AbpF abpF;
extern MotSelfAss motSA;
extern int *chAct, *chAbp, tglNeiAbpSC, tglNeiAbpDC;
extern int nChAc, nChAcX, nChAcY, nChAb, nActPerSeg;
extern double *rAct, *rAbp;
// Update neighboring list
extern ListInt neigh;
extern Cell cell;
extern double dispHeuSq, maxNeiLenSq;
// Arrays related to lists
extern ListInt sendAbpDyn, noAbpDyn, sendActDyn, noActDyn;
extern int *fixAct, *abpMotId;
// Domain (width, periodic boundary conditions, or repulsive force)
extern int pbc[NDIM], neiPbc[NDIM], dir2D;
extern double dimDom[NDIM], dimDomH[NDIM], minDimDomC;
extern Boundary bnd;
// For boundaries
extern double *rGridInit;

/*-------------------------- For parallel processing -------------------------*/
extern int nCpu;
extern int *iAct, *iAbp;
// Ranks of CPUs
extern int rank, *adjRank, *iRank, *cntAdjRank;
// For boundaries and indicies of subdomains
extern int nCell[NDIM], iCell[NDIM], nGrid[NDIM];
extern double edge[NDIM*2], **rGrid, neiEdge;
// Used in CopyParticles and MoveParticles
extern int *cntCpPar, *cntMvPar, modeActCh;
extern ListInt *cpPar, *mvPar, insNeiPar;
// Varilables related to messages
extern int sizeBufMsg, *mpiTestSendFlag, *mpiTestRecvFlag;
extern char **bufSendMsg, **bufRecvMsg;
extern MPI_Status status;
extern MPI_Request *sReq, *rReq;
// Variables related to longCh
extern ListInt longCh, longChExtMsg, longChIntMsg;
extern double *longChDist, maxDisp, maxActCh;
/*----------------------------------------------------------------------------*/

/*---------------- Dynamic behaviors of actin, ACP, and motor ----------------*/
// Dynamics of actin
extern AssDis actAss, actDis;
extern Nucle actNuc;
extern int tglActMoDyn, tglActFormDyn, gTglActDynNF;
extern int gTglActTherm, gTglAcpTherm, gTglMotTherm;
extern int durNoActDyn;
// Dynamics of ACPs and motors
extern Reb acpReb, motReb;
extern Unb acpUnb; 
extern MotUnbWalk motUnb, motWalk; 
extern InaUnbMbind acpInaUnb, motInaUnb, acpMoBind, motMoBind;
extern int tglAbpAcInaDyn, tglAbpInaMoDyn, tglAcpAcInaDyn, tglMotAcInaDyn;
extern int gTglAcpDynNF, gTglMotUnbRebNF, gTglMotWalkNF, gTglMotWalkSld;
extern int durNoAcpUnbReb, durNoMotUnbReb, durNoMotWalk;
extern int gTglImpAcpM, gTglImpMotM;
extern ListDbl unbLog, bindLog, toLog;
extern MotMechChem motMC;
extern int cntNucAss;
// Rho activity
extern Rho rho;
/*----------------------------------------------------------------------------*/
 
/*---------------------------- For data recording ----------------------------*/
extern FuncCont recProg, recConf;
// For general records
extern char fnOut[80], dataFold[80];
// Record the accumulated chain lengths and forces
extern RecLenForce recAct, recAbp;
// Record configuration for VMD
extern RecConfVmd recConfVmd;
// Record the turnover of ACPs or motors 
extern int tglRecAbpTurn;
extern double *abpTurn;
// Record (longitudinal) forces of ACPs or motors
extern FuncCont recLongF;
extern double *recLongSprFabp, *recInstSprFabp;
// Record etc
extern FuncCont confVmdInfo, recFilaL, recE, recInfo;
extern FuncCont recAbpUnb, recAbpBind, recAbpTurn;
/*----------------------------------------------------------------------------*/

// Check errors
extern int stopSig;
extern double magUnstF;
// Misc.
extern int seed, *allIntL;
extern double *arrAcos, *allDblL;

