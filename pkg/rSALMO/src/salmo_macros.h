/* file: salmo_macos.h

   preprocessor macros to improve code readability
   2008-01-01
   
   ToDo:
   - replace ixxx variables with constants in *final* #defines
   - remove remaining cB[....] calls from SalmoKern.cpp

*/

#define inumberOfInputs      0
#define inumberOfOutputs     1
#define inumberOfStates      2
#define inumberOfParameters  3
#define inumberOfPhyto       4
#define inumberOfLayers      5
#define inumberOfTributaries 6

// Order of states is still hard coded
// N 0
// P 1
// X 2 .. 2 + (nx-1)
// Z 2 + nx
// D 3 + nx
// O 4 + nx
// G 5 .. 5 + (nx-1)


#define ixin       0
#define iwx        1
#define irx        2
#define irxt       3
#define iphotx     4
#define iphoxt     5
#define iphoxl     6
#define iphoxn     7
#define iphoxp     8
#define iphoxns    9
#define ikzi      10
#define ihwgdb    11
#define ihwgdi    12
#define ihwg      13
#define igi       14
#define ipf       15
#define ixwa      16
#define ibx       17
#define ixsed     18
#define ixgraz    19
#define ixim      20
#define iminer    21
#define cB_rows   21  // how many different values are in cB for each j

#define	xin_j	    cB[j+nx*ixin]
#define	wx_j	    cB[j+nx*iwx]
#define	rx_j	    cB[j+nx*irx]
#define	rxt_j	    cB[j+nx*irxt]
#define	photx_j	  cB[j+nx*iphotx]
#define	phoxt_j	  cB[j+nx*iphoxt]
#define	phoxl_j	  cB[j+nx*iphoxl]
#define	phoxn_j	  cB[j+nx*iphoxn]
#define	phoxp_j	  cB[j+nx*iphoxp]
#define	phoxns_j	cB[j+nx*iphoxns]
#define	kzi_j	    cB[j+nx*ikzi]
#define	hwgdb_j	  cB[j+nx*ihwgdb]
#define	hwgdi_j	  cB[j+nx*ihwgdi]
#define	hwg_j	    cB[j+nx*ihwg]
#define	gi_j	    cB[j+nx*igi]
#define	pf_j	    cB[j+nx*ipf]
#define	xwa_j	    cB[j+nx*ixwa]
#define	bx_j	    cB[j+nx*ibx]
#define	xsed_j	  cB[j+nx*ixsed]
#define	xgraz_j	  cB[j+nx*ixgraz]
#define	xim_j	    cB[j+nx*ixim]
#define	miner_j	  cB[j+nx*iminer]

#define	EPSX_j	   pC[j+nx*pEPSX]
#define	BETA_j	   pC[j+nx*pBETA]
#define	GAMMAX_j	 pC[j+nx*pGAMMAX]
#define	GAMMIN_j	 pC[j+nx*pGAMMIN]
#define	KI_j	     pC[j+nx*pKI]
#define	KN_j	     pC[j+nx*pKN]
#define	KP_j	     pC[j+nx*pKP]
#define	KPF_j	     pC[j+nx*pKPF]
#define	NFIX_j	   pC[j+nx*pNFIX]
#define	PFC_j	     pC[j+nx*pPFC]
#define	PFX_j	     pC[j+nx*pPFX]
#define	PHOTXMAX_j pC[j+nx*pPHOTXMAX]
#define	PHOTXMIN_j pC[j+nx*pPHOTXMIN]
#define	RXTMIN_j	 pC[j+nx*pRXTMIN]
#define	RXTOPT_j	 pC[j+nx*pRXTOPT]
#define	RXTSLOPE_j pC[j+nx*pRXTSLOPE]
#define	TOPTX_j	   pC[j+nx*pTOPTX]
#define	UXZ_j	     pC[j+nx*pUXZ]
#define	VS_j	     pC[j+nx*pVS]
#define	YX_j	     pC[j+nx*pYX]
#define	ZETA_j	   pC[j+nx*pZETA]



