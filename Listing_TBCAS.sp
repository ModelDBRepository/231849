File: GME & GPE
**********************************************************************
*				by
* 	Giuseppe Massobrio - Sergio Martinoia - Paolo Massobrio                      
**********************************************************************
* GME: 		Gold Mushroom-Shaped MicroElectrode 
*
* GPE: 		Gold planar MicroElectrode	       
*
* NEURON: 	Behavioral model following the Hodgkin-Huxley formalism
*
* Note:		voltage units should be read as multiplied by 10**-3
**********************************************************************
.OPTIONS list node pathnum probe post method=gear numdgt=9 tnom=25.0 ingold=1
.OPTIONS converge=1 gmindc=1.0e-12 acout=0
*---------------------------
.TEMP 25.0
*---------------------------
IstimGME	11	0  PWL (0 0 1us 10 20s 7 20.000001s  0 80s 0)
IstimGPE	111	0  PWL (0 0 4s 10 20s 7 20.000001s  0 80s 0)
*---------------------------
.TRAN	1E-3	10	0	1E-3     UIC  
*---------------------------
*	Extracellular potential from GME/GPE
.PROBE	TRAN	V(0,840) V(0,940)
*---------------------------
.PARAM
*---------------------------------------------------------------------
*				CONSTANTS
*---------------------------------------------------------------------
* k		Boltzmann constant				       	[J/k]
* TREF		temp. at which model params are measured and extracted	[°C]
* temper	current circuit temperature. 				[°C]
* T		actual absolute temperature				[K]
* conv		conversion factor from mole/l to mole/m**3
* q		electronic charge					[Coulomb]
* eps0		dielectric permittivity	of vacuum			[F/m]
* epsw		relative permittivity of the diffuse layer
* epslys	relative permittivity of the insulation coating
* epsglyco	relative permittivity of the glycocalyx EDL
* NAv		Avogadro constant					[1/mole]
* Cbulk		electrolyte concentration				[mole/l]
* ros		electrolyte resistivity					[ohm*m]
* roAu		resistivity of Au					[ohm*m]
* roglyco	resistivity of the glycocalyx EDL			[ohm*m]
* dWaals	Van der Waals gap	 				[m]
*---------------------------------------------------------------------
+	conv='1.0e3'			
+	CTOK='273.15' 		T='(temper+CTOK)'	
+	k='1.380650e-23'	PI='(355/113)'		PI2='2.0*PI'
+	q='1.602176e-19'	ET='q/(k*T)'		
+	NAv='6.0221415e23'	dWaals='0.34n'
+	eps0='8.854187e-12'	epsw='78.5'	
+	epslys='3.0'		epsmed='eps0*epsw'	epsglyco='81'
+	ros='0.7'		roAu='2.27e-8'		roglyco='1.6e6'
+	Cbulk='150.0m'	
**********************************************************************
*				NEURON PARAMETERS
**********************************************************************
* laxont	Neuron axon length				[m]
* daxont	Neuron axon diameter				[m]
* dsomat	Neuron soma diameter				[m]
* Rax		axoplasmatic resistance per unit length		[ohm*cm]
* Ncomp		compartments number
+	Ncomp='5.0'
+	dsomat='100.0u'	daxont='10.0u'	laxont='500.0u'
+	Asomat='(2.0*PI2*(dsomat/2.0)*(dsomat/2.0))'
+	Aaxont='(PI2*(daxont/2.0)*laxont)'
+	ANEURt='(Asomat+Aaxont)'
+	dsoma='(dsomat/Ncomp)'	daxon='daxont'	laxon='(laxont/Ncomp)'
+	Asoma='(2.0*PI2*(dsoma/2.0)*(dsoma/2.0))'
+	Aaxon='(PI2*(daxon/2.0)*laxon)'
+	ANEUR='(Asoma+Aaxon)'
+	Rax='0.1'	
**********************************************************************
*				GME PARAMETERS
**********************************************************************
* phstalk	percentage of engulfed stalk height 		
* hstalk	GME stalk height 				[m]
* htstalk	GME engulfed stalk height 			[m]
* dmstalk	GME stalk diameter 				[m]
* dmhead	GME head diameter 				[m]
* Astalk	GME stalk area					[m**2]
* Ahead		GME half head area				[m**2]
* Arim		GME rim area					[m**2]
* AGME		GME total area					[m**2]
* tinsM		thickness of the insulation coating		[m]			
+	phstalk='1.0'			
+	dmstalk='1.0u'
+	dmhead='2.0u'
+	hstalk='1.0u'
+	htstalk='(hstalk*phstalk)'
+	cirhead='(PI2*(dmhead/2.0))'
+	Abstalk='(PI*(dmstalk/2.0)*(dmstalk/2.0))'
+	Astalk='(PI2*(dmstalk/2.0)*htstalk)'
+	Aheadt='(2.0*PI2*(dmhead/2.0)*(dmhead/2.0))'
+	Ahead='(Aheadt/2.0)'
+	Arim='(PI*((dmhead/2.0)*(dmhead/2.0)-(dmstalk/2.0)*(dmstalk/2.0)))'
+	AGME='(Astalk+Ahead+Arim)'
+	tinsM='1.2'
*---------------------------------------------------------------------
+ 	Rstalk='(roAu*(htstalk/Abstalk))'				
+ 	Rhead='(roAu*(1.0/Cirhead))'	
+ 	RGME='(Rstalk+Rhead)'
*---------------------------------------------------------------------
+	PHIM='(PI/2.0)'
+ 	CshGME='(((PI2*eps0*epslys)/log(tinsM))*(htstalk+(dmhead/2.0)*sin(PHIM)))'
**********************************************************************
*				NUMBER of GMEs 
**********************************************************************
* shxn		minimum distance between two GMEs head centers		[m]
* xGME		extra distance between two GMEs head centers 		[m]
* shn		effective inter-GMEs head distance 			[m]
* lambdaGME	Debye length						[m]
* nGME		number of vertically arranged GMEs in the soma compartment
*-----------------------------
+	lambdaGME='(sqrt((epsmed/(q*ET))*(1.0/(2.0*NAv*Cbulk*conv))))'
+	shxn='(dmhead+dWaals+2.0*lambdaGME)'	
+	xGME='6.0u'
+	shn='(shxn+xGME)'			
+	AmGME='(PI*shn*shn)'
+	AmSOMA='(PI*(dsoma/2.0)*(dsoma/2.0))'
+	nGME='(int(AmSOMA/AmGME))'
**********************************************************************
*				GPE PARAMETERS
**********************************************************************
* dGPE		GPE diameter 						[m]
* tGPE		GPE thickness						[m]
* tinsP		thickness of the insulation coating			[m]
*-----------------------------
+	dGPE='20u'	tGPE='200n'	AGPE='PI*(dGPE/2.0)*(dGPE/2.0)'
+	tinsP='1.2'
+	RGPE='roAu*(tGPE/AGPE)'
+ 	CshGPE='(((PI2*eps0*epslys)/log(tinsP))*tGPE)'
**********************************************************************
*				GME CLEFT PARAMETERS
**********************************************************************
* dsealjGME	cleft width at the neuron-GME junction			[m]
* alphahead	head starting angle for resistor concentric rings model
* etawm		covering factor of the cross-sectional area of the solution.
* deltawm	surface overlapping coefficient: percentage of microel.
*               sensitive area covered by the neuron
*-----------------------------
+		dsealjGME='5n'		
+		alphahead='0.4294*PI'		
+		Etawm='1.0'
*
+		Rspreadhead='(((ros/PI2)*(dsealjGME/((dmhead/2.0)*((dmhead/2.0)+dsealjGME))))*Etawm)'
+		H='(htstalk+(dmstalk/2.0)+dsealjGME)'
+		Rspreadstalk='(((ros/(PI*H))*log((4.0*H)/(PI2*(dmstalk/2.0))))*Etawm)'
+		RspreadGME='(Rspreadhead+Rspreadstalk)'
*
+		deltawm='1.0'
+		rodj='(Ros/(PI2*dsealjGME))'
*
+		Rsealhead='(rodj*(log(tan(PI/4.0))-log(tan(alphahead/2.0))))'
+		Rsealrim='(rodj*(log((dmhead/2.0)/(dmstalk/2.0))))'	
+		Rsealstalk='(rodj*(htstalk/(dmstalk/2.0)))'
+		RsealGME='((Rsealhead+Rsealrim+Rsealstalk)*deltawm)'
**********************************************************************
*				GPE CLEFT PARAMETERS
**********************************************************************
* dsealjGPE	cleft width at the neuron-GPE junction		[m]
* etawp		covering factor of the cross-sectional area of the solution
* deltawp	surface overlapping coefficient: percentage of microel.
*		sensitive area covered by the neuron
*-----------------------------
+		dsealjGPE='90n'		
+		etawp='1.0'
* 
+		RspreadGPE='((ros/AGPE)*dsealjGPE*etawp)'	
*
+		deltawp='1.0'
+		RsealGPE='((Ros/dsealjGPE)*deltawp)'
**********************************************************************
*		GME & GPE ELECTROLYTE INTERFACE PARAMETERS
**********************************************************************
* Chg		capacitance of EDL at the interface microel-solution	[F]
* Rhg		leakage resistance 					[ohm]
*-----------------------------
+	ChgGME='5.0p'				
+	RhgGME='1500G'				
*
+	ChgGPE='1.13n'
+	RhgGPE='0.141Meg'
**********************************************************************
*  		GME & GPE PROTEIN-GLYCOCALYX EDL PARAMETERS
**********************************************************************
* Chd		capacitance of EDL at interface of glycocalyx-solution 	[F]
* Rhd		resistance of the EDL					[ohm]
* tglyco	glycocalyx thickness					[m]
* gammaw	correction factor
* fscale	scale factor
*-----------------------------
+	tglyco='100.0n'
*
+	gammawm='1'
+	ChdGME='((epsglyco*eps0)*(AGME/tglyco)*gammawm)'
+	RhdGME='(roglyco*(tglyco/AGME)*gammawm)'
*
+	fscale='2.5'
+	gammawp='fscale*(AGPE/AGME)'
+	ChdGPE='((epsglyco*eps0)*(AGPE/tglyco)*gammawp)'
+	RhdGPE='(roglyco*(tglyco/AGPE)*gammawp)'
**********************************************************************
* 		H-H NEURON MODEL PARAMETERS
**********************************************************************
* ENa		Sodium equilibrium potential			[V]
* EK		Potassium equilibrium potential			[V]
* El		Leakage equilibrium potential			[V]
* am, ah	activation and inactivation degrees of Na channels
* an		activation degree of K channels
* amc,amd,ahc,	rates by which channels switch from a closed to an open state
* ahd,anc,and	rates by which channels switch from a closed to an open state
* bmc,bmd,bhc,	rates for reverse
* bhd,bnc,bnd	rates for reverse
* fcm,fch,fcn	fractions of channels
* Cmem		neuron membrane capacitance per unit area	[F/cm**2] 
*-----------------------------
* Na Channel 
*-----------------------------
+	ENa=-115  gNa0=120
+	amc=0.1   am=25   amd=10    bmc=4   bmd=18
+	ahc=0.07  ahd=20  bhc=1     bh=30   bhd=10
+	cm=1      ch=1
+	fcm=1     fch=1	
*-----------------------------
* K Channel 
*-----------------------------
+	EK=12     gK0=36
+	anc=0.01  an=10  and=10  bnc=0.125  bnd=80
+	cn=1
+	fcn=1
*-----------------------------
* Leakage  Channel 
*-----------------------------
+	El=-10.613         Rl=3.333
**********************************************************************
*		Solving H-H equations system    
**********************************************************************
.SUBCKT		HH_NEURON   	1     20
*				In  Out-Cmem
*-----------------------------
* Na-circuit block               
*-----------------------------
* GNa=INa=gNa0*(m**3)*h*(Vmem-ENa)
GNa     1	2  value='gNa0*(fcm**3)*(V(21,2)**3)*fch*V(22,2)*V(1,2)'
*
* Gm1=-alpham*(1-m)
Gm1   21	2  value='-1*(amc*(V(1,20)+am)*(1-fcm*V(21,2)))/(exp((V(1,20)+am)/amd)-1)'
*
* Gm2=betam*m
Gm2     21	2  value='bmc*exp(V(1,20)/bmd)*fcm*V(21,2)'
*
* Gh=-[alphah*(1-h)-betah*h]
Gh     22	2  value='-1*(ahc*exp(V(1,20)/ahd)*(1-fch*V(22,2))-(bhc*fch*V(22,2))/(1+exp((V(1,20)+bh)/bhd)))'
*
ENa	2	20     value='ENa'
Cm	21	2      C='Cm'		IC=52.932m
Ch	22	2      C='Ch'		IC=596.121m
Rm	21	20     1G
Rh	22	20     1G
*-----------------------------
* K-circuit block                 
*-----------------------------
* GK=IK=gK0*(n**4)*(Vmem-EK)
GK  1  3  value='gK0*(fcn**4)*(V(31,3)**4)*V(1,3)'
*
* Gn1=-alphan*(1-n)
Gn1    31   3  value='-1*anc*(V(1,20)+an)*(1-fcn*V(31,3))/(exp((V(1,20)+an)/and)-1)'
*
* Gn2=betan*n
Gn2    31   3  value='bnc*exp(V(1,20)/bnd)*fcn*V(31,3)'
*
EK	3	20     value='EK'
Cn	31	3      C='Cn'		IC=317.677m 
Rn	31	20     1G
*-----------------------------
* Leakage-circuit block           
*-----------------------------
*
El	4	20     value='El'
Rl	1	4      R='Rl'
*-----------------------------
* Membrane capacitance  
*-----------------------------
Cmem 	1	20	1		IC=0
*-----------------------------
.ENDS		HH_NEURON
**********************************************************************
*		GME protein-glycocalyx stage model
**********************************************************************
.SUBCKT		GlycoGME	30	40
*                     		Cmem  	Rspread
Rhd		30   40		R='RhdGME'
Chd		30   40		C='ChdGME'
.ENDS		GlycoGME
**********************************************************************
*		GPE protein-glycocalyx stage model
**********************************************************************
.SUBCKT		GlycoGPE	50	60
*                     		Cmem  	Rspread
Rhd		50   60		R='RhdGPE'
Chd		50   60		C='ChdGPE'
.ENDS		GlycoGPE
**********************************************************************
*		GME Cleft stage model
**********************************************************************
.SUBCKT		CleftGME	70    	80	81
*                     		Glyco 	Rseal  	gnd
RspreadG	70    80	R='RspreadGME'
RsealG		80    81	R='RsealGME'
.ENDS		CleftGME
**********************************************************************
*		GPE Cleft stage
**********************************************************************
.SUBCKT		CleftGPE     90    	 100     	101
*                    	     Glyco 	Rseal  		gnd
RspreadM    	90   100	R='RspreadGPE'
RsealM	    	100  101	R='RsealGPE'
.ENDS		CleftGPE
**********************************************************************
*		GME-electrolyte stage model
**********************************************************************
.SUBCKT		GME-el_yte     170    		 180
* 		               Cleft    	 GME
Rhg	170   	180  	R='RhgGME'	
Chg	170  	180  	C='ChgGME'    
.ENDS		GME-el_yte     
**********************************************************************
*		GPE-electrolyte stage model
**********************************************************************
.SUBCKT		GPE-el_yte     190    		 200
* 		               Cleft    	 GPE
Rhg	190   	200  	R='RhgGPE'
Chg	190  	200  	C='ChgGPE'
.ENDS		GPE-el_yte
**********************************************************************
*		GME stage model
**********************************************************************
.SUBCKT	GME	220		230		231
* 	      Rhg-Chg 	       Cshunt		gnd
RmGME	220	230	R='RGME'
CshGME	230	231  	C='CshGME'
.ENDS		GME   
**********************************************************************
*		GPE stage model
**********************************************************************
.SUBCKT		GPE     240    	    250		251
* 		       Rhg-Chg 	   Cshunt	gnd
RmGPE	240   	250  	R='RGPE'	
CshGPE	250  	251  	C='CshGPE'    
.ENDS		GPE
**********************************************************************
*	  	NEURON coupled to GME	
**********************************************************************
XHHM1	11	0			HH_NEURON
RaxM1	11	22	R='Rax'
*
XHHM2	22	800			HH_NEURON
RaxM2	22	33	R='Rax'
*
XHHM3	33	801			HH_NEURON
RaxM3	33	44	R='Rax'
*
XHHM4	44	802			HH_NEURON
RaxM4	44	55	R='Rax'
*
XHHM5	55	0			HH_NEURON
**********************************************************************
*		 NEURON-GLYCOCALYX coupling 	
**********************************************************************
XNGLYM2    	  800     810         		GlycoGME
*       	  Cmem  Rspread
XNGLYM3   	  801     811        		GlycoGME
XNGLYM4   	  802     812        		GlycoGME
**********************************************************************
*		 GLYCOCALYX-CLEFT coupling
**********************************************************************
XNCLEM2    	810     820     0   		CleftGME
*       	Glyco  Rspread Rseal
XNCLEM3  	811     821     0   		CleftGME
XNCLEM4  	812     822     0   		CleftGME
**********************************************************************
*	  	CLEFT-GME coupling	
**********************************************************************
XNELECM2   	820	    	830    		GME-el_yte
* 		Cleft       	GME
XNELECM3  	821	    	831    		GME-el_yte
XNELECM4 	822	      	832 		GME-el_yte
**********************************************************************
*		GME
**********************************************************************
XNGMEM2   830	    	840 		0	GME
* 	  GME       	Cshunt		gnd
XNGMEM3   831	    	841 		0	GME
XNGMEM4   832	    	842 		0	GME
**********************************************************************
*	  	NEURON coupled to GPE	
**********************************************************************
XHHP1	111	0			HH_NEURON
RaxP1	111	122	R='Rax'
*
XHHP2	122	900			HH_NEURON
RaxP2	122	133	R='Rax'
*
XHHP3	133	901			HH_NEURON
RaxP3	133	144	R='Rax'
*
XHHP4	144	902			HH_NEURON
RaxP4	144	155	R='Rax'
*
XHHP5	155	0			HH_NEURON
**********************************************************************
*		 NEURON-GLYCOCALYX coupling 	
**********************************************************************
XNGLYP2    	  900     910         		GlycoGPE
*       	 Cmem   Rspread
XNGLYP3   	  901     911        		GlycoGPE
XNGLYP4   	  902     912        		GlycoGPE
**********************************************************************
*		 GLYCOCALYX-CLEFT coupling
**********************************************************************
XNCLEP2    	910     920     0   		CleftGPE
*       	Glyco  Rspread Rseal
XNCLEP3  	911     921     0   		CleftGPE
XNCLEP4  	912     922     0   		CleftGPE
**********************************************************************
*	  	CLEFT-GPE coupling	
**********************************************************************
XNELECP2   920	    	930    			GPE-el_yte     
* 	  Cleft       	GPE
XNELECP3   921	    	931    			GPE-el_yte     
XNELECP4   922	      	932 			GPE-el_yte     
**********************************************************************
*		GPE
**********************************************************************
XNGMEP2   930	    	940 		0	GPE
* 	  GPE       	Cshunt		gnd
XNGMEP3   931	    	941 		0	GPE
XNGMEP4   932	    	942 		0	GPE
**********************************************************************
.END
