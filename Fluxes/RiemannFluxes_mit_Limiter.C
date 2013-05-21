/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    RiemannFluxes

Description
    generic Godunov flux class

Author
    Oliver Borm  All rights reserved.
    Florian Ettner

\*---------------------------------------------------------------------------*/

#include "RiemannFluxes.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Flux>
Foam::RiemannFlux<Flux>::RiemannFlux
(
    const volVectorField& rhoU,
    const volScalarField& rho,
    const volScalarField& rhoE,
    const UPtrList<volScalarField>& rhoPassiveScalars,
    const eCombustionThermo& thermophysicalModel,
    const compressible::turbulenceModel& turbulenceModel
)
:
    mesh_(rho.mesh()),
    thermophysicalModel_(thermophysicalModel),
    turbulenceModel_(turbulenceModel),
    Npassive_(rhoPassiveScalars.size()),
   
    rho_(rho),
    rhoE_(rhoE),
    rhoPassiveScalar_(rhoPassiveScalars),
   
    rhoU_(rhoU),
    
    rhoFlux_
    (
        IOobject
        (
            "rhoFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(rhoU_) & mesh_.Sf()         // only initialisation
    ),
    rhoEFlux_
    (
        IOobject
        (
            "rhoEFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	rhoFlux_*linearInterpolate(rhoE_/rho_)         // only initialisation
    ),
    rhoUFlux_
    (
        IOobject
        (
            "rhoUFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoFlux_*linearInterpolate(rhoU_/rho_)         // only initialisation
    ),

    p_(rho_/thermophysicalModel.psi()),
    gradp_(fvc::grad(p_,"grad(scalarSlope)")),
    gradrho_(fvc::grad(rho_,"grad(scalarSlope)")),   
    gradrhoE_(fvc::grad(rhoE_,"grad(scalarSlope)")), 
    gradrhoU_(fvc::grad(rhoU_,"grad(USlope)")),
    a_(sqrt(thermophysicalModel_.Cp()/thermophysicalModel_.Cv() / thermophysicalModel_.psi())),
    grada_(fvc::grad(a_,"grad(aSlope)")),
    minimumLimiter
    (
        IOobject
        (
            "minimumLimiter",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("one", dimless, 1.0),
        zeroGradientFvPatchScalarField::typeName
    )   

{
      Info << "setting " << endl;
      Info << "Number of passive scalars: " << Npassive_ << endl;

      rhoScalarFlux_.setSize(Npassive_);
      gradrhoScalar_.setSize(Npassive_);
    
  
// these are variables  => & operator required, otherwise a copy is created
      //rhoScalarFlux_.set(0,&rhoCFlux_);	// results in a lost pointer in the destructor
      //rhoScalarFlux_.set(1,&rhoTauFlux_);
      //rhoScalarFlux_.set(2,&rhofHFlux_);

      forAll(rhoPassiveScalar_,i)
      {	    
	    rhoScalarFlux_.set
	    (
	      i, 
	      new surfaceScalarField
	      (
		IOobject
		(
		  "rhoScalarFlux",
		  mesh_.time().timeName(),
		  mesh_,
		  IOobject::NO_READ,
		  IOobject::NO_WRITE
		),
		rhoFlux_*linearInterpolate(rhoPassiveScalar_[i]/rho_)         // only initialisation
	      )
	    );	    
	    
      }	

	    forAll(rhoPassiveScalar_,i)
	    {	
	    gradrhoScalar_.set
	    (
	      i,
	      new volVectorField
	      (
		fvc::grad(rhoPassiveScalar_[i],"grad(scalarSlope)")
	      )
	    );
	    }
// these are variables  => & operator required, otherwise a copy is created
      //gradrhoScalar_.set(0,&gradrhoC_);
      //gradrhoScalar_.set(1,&gradrhoTau_);
      //gradrhoScalar_.set(2,&gradrhofH_);
      
          multidimLimiterSwitch=false;  // only intialisation
	  limiterName="vanAlbadaSlope"; // only intialisation
	  epsilon="5";
	  Konstant=0.05;	// MaInf (AUSM+up) or entropy fix parameter (Roe)
	  
    // read riemann solver coeffs
    if(mesh_.solutionDict().found("Riemann"))   // system/fvSolution
    {
        dictionary riemann = mesh_.solutionDict().subDict("Riemann");
        if (riemann.found("multidimLimiter"))
        {
            multidimLimiterSwitch = Switch(riemann.lookup("multidimLimiter"));
        }
        if (riemann.found("limiterName"))
        {
            limiterName = word(riemann.lookup("limiterName"));
        }
        if (riemann.found("epsilon"))
        {
            epsilon = word(riemann.lookup("epsilon"));
        }
        Konstant = riemann.lookupOrDefault("RoeKonstant",Konstant);
    }


      
      Info << "end of constructor " << endl;
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Flux>
void Foam::RiemannFlux<Flux>::update(Switch secondOrder)
{   
 
    Info << "Calculating Riemann fluxes, second order = " << secondOrder << endl;
    //Info << "rhoE max. = " << max(rhoE_).value() << endl;
    
    // Get face-to-cell addressing: face area point from owner to neighbour
    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    // Get the face area vector
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();

    const volVectorField& cellCenter = mesh_.C();
    const surfaceVectorField& faceCenter = mesh_.Cf();
   
    a_ = sqrt(thermophysicalModel_.Cp()/thermophysicalModel_.Cv() / thermophysicalModel_.psi()); // speed of sound
    p_ = rho_/thermophysicalModel_.psi();


    // 2nd order correction
    IStringStream blendingFactor(epsilon);
    IStringStream blendingFactorV(epsilon); // cannot put the same blendingFactor twice, don't know why

    // MultiDimensional Limiters have to be selected here; 1D Limiters can be selected at runtime
    
    //BarthJespersenSlopeMultiLimiter scalarLimiter(blendingFactor);
    //BarthJespersenSlopeMultiLimiter vectorLimiter(blendingFactorV);

    //VenkatakrishnanSlopeMultiLimiter scalarLimiter(blendingFactor);
    //VenkatakrishnanSlopeMultiLimiter vectorLimiter(blendingFactorV);
    
    constantSlopeMultiLimiter scalarLimiter(1.0);
    constantSlopeMultiLimiter vectorLimiter(1.0);


    // Calculate and store the cell based limiter
    volVectorField rhoULimiter
    (
        IOobject
        (
            "rhoULimiter",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("rhoULimiter", dimless, vector::one),
        zeroGradientFvPatchVectorField::typeName
    );

    PtrList<volScalarField> rhoScalarLimiter(Npassive_);
    forAll(rhoScalarLimiter,i)	// create limiters for all passive scalar fields
    {
      rhoScalarLimiter.set
      (i,new volScalarField 
	(//rhoScalar_[i].name()+"Limiter",      
	  IOobject
	  (
            rhoPassiveScalar_[i].name()+"Limiter",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
	  ),
	  mesh(),
	  dimensionedScalar("one", dimless, 1.0),
	  zeroGradientFvPatchScalarField::typeName
	)
      );
    }
    
    volScalarField rhoLimiter
    (
        IOobject
        (
            "rhoLimiter",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("one", dimless, 1.0),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField rhoELimiter
    (
        IOobject
        (
            "rhoELimiter",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("one", dimless, 1.0),
        zeroGradientFvPatchScalarField::typeName
    );
    
    volScalarField minScalarLimiter	// compute a minimum limiter valid for all transported scalars: rho, rhoE, all passive scalars
    (					// not implemented on the boundaries yet
        IOobject
        (
            "minScalarLimiter",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("one", dimless, 1.0),
        zeroGradientFvPatchScalarField::typeName
    );
    
    


    volScalarField pLimiter
    (
        IOobject
        (
            "pLimiter",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("one", dimless, 1.0),
        zeroGradientFvPatchScalarField::typeName
    );
    
    volScalarField aLimiter
    (
        IOobject
        (
            "aLimiter",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("one", dimless, 1.0),
        zeroGradientFvPatchScalarField::typeName
    );    
    
    

    if (secondOrder)
    {
	// compute gradients
	forAll(rhoPassiveScalar_,i) 
	{
	    gradrhoScalar_[i] = fvc::grad(rhoPassiveScalar_[i],"grad(scalarSlope)");	
	}
	gradrho_  = fvc::grad(rho_,"grad(scalarSlope)");
	gradrhoE_ = fvc::grad(rhoE_,"grad(scalarSlope)");
	
	grada_ = fvc::grad(a_,"grad(scalarSlope)");
	gradp_ = fvc::grad(p_,"grad(scalarSlope)");
	
	gradrhoU_ = fvc::grad(rhoU_,"grad(USlope)");  
	
        if (multidimLimiterSwitch)
        {
	    // allocate fields to store the min/max values
	    PtrList<volScalarField> rhoScalarMinValue(Npassive_);
	    PtrList<volScalarField> rhoScalarMaxValue(Npassive_);
	    forAll(rhoPassiveScalar_,i) 
	    {
	      rhoScalarMinValue.set(i,new volScalarField(rhoPassiveScalar_[i].name()+"min",rhoPassiveScalar_[i]));
	      rhoScalarMaxValue.set(i,new volScalarField(rhoPassiveScalar_[i].name()+"max",rhoPassiveScalar_[i]));     
	    }
	    
	    volScalarField rhoMinValue("rhoMin",rho_);
            volScalarField rhoMaxValue("rhoMax",rho_);

            volScalarField rhoEMinValue("rhoEMin",rhoE_);
            volScalarField rhoEMaxValue("rhoEMax",rhoE_);
	    
            volScalarField pMinValue("pMin",p_);
            volScalarField pMaxValue("pMax",p_);

	    volScalarField aMinValue("aMin",a_);
            volScalarField aMaxValue("aMax",a_);

    	    volVectorField rhoUMinValue("rhoUMin",rhoU_);
            volVectorField rhoUMaxValue("rhoUMax",rhoU_);

	    

            // for BJ /VK find min/max of each value for each variable
            // internal face
            forAll(owner, faceI)
            {
                label own = owner[faceI];
                label nei = neighbour[faceI];

                // min values
                pMinValue[own]   = min(pMinValue[own], p_[nei]);
                pMinValue[nei]   = min(pMinValue[nei], p_[own]);
                rhoUMinValue[own]   = min(rhoUMinValue[own], rhoU_[nei]);
                rhoUMinValue[nei]   = min(rhoUMinValue[nei], rhoU_[own]);
                rhoMinValue[own] = min(rhoMinValue[own], rho_[nei]);
                rhoMinValue[nei] = min(rhoMinValue[nei], rho_[own]);
                aMinValue[own]   = min(aMinValue[own], a_[nei]);
                aMinValue[nei]   = min(aMinValue[nei], a_[own]);		
		rhoEMinValue[own]   = min(rhoEMinValue[own], rhoE_[nei]);
                rhoEMinValue[nei]   = min(rhoEMinValue[nei], rhoE_[own]);		
		
                pMaxValue[own]   = max(pMaxValue[own], p_[nei]);
                pMaxValue[nei]   = max(pMaxValue[nei], p_[own]);
                rhoUMaxValue[own]   = max(rhoUMaxValue[own], rhoU_[nei]);
                rhoUMaxValue[nei]   = max(rhoUMaxValue[nei], rhoU_[own]);
                rhoMaxValue[own] = max(rhoMaxValue[own], rho_[nei]);
                rhoMaxValue[nei] = max(rhoMaxValue[nei], rho_[own]);
		aMaxValue[own]   = max(aMaxValue[own], a_[nei]);
                aMaxValue[nei]   = max(aMaxValue[nei], a_[own]);
                rhoEMaxValue[own]   = max(rhoEMaxValue[own], rhoE_[nei]);
                rhoEMaxValue[nei]   = max(rhoEMaxValue[nei], rhoE_[own]);

		forAll(rhoPassiveScalar_,i) 
		{
		  const volScalarField& rhoScalari=rhoPassiveScalar_[i];
		  volScalarField& rhoScalarMinValuei=rhoScalarMinValue[i];
		  volScalarField& rhoScalarMaxValuei=rhoScalarMaxValue[i];
	      
		  rhoScalarMinValuei[own]=min(rhoScalarMinValuei[own], rhoScalari[nei]);
		  rhoScalarMinValuei[nei]=min(rhoScalarMinValuei[nei], rhoScalari[own]);
		
		  rhoScalarMaxValuei[own]=max(rhoScalarMaxValuei[own], rhoScalari[nei]);
		  rhoScalarMaxValuei[nei]=max(rhoScalarMaxValuei[nei], rhoScalari[own]);
		}

            }

            // Update coupled boundary min/max values of conservative variables
            forAll(rho_.boundaryField(), patchi)
            {	
                const fvPatchScalarField& patchp 	= p_.boundaryField()[patchi];
                const fvPatchVectorField& patchrhoU 	= rhoU_.boundaryField()[patchi];
                const fvPatchScalarField& patchrho 	= rho_.boundaryField()[patchi];
		const fvPatchScalarField& patcha 	= a_.boundaryField()[patchi];
                const fvPatchScalarField& patchrhoE 	= rhoE_.boundaryField()[patchi];

                const unallocLabelList& faceCells = patchrho.patch().faceCells();

                if (patchrho.coupled()) // domain boundaries in parallel runs
                {
                    const scalarField patchpRight 	= patchp.patchNeighbourField();
                    const vectorField patchrhoURight 	= patchrhoU.patchNeighbourField();
                    const scalarField patchrhoRight 	= patchrho.patchNeighbourField();
                    const scalarField patchaRight 	= patcha.patchNeighbourField();
		    const scalarField patchrhoERight 	= patchrhoE.patchNeighbourField();
		    
                    forAll(patchrho, facei)
                    {
                        label own = faceCells[facei];

                        // min values at coupled boundary faces
                        pMinValue[own]   	=  min(pMinValue[own], patchpRight[facei]);
                        rhoUMinValue[own]   	=  min(rhoUMinValue[own], patchrhoURight[facei]);
                        rhoMinValue[own] 	=  min(rhoMinValue[own], patchrhoRight[facei]);
			aMinValue[own]   	=  min(aMinValue[own], patchaRight[facei]);
   		        rhoEMinValue[own]   	=  min(rhoEMinValue[own], patchrhoERight[facei]);
			
                        // max values at coupled boundary faces
                        pMaxValue[own]   	= max(pMaxValue[own], patchpRight[facei]);
                        rhoUMaxValue[own]   	= max(rhoUMaxValue[own], patchrhoURight[facei]);
                        rhoMaxValue[own] 	= max(rhoMaxValue[own], patchrhoRight[facei]);
			aMaxValue[own]   	= max(aMaxValue[own], patchaRight[facei]);
                        rhoEMaxValue[own]   	= max(rhoEMaxValue[own], patchrhoERight[facei]);	
                    }
                    //Info << rhoPassiveScalar_[0].boundaryField()[0].patchNeighbourField() << endl;
		    forAll(rhoPassiveScalar_,i) 
		    {
		      //Info << "boundaryField[patchi]: " << rhoPassiveScalar_[i].boundaryField()[patchi] << endl;
		      //Info << "boundaryField[patchi].neighbour: " << rhoPassiveScalar_[i].boundaryField()[patchi].patchNeighbourField() << endl;				
			const scalarField patchrhoScalarRighti = rhoPassiveScalar_[i].boundaryField()[patchi].patchNeighbourField();	
			
			volScalarField& rhoScalarMinValuei=rhoScalarMinValue[i];
			volScalarField& rhoScalarMaxValuei=rhoScalarMaxValue[i];
			forAll(patchrho,facei)
			{
			    label own = faceCells[facei];
			    rhoScalarMinValuei[own] = min(rhoScalarMinValuei[own],patchrhoScalarRighti[facei]);
			    rhoScalarMaxValuei[own] = max(rhoScalarMaxValuei[own],patchrhoScalarRighti[facei]);
			}
			
		    }
                    
                }
               
            }

	    volScalarField cellVolume
            (
                IOobject
                (
                    "cellVolume",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimVolume,
                zeroGradientFvPatchScalarField::typeName
            );

            cellVolume.internalField() = mesh().V();
            cellVolume.correctBoundaryConditions();

            // compute for each cell a limiter
            // Loop over all faces with different deltaR vector
            forAll(owner, faceI)  
            {
                label own = owner[faceI];
                label nei = neighbour[faceI];

                vector deltaRLeft  = faceCenter[faceI] - cellCenter[own];
                vector deltaRRight = faceCenter[faceI] - cellCenter[nei];

                // flux dummy
                scalar upwindFlux = 1.0;
		
		// Blazek p.167 letzte Zeile: jede Variable braucht ihren eigenen Limiter
		
                vector rhoUOwnerLimiterV = vectorLimiter.limiter
                (
                    cellVolume[own],
                    upwindFlux,
                    rhoUMaxValue[own]-rhoU_[own],
                    rhoUMinValue[own]-rhoU_[own],
                    gradrhoU_[own],
                    gradrhoU_[own],
                    deltaRLeft
                );
                vector rhoUNeighbourLimiterV = vectorLimiter.limiter
                (
                    cellVolume[nei],
                    upwindFlux,
                    rhoUMaxValue[nei]-rhoU_[nei],
                    rhoUMinValue[nei]-rhoU_[nei],
                    gradrhoU_[nei],
                    gradrhoU_[nei],
                    deltaRRight
                );

                scalar pOwnerLimiter = scalarLimiter.limiter
                (
                    cellVolume[own],
                    upwindFlux,
                    pMaxValue[own]-p_[own],
                    pMinValue[own]-p_[own],
                    gradp_[own],
                    gradp_[own],
                    deltaRLeft
                );
                scalar pNeighbourLimiter = scalarLimiter.limiter
                (
                    cellVolume[nei],
                    upwindFlux,
                    pMaxValue[nei]-p_[nei],
                    pMinValue[nei]-p_[nei],
                    gradp_[nei],
                    gradp_[nei],
                    deltaRRight
                );
		
		scalar aOwnerLimiter = scalarLimiter.limiter
                (
                    cellVolume[own],
                    upwindFlux,
                    aMaxValue[own]-a_[own],
                    aMinValue[own]-a_[own],
                    grada_[own],
                    grada_[own],
                    deltaRLeft
                );
                scalar aNeighbourLimiter = scalarLimiter.limiter
                (
                    cellVolume[nei],
                    upwindFlux,
                    aMaxValue[nei]-a_[nei],
                    aMinValue[nei]-a_[nei],
                    grada_[nei],
                    grada_[nei],
                    deltaRRight
                );

                scalar rhoOwnerLimiter = scalarLimiter.limiter
                (
                    cellVolume[own],
                    upwindFlux,		// dummy
                    rhoMaxValue[own]-rho_[own],	// deltaOneMax
                    rhoMinValue[own]-rho_[own],	// deltaOneMin
                    gradrho_[own],	// gradPhiP
                    gradrho_[own],	// dummy
                    deltaRLeft		// d = vector(faceCenter-cellCenter)
                );
                scalar rhoNeighbourLimiter = scalarLimiter.limiter
                (
                    cellVolume[nei],
                    upwindFlux,
                    rhoMaxValue[nei]-rho_[nei],
                    rhoMinValue[nei]-rho_[nei],
                    gradrho_[nei],
                    gradrho_[nei],
                    deltaRRight
                );

                scalar rhoEOwnerLimiter = scalarLimiter.limiter
                (
                    cellVolume[own],
                    upwindFlux,
                    rhoEMaxValue[own]-rhoE_[own],
                    rhoEMinValue[own]-rhoE_[own],
                    gradrhoE_[own],
                    gradrhoE_[own],
                    deltaRLeft
                );
                scalar rhoENeighbourLimiter = scalarLimiter.limiter
                (
                    cellVolume[nei],
                    upwindFlux,
                    rhoEMaxValue[nei]-rhoE_[nei],
                    rhoEMinValue[nei]-rhoE_[nei],
                    gradrhoE_[nei],
                    gradrhoE_[nei],
                    deltaRRight
                );
		
		scalarList rhoScalarOwnerLimiter(Npassive_);
		scalarList rhoScalarNeighbourLimiter(Npassive_);
		
		forAll(rhoPassiveScalar_,i)
		{		
		  const volScalarField& rhoScalari = rhoPassiveScalar_[i];
		  const volVectorField& gradrhoScalari = gradrhoScalar_[i];
		  const volScalarField& rhoScalarMinValuei = rhoScalarMinValue[i];
		  const volScalarField& rhoScalarMaxValuei = rhoScalarMaxValue[i];
		    rhoScalarOwnerLimiter[i] = scalarLimiter.limiter
		    (
                    cellVolume[own],
                    upwindFlux,
                    rhoScalarMaxValuei[own]-rhoScalari[own],
                    rhoScalarMinValuei[own]-rhoScalari[own],
                    gradrhoScalari[own],
                    gradrhoScalari[own],
                    deltaRLeft
		    );
		    
    		    rhoScalarNeighbourLimiter[i] = scalarLimiter.limiter
		    (
                    cellVolume[nei],
                    upwindFlux,
                    rhoScalarMaxValuei[nei]-rhoScalari[nei],
                    rhoScalarMinValuei[nei]-rhoScalari[nei],
                    gradrhoScalari[nei],
                    gradrhoScalari[nei],
                    deltaRRight
		    );
		}

                // find minimal limiter value in each cell
		
                rhoULimiter[own] = min(rhoULimiter[own], rhoUOwnerLimiterV);
                rhoULimiter[nei] = min(rhoULimiter[nei], rhoUNeighbourLimiterV);

                pLimiter[own]   = min(pLimiter[own], pOwnerLimiter);
                pLimiter[nei]   = min(pLimiter[nei], pNeighbourLimiter);

		aLimiter[own]   = min(aLimiter[own], aOwnerLimiter);
                aLimiter[nei]   = min(aLimiter[nei], aNeighbourLimiter);

                rhoLimiter[own] = min(rhoLimiter[own], rhoOwnerLimiter);
                rhoLimiter[nei] = min(rhoLimiter[nei], rhoNeighbourLimiter);

                rhoELimiter[own]   = min(rhoELimiter[own], rhoEOwnerLimiter);
                rhoELimiter[nei]   = min(rhoELimiter[nei], rhoENeighbourLimiter);
		
		minimumLimiter[own]=min(minimumLimiter[own],pLimiter[own]);
		minimumLimiter[own]=min(minimumLimiter[own],aLimiter[own]);
		minimumLimiter[own]=min(minimumLimiter[own],rhoLimiter[own]);
		minimumLimiter[own]=min(minimumLimiter[own],rhoELimiter[own]);
		minimumLimiter[nei]=min(minimumLimiter[nei],pLimiter[nei]);
		minimumLimiter[nei]=min(minimumLimiter[nei],aLimiter[nei]);
		minimumLimiter[nei]=min(minimumLimiter[nei],rhoLimiter[nei]);
		minimumLimiter[nei]=min(minimumLimiter[nei],rhoELimiter[nei]);
		
		
		forAll(rhoPassiveScalar_,i)
		{
		    volScalarField& rhoScalarLimiteri = rhoScalarLimiter[i];
		    rhoScalarLimiteri[own] = min(rhoScalarLimiteri[own], rhoScalarOwnerLimiter[i]);
		    rhoScalarLimiteri[nei] = min(rhoScalarLimiteri[nei], rhoScalarNeighbourLimiter[i]);
		    minimumLimiter[own]=min(minimumLimiter[own],rhoScalarLimiteri[own]); minimumLimiter[nei]=min(minimumLimiter[nei],rhoScalarLimiteri[nei]);
  
		}
		
            }

            // Update coupled boundary limiters
            forAll(rho_.boundaryField(), patchi)	
            {
                const fvPatchScalarField& patchrho = rho_.boundaryField()[patchi];

//                 const fvPatchVectorField& pCellCenter = cellCenter.boundaryField()[patchi];
                const unallocLabelList& faceCells = patchrho.patch().faceCells();

                if (patchrho.coupled())
                {
                    // cell and face centers
                    const vectorField delta =  patchrho.patch().delta();

                    forAll(patchrho, facei)
                    {
                        label own = faceCells[facei];

                        vector deltaRLeft  = delta[facei];

                        // flux dummy
                        scalar upwindFlux = 1.0;
			
                        vector rhoUOwnerLimiterV = vectorLimiter.limiter
                        (
                            cellVolume[own],
                            upwindFlux,
                            rhoUMaxValue[own]-rhoU_[own],
                            rhoUMinValue[own]-rhoU_[own],
                            gradrhoU_[own],
                            gradrhoU_[own],
                            deltaRLeft
                        );
                        scalar pOwnerLimiter = scalarLimiter.limiter
                        (
                            cellVolume[own],
                            upwindFlux,
                            pMaxValue[own]-p_[own],
                            pMinValue[own]-p_[own],
                            gradp_[own],
                            gradp_[own],
                            deltaRLeft
                        );
                        scalar aOwnerLimiter = scalarLimiter.limiter
                        (
                            cellVolume[own],
                            upwindFlux,
                            aMaxValue[own]-a_[own],
                            aMinValue[own]-a_[own],
                            grada_[own],
                            grada_[own],
                            deltaRLeft
                        );		
                        scalar rhoOwnerLimiter = scalarLimiter.limiter
                        (
                            cellVolume[own],
                            upwindFlux,
                            rhoMaxValue[own]-rho_[own],
                            rhoMinValue[own]-rho_[own],
                            gradrho_[own],
                            gradrho_[own],
                            deltaRLeft
                        );
                        scalar rhoEOwnerLimiter = scalarLimiter.limiter
                        (
                            cellVolume[own],
                            upwindFlux,
                            rhoEMaxValue[own]-rhoE_[own],
                            rhoEMinValue[own]-rhoE_[own],
                            gradrhoE_[own],
                            gradrhoE_[own],
                            deltaRLeft
                        );
			
			scalarList rhoScalarOwnerLimiter(Npassive_);
			forAll(rhoPassiveScalar_,i)
			{
			  const volScalarField& rhoScalari = rhoPassiveScalar_[i];
			  const volScalarField& rhoScalarMaxValuei = rhoScalarMaxValue[i];
			  const volScalarField& rhoScalarMinValuei = rhoScalarMinValue[i];
			  const volVectorField& gradrhoScalari = gradrhoScalar_[i];
			    rhoScalarOwnerLimiter[i] = scalarLimiter.limiter
			    (
			    cellVolume[own],
			    upwindFlux,
			    rhoScalarMaxValuei[own]-rhoScalari[own],
			    rhoScalarMinValuei[own]-rhoScalari[own],
			    gradrhoScalari[own],
			    gradrhoScalari[own],
			    deltaRLeft
			    );
			}
			
			
                        rhoULimiter[own] 	= min(rhoULimiter[own], rhoUOwnerLimiterV);
                        pLimiter[own] 		= min(pLimiter[own], pOwnerLimiter);
                        aLimiter[own] 		= min(aLimiter[own], aOwnerLimiter);
                        rhoLimiter[own] 	= min(rhoLimiter[own],rhoOwnerLimiter);
                        rhoELimiter[own] 	= min(rhoELimiter[own], rhoEOwnerLimiter);	
			    
			forAll(rhoPassiveScalar_,i)
			{
			    volScalarField& rhoScalarLimiteri = rhoScalarLimiter[i];
			    rhoScalarLimiteri[own] = min(rhoScalarLimiteri[own], rhoScalarOwnerLimiter[i]);		
			}
			    
                    }
                }
            }
        }
        else  // 1D limiters, limiterName from dictionary
        {
	  
	    //Info << "limiter Name: " << limiterName << endl;
	    
            IStringStream USchemeData(limiterName+"V");		// the streams need to be individual or each variable - don't know why
	    IStringStream pSchemeData(limiterName);
	    IStringStream aSchemeData(limiterName);
	    IStringStream rhoSchemeData(limiterName);
	    IStringStream eSchemeData(limiterName);	    

            surfaceScalarField rhoULimiterSurface =
                limitedSurfaceInterpolationScheme<vector>::New
                (
                    mesh_,
                    rhoFlux_,
                    USchemeData
                )().limiter(rhoU_);
	    
	    //Info << " rhoU " ;
	    
            surfaceScalarField pLimiterSurface =
                limitedSurfaceInterpolationScheme<scalar>::New
                (
                    mesh_,
                    rhoFlux_,
                    pSchemeData
                )().limiter(p_);
		
	    //Info << " / p " ;
	    
            surfaceScalarField aLimiterSurface =
                limitedSurfaceInterpolationScheme<scalar>::New
                (
                    mesh_,
                    rhoFlux_,
                    aSchemeData
                )().limiter(a_);
	    //Info << " / a " ;
		
            surfaceScalarField rhoLimiterSurface =
                limitedSurfaceInterpolationScheme<scalar>::New
                (
                    mesh_,
                    rhoFlux_,
                    rhoSchemeData
                )().limiter(rho_);
	    //Info << " / rho " ;
            surfaceScalarField rhoELimiterSurface =
                limitedSurfaceInterpolationScheme<scalar>::New
                (
                    mesh_,
                    rhoFlux_,
                    eSchemeData
                )().limiter(rhoE_);
		
	    //Info << " / rhoE " ;

	    PtrList<surfaceScalarField> rhoScalarLimiterSurface(Npassive_);

	    forAll(rhoPassiveScalar_,i)
	    {	      
	      IStringStream scalarSchemeData(limiterName);
	      rhoScalarLimiterSurface.set
	      (i,new surfaceScalarField
		(
		    limitedSurfaceInterpolationScheme<scalar>::New
		    (
                    mesh_,
                    rhoFlux_,
                    scalarSchemeData
		    )().limiter(rhoPassiveScalar_[i])
		)
	      );
	      	    //Info << " / scalar " ;
	      
	    }
	    Info << endl;
	    
	    
            
            forAll(owner, faceI)
            {
                label own = owner[faceI];
                label nei = neighbour[faceI];

                // find minimal limiter value in each cell
                rhoULimiter[own]   = min(rhoULimiter[own], rhoULimiterSurface[faceI]*vector::one);
                rhoULimiter[nei]   = min(rhoULimiter[nei], rhoULimiterSurface[faceI]*vector::one);
		
                pLimiter[own]   = min(pLimiter[own], pLimiterSurface[faceI]);
                pLimiter[nei]   = min(pLimiter[nei], pLimiterSurface[faceI]);

		aLimiter[own]   = min(aLimiter[own], aLimiterSurface[faceI]); 
                aLimiter[nei]   = min(aLimiter[nei], aLimiterSurface[faceI]);
		
                rhoLimiter[own] = min(rhoLimiter[own], rhoLimiterSurface[faceI]);
                rhoLimiter[nei] = min(rhoLimiter[nei], rhoLimiterSurface[faceI]);

		rhoELimiter[own]   = min(rhoELimiter[own], rhoELimiterSurface[faceI]); 
                rhoELimiter[nei]   = min(rhoELimiter[nei], rhoELimiterSurface[faceI]);

		forAll(rhoPassiveScalar_,i)
		{
		  volScalarField& rhoScalarLimiteri = rhoScalarLimiter[i];
		  const surfaceScalarField& rhoScalarLimiterSurfacei = rhoScalarLimiterSurface[i];
		    rhoScalarLimiteri[own] = min(rhoScalarLimiteri[own], rhoScalarLimiterSurfacei[faceI]);
		    rhoScalarLimiteri[nei] = min(rhoScalarLimiteri[nei], rhoScalarLimiterSurfacei[faceI]);    
		}
		    
            }

            // Update coupled boundary limiters
            forAll(rho_.boundaryField(), patchi)
            {
                const fvPatchScalarField& patch = rho_.boundaryField()[patchi];
                const unallocLabelList& faceCells = patch.patch().faceCells();

                if (patch.coupled())
                {
                    forAll(patch, facei)
                    {
                        label own = faceCells[facei];

                        pLimiter[own]   = min(pLimiter[own], pLimiterSurface[facei]);
                        rhoULimiter[own]   = min(rhoULimiter[own], rhoULimiterSurface[facei]*vector::one);
                        rhoLimiter[own] = min(rhoLimiter[own], rhoLimiterSurface[facei]);
                        aLimiter[own]   = min(aLimiter[own], aLimiterSurface[facei]);			
                        rhoELimiter[own]   = min(rhoELimiter[own], rhoELimiterSurface[facei]);
			forAll(rhoPassiveScalar_,i)
			{
			  volScalarField& rhoScalarLimiteri = rhoScalarLimiter[i];
			  const surfaceScalarField& rhoScalarLimiterSurfacei = rhoScalarLimiterSurface[i];
			  rhoScalarLimiteri[own] = min(rhoScalarLimiteri[own], rhoScalarLimiterSurfacei[facei]);
			
			}
                    }

                }
            }
        }	
	
	
    forAll(minScalarLimiter,face)
    {
	minScalarLimiter[face]=min(minScalarLimiter[face],rhoLimiter[face]);
	forAll(rhoPassiveScalar_,i)
	{
	  minScalarLimiter[face]=min(minScalarLimiter[face],rhoScalarLimiter[i][face]);
	}
    }

	
      
    } // end: second Order - condition


    scalarList rhoScalarOwns(Npassive_);
    scalarList rhoScalarNeis(Npassive_);
    scalarList rhoScalarFluxList(Npassive_);

    

  
    forAll(owner, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        vector deltaRLeft  = faceCenter[faceI] - cellCenter[own];
        vector deltaRRight = faceCenter[faceI] - cellCenter[nei];
	
	forAll(rhoScalarOwns,i)
	{
	  const volScalarField& rhoScalari = rhoPassiveScalar_[i];
	  const volVectorField& gradrhoScalari = gradrhoScalar_[i];
	  const surfaceScalarField& rhoScalarFluxi = rhoScalarFlux_[i];

	    //const volScalarField& rhoScalarLimiteri = rhoScalarLimiter[i];
	    //rhoScalarOwns[i] = rhoScalari[own] + secondOrder*rhoScalarLimiteri[own]*(deltaRLeft  & gradrhoScalari[own]); 
	    //rhoScalarNeis[i] = rhoScalari[nei] + secondOrder*rhoScalarLimiteri[nei]*(deltaRRight & gradrhoScalari[nei]); 
	    
	    rhoScalarOwns[i] = rhoScalari[own] + secondOrder*rhoLimiter[own]*(deltaRLeft  & gradrhoScalari[own]); 
	    rhoScalarNeis[i] = rhoScalari[nei] + secondOrder*rhoLimiter[nei]*(deltaRRight & gradrhoScalari[nei]); 
	    
	    //rhoScalarOwns[i] = rhoScalari[own] + secondOrder*min(rhoLimiter[own],rhoScalarLimiteri[own])*(deltaRLeft  & gradrhoScalari[own]); 
	    //rhoScalarNeis[i] = rhoScalari[nei] + secondOrder*min(rhoLimiter[nei],rhoScalarLimiteri[nei])*(deltaRRight & gradrhoScalari[nei]); 
	    
    	    //rhoScalarOwns[i] = rhoScalari[own] + secondOrder*minScalarLimiter[own]*(deltaRLeft  & gradrhoScalari[own]); 
	    //rhoScalarNeis[i] = rhoScalari[nei] + secondOrder*minScalarLimiter[nei]*(deltaRRight & gradrhoScalari[nei]); 

	  
	    rhoScalarFluxList[i] = rhoScalarFluxi[faceI];
	}
	
	//Info << faceI << ":\t" << rhoE_[own] << " --> " << rhoE_[own] + secondOrder*rhoELimiter[own]*(deltaRLeft  & gradrhoE_[own]) << "|\t" << rhoE_[nei] + secondOrder*rhoELimiter[nei]*(deltaRRight & gradrhoE_[nei]) << " <-- " << rhoE_[nei] << "\t (" << rhoELimiter[own] << " | " << rhoELimiter[nei] << ")" << endl;
	

	//Info << " evaluate face " << faceI << "\t" ;
        // calculate fluxes with reconstructed conservative variables at faces 
        Flux::evaluateFlux // compute new values of rhoUFlux_ and rhoScalarFluxList
        (
            rhoFlux_[faceI],
            rhoUFlux_[faceI],
            rhoEFlux_[faceI],
	    rhoScalarFluxList,

/*	 
	    p_[own] + secondOrder*pLimiter[own]*(deltaRLeft  & gradp_[own]), // try to express via rho_
            p_[nei] + secondOrder*pLimiter[nei]*(deltaRRight & gradp_[nei]),// 
   
	    rhoU_[own] + secondOrder*cmptMultiply(rhoULimiter[own],(deltaRLeft & gradrhoU_[own])), // bugfix ssaegeler
            rhoU_[nei] + secondOrder*cmptMultiply(rhoULimiter[nei],(deltaRRight& gradrhoU_[nei])), // bugfix ssaegeler
	 
	    rho_[own] + secondOrder*rhoLimiter[own]*(deltaRLeft  & gradrho_[own]), 
            rho_[nei] + secondOrder*rhoLimiter[nei]*(deltaRRight & gradrho_[nei]),
      
	    a_[own] + secondOrder*aLimiter[own]*(deltaRLeft  & grada_[own]),       
            a_[nei] + secondOrder*aLimiter[nei]*(deltaRRight & grada_[nei]),      

	    rhoE_[own] + secondOrder*rhoELimiter[own]*(deltaRLeft  & gradrhoE_[own]),	
	    rhoE_[nei] + secondOrder*rhoELimiter[nei]*(deltaRRight & gradrhoE_[nei]),
	    //rhoE_[own] + secondOrder*rhoLimiter[own]*(deltaRLeft  & gradrhoE_[own]),	// rhoLimiter
	    //rhoE_[nei] + secondOrder*rhoLimiter[nei]*(deltaRRight & gradrhoE_[nei]),	// rhoLimiter
*/
	 // one limiter for all variables:
	    p_[own] + secondOrder*minimumLimiter[own]*(deltaRLeft  & gradp_[own]),
            p_[nei] + secondOrder*minimumLimiter[nei]*(deltaRRight & gradp_[nei]),   
	    rhoU_[own] + secondOrder*minimumLimiter[own]*(deltaRLeft & gradrhoU_[own]), 
            rhoU_[nei] + secondOrder*minimumLimiter[nei]*(deltaRRight& gradrhoU_[nei]),

	    rho_[own] + secondOrder*minimumLimiter[own]*(deltaRLeft  & gradrho_[own]), 
            rho_[nei] + secondOrder*minimumLimiter[nei]*(deltaRRight & gradrho_[nei]),
      
	    a_[own] + secondOrder*minimumLimiter[own]*(deltaRLeft  & grada_[own]),       
            a_[nei] + secondOrder*minimumLimiter[nei]*(deltaRRight & grada_[nei]),      

	    rhoE_[own] + secondOrder*minimumLimiter[own]*(deltaRLeft  & gradrhoE_[own]),	
	    rhoE_[nei] + secondOrder*minimumLimiter[nei]*(deltaRRight & gradrhoE_[nei]),
	
	    rhoScalarOwns,
	    rhoScalarNeis,

            Sf[faceI],      // face vector
            magSf[faceI],   // face area
            Konstant	    // Roe Konstant (from dictionary, default 0.05)
        );
	
	forAll(rhoScalarOwns,i)	// write results back to rhoScalarFlux_ fields
	{
	  surfaceScalarField& rhoScalarFluxi = rhoScalarFlux_[i];
	    rhoScalarFluxi[faceI] = rhoScalarFluxList[i];
	    //Info << " [" << i << "]: \t" << rhoScalarFluxi.name() << ": " << rhoScalarFluxi[faceI] << "\t <-- " << rhoScalarFluxList[i] << endl;
	}
    }
    
   
    // Update boundary field and values
    forAll(rho_.boundaryField(), patchi)
    {      
        fvsPatchScalarField& patchrhoFlux  = rhoFlux_.boundaryField()[patchi];
        fvsPatchVectorField& patchrhoUFlux = rhoUFlux_.boundaryField()[patchi];
        fvsPatchScalarField& patchrhoEFlux = rhoEFlux_.boundaryField()[patchi];

        const fvPatchScalarField& patchp = p_.boundaryField()[patchi];
        const fvPatchVectorField& patchrhoU = rhoU_.boundaryField()[patchi];
        const fvPatchScalarField& patchrho = rho_.boundaryField()[patchi];
        const fvPatchScalarField& patcha = a_.boundaryField()[patchi];
	const fvPatchScalarField& patchrhoE = rhoE_.boundaryField()[patchi];
	

        const fvPatchVectorField& patchGradp = gradp_.boundaryField()[patchi];
        const fvPatchTensorField& patchGradrhoU = gradrhoU_.boundaryField()[patchi];
        const fvPatchVectorField& patchGradrho = gradrho_.boundaryField()[patchi];
	const fvPatchVectorField& patchGrada = grada_.boundaryField()[patchi];
	const fvPatchVectorField& patchGradrhoE = gradrhoE_.boundaryField()[patchi];

        const fvsPatchVectorField& patchSf = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& patchMagSf = magSf.boundaryField()[patchi];
        //const fvsPatchVectorField& patchDotX = dotX_.boundaryField()[patchi];

        const fvPatchVectorField& patchCellCenter = cellCenter.boundaryField()[patchi];

        const fvPatchScalarField& patchpLimiter = pLimiter.boundaryField()[patchi];
        const fvPatchVectorField& patchrhoULimiter = rhoULimiter.boundaryField()[patchi];
        const fvPatchScalarField& patchrhoLimiter = rhoLimiter.boundaryField()[patchi];
	const fvPatchScalarField& patchaLimiter = aLimiter.boundaryField()[patchi];
	const fvPatchScalarField& patchrhoELimiter = rhoELimiter.boundaryField()[patchi];	
	const fvPatchScalarField& patchminimumLimiter = minimumLimiter.boundaryField()[patchi];
	
        if (patchrho.coupled())
        {	//Info << "Patch " << patchi << " is coupled." << endl;
            const scalarField patchpLeft  = patchp.patchInternalField();
            const scalarField patchpRight = patchp.patchNeighbourField();

            const vectorField patchrhoULeft  = patchrhoU.patchInternalField();
            const vectorField patchrhoURight = patchrhoU.patchNeighbourField();

            const scalarField patchrhoLeft  = patchrho.patchInternalField();
            const scalarField patchrhoRight = patchrho.patchNeighbourField();

            const scalarField patchaLeft  = patcha.patchInternalField();
            const scalarField patchaRight = patcha.patchNeighbourField();

            const scalarField patchrhoELeft  = patchrhoE.patchInternalField();
            const scalarField patchrhoERight = patchrhoE.patchNeighbourField();

            // cell gradients
            const vectorField patchGradpLeft  = patchGradp.patchInternalField();
            const vectorField patchGradpRight = patchGradp.patchNeighbourField();

            const tensorField patchGradrhoULeft  = patchGradrhoU.patchInternalField();
            const tensorField patchGradrhoURight = patchGradrhoU.patchNeighbourField();

            const vectorField patchGradrhoLeft  = patchGradrho.patchInternalField();
            const vectorField patchGradrhoRight = patchGradrho.patchNeighbourField();

	    const vectorField patchGradaLeft  = patchGrada.patchInternalField();
            const vectorField patchGradaRight = patchGrada.patchNeighbourField();
	    
            const vectorField patchGradrhoELeft  = patchGradrhoE.patchInternalField();
            const vectorField patchGradrhoERight = patchGradrhoE.patchNeighbourField();
	    

            // cell limiters
            scalarField patchpLimiterLeft  = patchpLimiter.patchInternalField();
            scalarField patchpLimiterRight = patchpLimiter.patchNeighbourField();

            vectorField patchrhoULimiterLeft  = patchrhoULimiter.patchInternalField();
            vectorField patchrhoULimiterRight = patchrhoULimiter.patchNeighbourField();

            scalarField patchrhoLimiterLeft  = patchrhoLimiter.patchInternalField();
            scalarField patchrhoLimiterRight = patchrhoLimiter.patchNeighbourField();

            scalarField patchaLimiterLeft  = patchaLimiter.patchInternalField();
            scalarField patchaLimiterRight = patchaLimiter.patchNeighbourField();
	    
            scalarField patchrhoELimiterLeft  = patchrhoELimiter.patchInternalField();
            scalarField patchrhoELimiterRight = patchrhoELimiter.patchNeighbourField();
	    
	    scalarField patchminimumLimiterLeft  = patchminimumLimiter.patchInternalField();
	    scalarField patchminimumLimiterRight  = patchminimumLimiter.patchNeighbourField();
	    
            // cell and face centers
            const vectorField faceCenter =  patchrho.patch().Cf();
            const vectorField patchCellCenterLeft  =  patchCellCenter.patchInternalField();
            const vectorField patchCellCenterRight =  patchCellCenter.patchNeighbourField();

            forAll(patchrho, facei)
            {
		//Info << " face " << facei << endl;	      
                vector deltaRLeft  = faceCenter[facei] - patchCellCenterLeft[facei];
                vector deltaRRight = faceCenter[facei] - patchCellCenterRight[facei];

                // bound the limiters between 0 and 1 in order to prevent problems due to interpolation

                patchpLimiterLeft[facei]   = max(min(patchpLimiterLeft[facei],1.0),0.0);
                patchpLimiterRight[facei]  = max(min(patchpLimiterRight[facei],1.0),0.0);

                patchrhoULimiterLeft[facei]  = max(min(patchrhoULimiterLeft[facei],vector::one),vector::zero);
                patchrhoULimiterRight[facei] = max(min(patchrhoULimiterRight[facei],vector::one),vector::zero);

                patchrhoLimiterLeft[facei]   = max(min(patchrhoLimiterLeft[facei],1.0),0.0);
                patchrhoLimiterRight[facei]  = max(min(patchrhoLimiterRight[facei],1.0),0.0);

                patchaLimiterLeft[facei]   = max(min(patchaLimiterLeft[facei],1.0),0.0);
                patchaLimiterRight[facei]  = max(min(patchaLimiterRight[facei],1.0),0.0);

                patchrhoELimiterLeft[facei]   = max(min(patchrhoELimiterLeft[facei],1.0),0.0);
                patchrhoELimiterRight[facei]  = max(min(patchrhoELimiterRight[facei],1.0),0.0);
		
		patchminimumLimiterLeft[facei] = min(patchpLimiterLeft[facei],patchminimumLimiterLeft[facei]);
		patchminimumLimiterLeft[facei] = min(patchrhoLimiterLeft[facei],patchminimumLimiterLeft[facei]);
		patchminimumLimiterLeft[facei] = min(patchaLimiterLeft[facei],patchminimumLimiterLeft[facei]);
		patchminimumLimiterLeft[facei] = min(patchrhoELimiterLeft[facei],patchminimumLimiterLeft[facei]);
		
		patchminimumLimiterRight[facei] = min(patchpLimiterRight[facei],patchminimumLimiterRight[facei]);
		patchminimumLimiterRight[facei] = min(patchrhoLimiterRight[facei],patchminimumLimiterRight[facei]);
		patchminimumLimiterRight[facei] = min(patchaLimiterRight[facei],patchminimumLimiterRight[facei]);
		patchminimumLimiterRight[facei] = min(patchrhoELimiterRight[facei],patchminimumLimiterRight[facei]);
		

		
		scalarList patchrhoScalarOwns(Npassive_);
		scalarList patchrhoScalarNeis(Npassive_);
		scalarList patchrhoScalarFluxList(Npassive_);
		
		forAll(patchrhoScalarOwns,i)
		{
		      const fvPatchScalarField& patchrhoScalari = rhoPassiveScalar_[i].boundaryField()[patchi];
		      const scalarField patchrhoScalariLeft = patchrhoScalari.patchInternalField();	// references don't work here (first value of scalarField& patchrhoScalariLeft is always 0) - why ??? patchInternalField() returns a temporary field (?)
		      const scalarField patchrhoScalariRight= patchrhoScalari.patchNeighbourField();
		      
/*
		      Info << "i = " << i << endl;
		      Info << " patchrhoScalari: " << patchrhoScalari << endl;	      		      
		      Info << " patchrhoScalariLeft: " << patchrhoScalariLeft << endl;	      
		      Info << " patchrhoScalariRight: " << patchrhoScalariRight << endl;	      
*/

/*
		      const fvPatchScalarField& patchrhoScalariLimiter = rhoScalarLimiter[i].boundaryField()[patchi];
		      scalarField patchrhoScalariLimiterLeft  = patchrhoScalariLimiter.patchInternalField();  // rvalue, reference not possible
		      scalarField patchrhoScalariLimiterRight = patchrhoScalariLimiter.patchNeighbourField(); // rvalue, reference not possible
		      patchrhoScalariLimiterLeft[facei]   = max(min(patchrhoScalariLimiterLeft[facei],1.0),0.0);
		      patchrhoScalariLimiterRight[facei]  = max(min(patchrhoScalariLimiterRight[facei],1.0),0.0);
*/		      
		      const fvPatchVectorField& patchGradrhoScalari = gradrhoScalar_[i].boundaryField()[patchi];
		      const vectorField patchGradrhoScalariLeft  = patchGradrhoScalari.patchInternalField();
		      const vectorField patchGradrhoScalariRight = patchGradrhoScalari.patchNeighbourField();

		      // noch kein minimumLimiter implementiert
		      
		      // viel weniger Schwund/Produktion mit rhoLimiter statt eigenem scalarLimiter:
		      patchrhoScalarOwns[i]=patchrhoScalariLeft[facei]  + secondOrder * patchrhoLimiterLeft[facei] *(patchGradrhoScalariLeft[facei] & deltaRLeft);
		      patchrhoScalarNeis[i]=patchrhoScalariRight[facei] + secondOrder * patchrhoLimiterRight[facei]*(patchGradrhoScalariRight[facei] & deltaRRight);
		      


		      fvsPatchScalarField& rhoScalariFluxOnBoundary = rhoScalarFlux_[i].boundaryField()[patchi];
		      patchrhoScalarFluxList[i] = rhoScalariFluxOnBoundary[facei];	      
		}
	
                // Calculate fluxes at coupled boundary faces
                Flux::evaluateFlux
                (
                    patchrhoFlux[facei],
                    patchrhoUFlux[facei],
                    patchrhoEFlux[facei],
	 	    patchrhoScalarFluxList,
/*
		    patchpLeft[facei]  + secondOrder*patchpLimiterLeft[facei]*(patchGradpLeft[facei] & deltaRLeft),         
		    patchpRight[facei] + secondOrder*patchpLimiterRight[facei]*(patchGradpRight[facei] & deltaRRight),      
		 
                    patchrhoULeft[facei]  + secondOrder*cmptMultiply(patchrhoULimiterLeft[facei],(deltaRLeft & patchGradrhoULeft[facei])), 
                    patchrhoURight[facei] + secondOrder*cmptMultiply(patchrhoULimiterRight[facei],(deltaRRight & patchGradrhoURight[facei])),
		 
                    patchrhoLeft[facei]  + secondOrder*patchrhoLimiterLeft[facei] *(patchGradrhoLeft[facei]  & deltaRLeft),        
                    patchrhoRight[facei] + secondOrder*patchrhoLimiterRight[facei]*(patchGradrhoRight[facei] & deltaRRight),       
		 
                    patchaLeft[facei]  + secondOrder*patchaLimiterLeft[facei]* (patchGradaLeft [facei] & deltaRLeft),       
                    patchaRight[facei] + secondOrder*patchaLimiterRight[facei]*(patchGradaRight[facei] & deltaRRight),      
		 
                    patchrhoELeft[facei]  + secondOrder*patchrhoLimiterLeft[facei]* (patchGradrhoELeft [facei] & deltaRLeft), // rhoLimiter!
                    patchrhoERight[facei] + secondOrder*patchrhoLimiterRight[facei]*(patchGradrhoERight[facei] & deltaRRight),		 
*/
		    patchpLeft[facei]  + secondOrder*patchpLimiterLeft[facei]*(patchGradpLeft[facei] & deltaRLeft),         
		    patchpRight[facei] + secondOrder*patchpLimiterRight[facei]*(patchGradpRight[facei] & deltaRRight),      
		 
                    patchrhoULeft[facei]  + secondOrder*cmptMultiply(patchrhoULimiterLeft[facei],(deltaRLeft & patchGradrhoULeft[facei])), 
                    patchrhoURight[facei] + secondOrder*cmptMultiply(patchrhoULimiterRight[facei],(deltaRRight & patchGradrhoURight[facei])),
		 
                    patchrhoLeft[facei]  + secondOrder*patchrhoLimiterLeft[facei] *(patchGradrhoLeft[facei]  & deltaRLeft),        
                    patchrhoRight[facei] + secondOrder*patchrhoLimiterRight[facei]*(patchGradrhoRight[facei] & deltaRRight),       
		 
                    patchaLeft[facei]  + secondOrder*patchaLimiterLeft[facei]* (patchGradaLeft [facei] & deltaRLeft),       
                    patchaRight[facei] + secondOrder*patchaLimiterRight[facei]*(patchGradaRight[facei] & deltaRRight),      
		 
                    patchrhoELeft[facei]  + secondOrder*patchrhoLimiterLeft[facei]* (patchGradrhoELeft [facei] & deltaRLeft), // rhoLimiter!
                    patchrhoERight[facei] + secondOrder*patchrhoLimiterRight[facei]*(patchGradrhoERight[facei] & deltaRRight),		 


		    patchrhoScalarOwns,
		    patchrhoScalarNeis,
		 
                    patchSf[facei],       // face vector
                    patchMagSf[facei],    // face area
                    //patchDotX[facei],      // face velocity
                    Konstant
                );
		
		forAll(rhoScalarOwns,i)  // writing the resulting flux from the temporary to the permanent field
		{
		  fvsPatchScalarField& rhoScalariFluxOnBoundary = rhoScalarFlux_[i].boundaryField()[patchi];
		  rhoScalariFluxOnBoundary[facei] = patchrhoScalarFluxList[i];
		}
	
            }
        }
        else  // patch is not coupled
        {
		scalarList patchrhoScalarOwns(Npassive_);
		scalarList patchrhoScalarNeis(Npassive_);
		scalarList patchrhoScalarFluxList(Npassive_);
		
            forAll(patchrho, facei)
            {
	
		forAll(patchrhoScalarOwns,i)
		{
		      const fvPatchScalarField& rhoScalari = rhoPassiveScalar_[i].boundaryField()[patchi];
		      patchrhoScalarOwns[i] = rhoScalari[facei];
		      
      		      fvsPatchScalarField& rhoScalariFluxOnBoundary = rhoScalarFlux_[i].boundaryField()[patchi];
		      patchrhoScalarFluxList[i] = rhoScalariFluxOnBoundary[facei];
		}
		// calculate fluxes at not-coupled boundaries
                Flux::evaluateFlux
                (
                    patchrhoFlux[facei],
                    patchrhoUFlux[facei],
                    patchrhoEFlux[facei],
		    patchrhoScalarFluxList,

		    patchp[facei],        
                    patchp[facei],       
		 
                    patchrhoU[facei],     
                    patchrhoU[facei],     

		    patchrho[facei],      
                    patchrho[facei],       

                    patcha[facei],      
                    patcha[facei],      

                    patchrhoE[facei],
                    patchrhoE[facei],
		 
		    patchrhoScalarOwns,		// scalar list containing the values for all species
		    patchrhoScalarOwns,		// scalar list containing the values for all species

		    patchSf[facei],       
                    patchMagSf[facei],    
                    Konstant
                );
		
		forAll(rhoScalarOwns,i)  // write results from rhoScalarFluxList back to the correct position in rhoScalarFlux_
		{
	  	  fvsPatchScalarField& rhoScalariFluxOnBoundary = rhoScalarFlux_[i].boundaryField()[patchi];
		  rhoScalariFluxOnBoundary[facei] = patchrhoScalarFluxList[i];
		}
		
	    }
        }
    }
	

 
}

// ************************************************************************* //
