
  #include "compute_sL.H"
  

    alphaUnburned = oneKGMS * (-5.47712e-01*pow(fH,8)+1.44848e+00*pow(fH,7)-1.59433e+00*pow(fH,6) +9.49494e-01*pow(fH,5)-3.33735e-01*pow(fH,4)+7.10672e-02*pow(fH,3)-9.01563e-03*pow(fH,2)+5.74169e-04*fH+2.59170e-05);
    alphaUnburned *= pow(p/pref,0.200);

    const volScalarField& rhoUnburned = rhou;
    
    aUnburned = alphaUnburned / rhoUnburned;
    aUnburned.max(1e-10); // avoid division by zero in calculation of gCr

    dimensionedScalar eps0("eps0",dimensionSet(0,2,-3,0,0),1e-9);
    volScalarField epsilon = max(turbulence->epsilon(),eps0);

    dimensionedScalar k0("k0",dimensionSet(0,2,-2,0,0),1e-9);
    volScalarField k = max(turbulence->k(),k0);
    volScalarField up = uPrimeCoef*sqrt(2.0/3.0*k);
    up.max(1e-9);

    
    volScalarField lInt_Peters = 0.37*pow(up,3)/epsilon;	// Peters (CFX)

volScalarField b = scalar(1.0) - c;

muUnburned = oneKGMS * (-1.08702e-04*pow(fH,2)+1.70856e-04*fH+1.59871e-05);
p.max(1.0); // avoid p=0 in the following expression;
muUnburned *= pow(p/pref,0.205);
muUnburned.max(1e-8);	// bounding


// quenching:
    nu = thermo.mu()/rho;
    gCr = pow(sL,2)/aUnburned;
    epsilonCr = 15*nu*pow(gCr,2);
    epsilonCr.max(0.1); // can become zero on BCs
    volScalarField& lInt = lInt_Peters;
    volScalarField lKolmogorov = pow(pow(nu,3)/epsilon,0.25); 
    volScalarField sigma = max(0.28*Foam::log(lInt/lKolmogorov), scalar(1e-15));
    G = 0.5 * erfc(-pow(2*sigma,-0.5)*(log(epsilonCr/epsilon) + sigma/2.0));
    

// Calculate flame normal etc.
volVectorField n = -fvc::grad(c);

volScalarField mgb = mag(n);

dimensionedScalar dMgb = 1.0e-3*
    (b*c*mgb)().weightedAverage(mesh.V())
   /((b*c)().weightedAverage(mesh.V()) + SMALL)
  + dimensionedScalar("ddMgb", mgb.dimensions(), SMALL);

mgb += dMgb;

surfaceVectorField SfHat = mesh.Sf()/mesh.magSf();
surfaceVectorField nfVec = fvc::interpolate(n);
nfVec -= SfHat*(fvc::snGrad(b) + (SfHat & nfVec));
nfVec /= (mag(nfVec) + dMgb);
surfaceScalarField nf = (mesh.Sf() & nfVec);
n /= mgb;



// Calculate turbulent flame speed flux
surfaceScalarField phiSt = fvc::interpolate(rhoUnburned*sL*Xi)*nf;


// Xi Equation from XiFoam
b = scalar(1.0)-c;

volScalarField epsilonTilde = pow(uPrimeCoef, 3)*turbulence->epsilon();

volScalarField tauEta = sqrt(muUnburned/(rhoUnburned*epsilonTilde)); // Zeitmass
volScalarField Reta = up/
(
    sqrt(epsilonTilde*tauEta) + dimensionedScalar("1e-8", up.dimensions(), 1e-8)
);

    // Calculate mean and turbulent strain rates
    volVectorField Ut = U + sL*Xi*n;
    volScalarField Xisigmat = (n & n)*fvc::div(Ut) - (n & fvc::grad(Ut) & n);
    volScalarField Xisigmas =
      ((n & n)*fvc::div(U) - (n & fvc::grad(U) & n))/Xi
      + ( (n & n)*fvc::div(sL*n) - (n & fvc::grad(sL*n) & n))*(Xi + scalar(1))/(2*Xi);



// Calculate Xi according to the selected flame wrinkling model
if (XiModel == "fixed")
{
    // Do nothing, Xi is fixed!
}
else if (XiModel == "algebraic")
{
    // Simple algebraic model for Xi
    // with a linear correction function to give a plausible profile for Xi
    Xi == scalar(1) +
        (scalar(1) + (2*XiShapeCoef)*(c - scalar(0.5))) * XiCoef * sqrt(up/(sL + sLMin)) * Reta;
} 
else if (XiModel == "transport")
{

    // Calculate Xi flux
    surfaceScalarField phiXi =
      phiSt
      + fvc::interpolate(fvc::laplacian(turbulence->alphaEff(), c)/mgb)*nf
      + fvc::interpolate(rho)*fvc::interpolate(sL*(1.0/Xi - Xi))*nf;
      
  
    // Calculate Xi transport coefficients based on Gulders correlation
    // and DNS data for the rate of generation
    // with a linear correction function to give a plausible profile for Xi
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    volScalarField XiEqStar =
        scalar(1.001) + XiCoef*sqrt(up/(sL + sLMin))*Reta;

    volScalarField XiEq =
        scalar(1.001)
      + (scalar(1) + (2*XiShapeCoef)*(c - scalar(0.5)))
       *(XiEqStar - scalar(1.001));

    volScalarField Gstar = 0.28/tauEta;
    volScalarField XiR = Gstar*XiEqStar/(XiEqStar - scalar(1));
    volScalarField XiG = XiR*(XiEq - scalar(1.001))/XiEq;
  

    // Solve for the flame wrinkling
    fvScalarMatrix XiEqn
    (
        fvm::ddt(rho, Xi)
      + fvm::div(phi + phiXi, Xi, "div(phiXi,Xi)")
      - fvm::Sp(fvc::div(phiXi), Xi)
     ==
        rho*XiR			// XiR = Gstar bei hohen XiEqStar
      - fvm::Sp(rho*(XiR - XiG), Xi)	// = 0 bei hohen XiEq
      - fvm::Sp
        (
            rho*max
            (
                Xisigmat - Xisigmas,	// 
                dimensionedScalar("0", Xisigmat.dimensions(), 0)
            ),
            Xi
        )
    );

    XiEqn.relax();
    XiEqn.solve();

    // Xi should not fall below 1
    Xi.max(1.0);
    //Info<< "max(Xi) = " << max(Xi).value() << endl;
    //Info<< "max(XiEq) = " << max(XiEq).value() << endl;

}
else
{
    FatalError
        << args.executable() << " : Unknown Xi model " << XiModel
        << abort(FatalError);
}

Info<< "Combustion progress = "
    << 100*c.weightedAverage(mesh.V()).value() << "%" << endl;


sT = Xi*sL;
sT.correctBoundaryConditions(); // to avoid 0 on BCs, required for divisions by sT

