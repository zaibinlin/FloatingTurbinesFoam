    //U = UR;
        	        
    MRF.correctBoundaryVelocity(U);

    surfaceScalarField muEff
    (
        "muEff",
        mixture.muf()
      + fvc::interpolate(rho*turbulence->nut())
    );

    fvVectorMatrix UEqn
    (
        //original:
        //(1.0 + cPorField) / porosity * fvm::ddt(rho, U)
        //+ (1.0 + cPorField) / porosity * MRF.DDt(rho, U)
        
        //modify:
        1.0 / porosity * fvm::ddt(rho, U) + cPorField / porosity * fvm::ddt(rho, U)
        + 1.0 / porosity * MRF.DDt(rho, U) + cPorField / porosity * MRF.DDt(rho, U)
        
        //original:
        + 1.0/porosity * fvm::div(rhoPhi/porosityF, U)
        // + turbulence->divDevRhoReff(rho, U)
        - fvm::laplacian(muEff/porosityF , U)
        - 1.0/porosity * ( fvc::grad(U) & fvc::grad(muEff) )
        
        
        // Closure Terms
        + fvc::Sp(aPorField * pow(1.0 - porosity, 2) / pow(porosity,3)
            * mixture.mu() / pow(D50Field,2), U)   //UR
        + fvc::Sp(bPorField * rho * (1.0 - porosity) / pow(porosity,3) / D50Field
            * mag(U) *  //UR
        // Transient formulation
        (1.0 + useTransMask * 7.5 / KCPorField), U)  //UR
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
