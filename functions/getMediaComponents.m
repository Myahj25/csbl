

function media = getMediaComponents(modelnames)

    % List from Human1 Publication：Robinson, J. L. et al. 
    % An atlas of human metabolism. Sci Signal 13, doi:10.1126/scisignal.aaz1482 (2020).
    humangem.hamF12 = ["alanine";
                    "alpha-tocopherol";
                    "aquacob(III)alamin";
                    "arginine";
                    "asparagine";
                    "aspartate";
                    "biotin";
                    "choline";
                    "cysteine";
                    "Fe2+";
                    "folate";
                    "gamma-tocopherol";
                    "glucose";
                    "glutamate";
                    "glutamine";
                    "glycine";
                    "H2O";
                    "histidine";
                    "hypoxanthine";
                    "inositol";
                    "isoleucine";
                    "leucine";
                    "linoleate";
                    "linolenate";
                    "lipoic acid";
                    "lysine";
                    "methionine";
                    "nicotinamide";
                    "O2";
                    "pantothenate";
                    "phenylalanine";
                    "Pi";
                    "proline";
                    "pyridoxine";
                    "retinoate";
                    "riboflavin";
                    "serine";
                    "sulfate";
                    "thiamin";
                    "threonine";
                    "thymidine";
                    "tryptophan";
                    "tyrosine";
                    "valine"];
    
    humangem.dmem =  [ % Aminoacids
                    "glycine";
                    "alanine";
                    "arginine";
                    "asparagine";
                    "aspartate";
                    "cysteine"; 
                    "cystine" ;
                    "glutamate";
                    "histidine";
                    "isoleucine";
                    "leucine";
                    "lysine";
                    "methionine";
                    "phenylalanine";
                    "proline";
                    "serine";
                    "threonine";
                    "tryptophan";
                    "tyrosine";
                    "valine"; 
                    % Vitamins
                    "ascorbate";
                    "biotin";
                    "choline";
                    "pantothenate";
                    "folate";
                    "nicotinamide";
                    "pyridoxine";
                    "riboflavin";
                    "thiamin";
                    "cobamide-coenzyme"; % "vitamin b12"
                    "inositol";
                    %Inorganic Salts
                    "O2";
                    "H2O";
                    "H+";
                    "Pi";
                    "Ca2+";
                    "Cu2+";
                    "Fe2+";
                    "Fe3+";
                    "Na+";
                    "chloride";
                    "K+";
                    "Li+";
                    "Mg2+";
                    "NH3";
                    "HCO3-";
                    "formate";
                    "guanidine";
                    "iodide";
                    "H2S";
                    "NO";
                    "zinc";
                    "sulfate";
                    "nitrite";
                    "sulfite";
                    "methanol";
                    % Proteins
                    "Very Low Density Lipoprotein";
                    "Low Density Lipoprotein";
                    "Intermediate Density Lipoprotein";
                    "High Density Lipoprotein";
                    "fatty acid-chylomicron pool";
                    "fatty acid-uptake pool";
                    "fatty acid-VLDL pool";
                    "fatty acid pool";
                    % Trace Elements: Mn and V is missing.
                    "selenate";
                    "glucose";
                    "ethanolamine";
                    "hypoxanthine";
                    "linoleate";
                    "lipoic acid";
                    "putrescine";
                    "pyruvate";
                    "thymidine";
                    ];
    
    humangem.b27 = [ % Vitamins
                    "biotin";
                    "alpha-tocopherol";
                    "acetate";
                    % Proteins. missings: "catalase";  "human recombinant insulin"; "superoxide dismutase";
                    "fatty acid-chylomicron pool";
                    "fatty acid-uptake pool";
                    "fatty acid-VLDL pool";
                    "fatty acid pool";
                    "[apotransferin]";
                    %Other Components
                    "corticosterone";
                    "galactose";
                    "ethanolamine";
                    "cysteine"; "glycine"; "glutamate"; %glutathione (reduced)";
                    "L-carnitine";
                    "linolenate";
                    "linoleate";
                    "progesterone";
                    "putrescine";
                    "selenate";
                    "triiodothyronine"];
    

%%

    recon3d.hamF12 = ["ala_D";
                    "ala_L"; %"alanine";
                    "avite1"; % "alpha-tocopherol";
                    "aqcobal"; % "aquacob(III)alamin";
                    "arg_D";
                    "arg_L"; % "arginine";
                    "asn_L"; % "asparagine";
                    "asp_D"; 
                    "asp_L";% "aspartate";
                    "btn"; % "biotin";
                    "chol"; % "choline";
                    "cys_L"; % "cysteine";
                    "fe2"; % "Fe2+";
                    "fol"; % "folate";
                    "CE5655"; % "gamma-tocopherol";
                    "glc_D"; % "glucose";
                    "glu_L"; % "glutamate";
                    "gln_L"; % "glutamine";
                    "gly"; % "glycine";
                    "h2o"; % "H2O";
                    "his_L"; % "histidine";
                    "hxan"; % "hypoxanthine";
                    "inost"; % "inositol";
                    "ile_L"; % "isoleucine";
                    "leu_L"; % "leucine";
                    "lnlc"; % "linoleate";
                    "lnlnca"; 
                    "lnlncg"; % "linolenate";
                    "lipoate"; % "lipoic acid";
                    "lys_L"; % "lysine";
                    "met_L"; % "methionine";
                    "nmn"; 
                    "rnam";
                    "ncam"; % "nicotinamide";
                    "o2"; % "O2";
                    "pnto_R"; % "pantothenate";
                    "phe_L"; % "phenylalanine";
                    "pi"; % "Pi";
                    "pro_D"; % "proline";
                    "pro_L"; % "proline";
                    "pydxn"; % "pyridoxine";
                    "retn"; % "retinoate";
                    "ribflv"; % "riboflavin";
                    "ser_D"; % "serine";
                    "ser_L"; % "serine";
                    "so4"; % "sulfate";
                    "thm"; % "thiamin";
                    "thr_L"; % "threonine";
                    "thymd"; % "thymidine";
                    "trp_L"; % "tryptophan";
                    "tyr_L"; % "tyrosine";
                    "val_L"; % "valine"];
                    ];

    recon3d.dmem =  [ % Aminoacids
                    "gly"; % "glycine";
                    "ala_L";"ala_D"; % "alanine";
                    "arg_L"; % "arginine";
                    "asn_L"; % "asparagine";
                    "asp_L";% "aspartate";
                    "cys_L"; % "cysteine";
                    "Lcystin"; % "cystine" ;
                    "glu_L"; % "glutamate";
                    "his_L"; % "histidine";
                    "ile_L"; % "isoleucine";
                    "leu_L"; % "leucine";
                    "lys_L"; % "lysine";
                    "met_L"; % "methionine";
                    "phe_L"; % "phenylalanine";
                    "pro_D"; "pro_L"; % "proline";
                    "ser_D"; "ser_L"; % "serine";
                    "thr_L"; % "threonine";
                    "trp_L"; % "tryptophan";
                    "tyr_L"; % "tyrosine";
                    "val_L"; % "valine"];
                    % Vitamins
                    "ascb_L"; % "ascorbate";
                    "btn"; % "biotin";
                    "chol"; % "choline";
                    "pnto_R"; % "pantothenate";
                    "fol"; % "folate";
                    "ncam"; % "nicotinamide";
                    "pydxn"; % "pyridoxine";
                    "ribflv"; % "riboflavin";
                    "thm"; % "thiamin";
                    "cbl1"; "cbl2"; % "cobamide-coenzyme"; % "vitamin b12"
                    "inost"; % "inositol";
                    %Inorganic Salts
                    "o2"; % "O2";
                    "h2o"; % "H2O";
                    "h"; % "H+";
                    "pi"; % "Pi";
                    "ca2"; % "Ca2+";
                    "cu2"; % "Cu2+";
                    "fe2"; % "Fe2+";
                    "fe3"; % "Fe3+";
                    "na1"; % "Na+";
                    "cl"; % "chloride";
                    "k"; % "K+";
                    "M02382"; % "Li+";
                    "mg2"; % "Mg2+";
                    "nh4"; % "NH3";
                    "hco3"; % "HCO3-";
                    "for"; % "formate";
                    "M02035"; % "guanidine";
                    "i"; % "iodide";
                    "HC00250"; % "H2S";
                    "no"; % "NO";
                    "zn2"; % "zinc";
                    "so4"; % "sulfate";
                    "no2"; % "nitrite";
                    "so3"; % "sulfite";
                    "meoh"; % "methanol";
                    % Proteins
                    "vldl_hs"; % "Very Low Density Lipoprotein";
                    "ldl_hs"; % "Low Density Lipoprotein";
                    "idl_hs"; % "Intermediate Density Lipoprotein";
                    "hdl_hs"; % "High Density Lipoprotein";
                    "M01807"; % "fatty acid-chylomicron pool";
                    "M01819"; % "fatty acid-uptake pool";
                    "M01820"; % "fatty acid-VLDL pool";
                    "M01819"; % "fatty acid pool";
                    % Trace Elements: Mn and V is missing.
                    "sel"; % "selenate";
                    "glc_D"; % "glucose";
                    "etha"; % "ethanolamine";
                    "hxan"; % "hypoxanthine";
                    "lnlc"; % "linoleate";
                    "lipoate"; % "lipoic acid";
                    "ptrc"; % "putrescine";
                    "pyr"; % "pyruvate";
                    "thymd"; % "thymidine";
                    ];
    
    recon3d.b27 = [   % Vitamins
                    "btn"; % "biotin";
                    "avite1"; "avite2"; % "alpha-tocopherol";
                    "ac"; % "acetate";
                    % Proteins. missings: "catalase";  "human recombinant insulin"; "superoxide dismutase";
                    "M01807"; % "fatty acid-chylomicron pool";
                    "M01819"; % "fatty acid-uptake pool";
                    "M01820"; % "fatty acid-VLDL pool";
                    "M01819"; % "fatty acid pool";
                    "HC01944"; % "[apotransferin]";
                    %Other Components
                    "crtstrn"; % "corticosterone";
                    "gal"; % "galactose";
                    "etha"; % "ethanolamine";
                    "cys_L"; % "cysteine"; 
                    "gly"; % "glycine";
                    "glu_L"; % "glutamate";
                    "gthrd"; %glutathione (reduced)";
                    "crn"; % "L-carnitine";
                    "lnlncg"; % "linolenate";
                    "lnlc"; % "linoleate";
                    "prgstrn"; % "progesterone";
                    "ptrc"; % "putrescine";
                    "sel"; % "selenate";
                    "CE2754"; "HC02187"; % "triiodothyronine":
                    ];
    
    if find(contains(modelnames, 'humangem'))
        media.(modelnames(contains(modelnames, 'humangem'))) = unique([humangem.hamF12; humangem.dmem; humangem.b27]);
    end
    if find(contains(modelnames, 'recon'))
       media.(modelnames(contains(modelnames, 'recon'))) = unique([recon3d.hamF12; recon3d.dmem; recon3d.b27]);
    end


end