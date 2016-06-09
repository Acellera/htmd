
from htmd.parameterization.step import Step, ScriptStep, StepSet, QMScriptStep

class NonpolarSet(StepSet):
    def __init__(self):
        super(NonpolarSet, self).__init__()

        self.add(ScriptStep( 0, "nonpolar", "input",
            []
             ))

        self.add(QMScriptStep( 10, "nonpolar", "preoptimize",
            [  "0:ligand.xyz" ]
             ))

        self.add(ScriptStep( 20,"nonpolar", "initial_parameters",
            [ "10:mol-opt.xyz"]
             ))

        self.add(ScriptStep( 30,"nonpolar", "generate_structure",
            [ "20:mol.rtf", "20:mol.prm", "10:mol-opt.xyz",
              "0:equiv-org.txt", "0:fixq.txt", "0:neutral.txt", "0:soft-dih-list.txt"
            ]

            ))

        self.add(QMScriptStep( 110,"nonpolar", "esp",
            [ "30:mol.xpsf", "10:mol-opt.xyz" , "30:mol.prm" ]
            ))

        self.add(ScriptStep( 120,"nonpolar", "fit_esp",
            [ "110:mol.xpsf", "110:mol-esp.dat", "110:mol-opt.xyz", "110:elem-list.txt", "110:mol-cgrid.out:cal-esp.out",
              "30:mol.prm", "30:mol.rtf", "30:equiv-org.txt", "30:neutral.txt", "30:fixq.txt" ]
         ))

        self.add(QMScriptStep( 150,"nonpolar", "fit_water",
            [ "30:equiv-org.txt", "30:mol.xpsf", "30:mol-wat.xpsf",
              "30:soft-dih-list.txt",
              "120:cg-list.txt", "120:elem-list.txt",
              "120:mol-opt.xyz",
              #"120:mol-esp.rtf", # this well get overwritten if there are acc/don
              "120:new-mol.xpsf:mol-esp.xpsf",
              "120:fit-mol.conf",
              "120:mol.prm", "120:mol.rtf", "120:mol-esp.dat", "120:fit-mol.conf",
              "120:para-check.dat:para-opt-start.dat" ]
            ))

#        self.add(ScriptStep( 210,"nonpolar", "detect_soft_torsion",
#            [ "150:mol.prm",
#              "150:new-mol.xpsf:mol.xpsf",
#              "150:mol-opt.crd" ]
#             ))


        self.add(ScriptStep( 220,"nonpolar", "pes",
            [ "150:soft-dih-list.txt", "150:mol-opt.xyz",
                "150:mol.prm", "150:new-mol.xpsf:mol.xpsf"
            ]
             ))

        self.add(QMScriptStep( 230,"nonpolar", "qm_1d_scan", [
            "120:elem-list.txt",
            "220:mol.prm", "220:mol.xpsf", "220:mol-opt.xyz", "220:soft-dih-list-new.txt:soft-dih-list.txt"
        ] ))

        self.add(ScriptStep( 240, "nonpolar", "1d-fitting", [
            "230:mol.prm", "230:mol.xpsf", "230:tor-1D-idx-*.dat", "230:soft-dih-list.txt",
            "150:mol-esp.rtf:mol-tor.rtf","150:mol-esp.rtf:mol.rtf", "230:mol-opt.xyz",
            "30:mol.rtf:org-mol.rtf", "30:mol.prm:org-mol.prm" , "30:mol.xpsf:org-mol.xpsf"
        ] ))

        self.add(ScriptStep( 250, "nonpolar", "refitcharges", [
            "150:para-check.dat:para-check-start.dat", "150:fit-mol.conf", "150:cg-list.txt",
            "240:soft-dih-list.txt", "240:saved-para.dat", "240:mol.xpsf", "240:mol.prm",
            "240:mol-tor.rtf:mol.rtf", "240:mol-opt.xyz", "150:mol-esp.dat", "150:mol-wat.xpsf",
            "150:E-mol-wat.txt"
        ] ))

        self.add(ScriptStep( 260, "nonpolar", "refittorsions", [
            "230:mol-opt.xyz",
            "250:new-mol.xpsf:mol.xpsf", "240:mol.prm", "250:soft-dih-list.txt",
            "240:tor-1D-idx-*.dat", "250:mol-esp.rtf:mol-tor.rtf", "250:mol-esp.rtf:mol.rtf"
        ] ))

        self.add(QMScriptStep( 270, "nonpolar", "rotamer", [
            "260:torsion-para-*.dat",
            "230:qm-1d-states.dat", "260:mol.prm", "260:mol.xpsf", "230:mol-opt.xyz", "230:n_rotamer.txt",
            "110:elem-list.txt"
        ] ))

        self.add(ScriptStep( 280, "nonpolar", "fit_rotamer", [
            "260:mm-tor-1D-idx-*.dat", "260:torsion-para-*.dat", "260:soft-dih-list.txt",
            "250:fit-mol.conf", "260:mol.prm", "260:mol.xpsf", "270:all-rotamer.dat",
            "260:mol-tor.rtf:mol.rtf", "260:mol-tor.rtf", "270:mol-opt.xyz" , "260:saved-para.dat"
        ] ))

        self.add(ScriptStep( 999, "nonpolar", "results", [
            "280:mol.rtf", "280:mol.prm", "270:QM-min.xyz:mol.xyz",
#            "260:fitting-1d-*.dat", "240:org-1d-*.dat",
            "280:1d-qm-mm-*.dat", "240:org-1d-*.dat",
            "260:soft-dih-list.txt"
        ], outputfiles=[ "mol.rtf", "mol.prm", "mol.xyz", "soft-dih-list.txt", "1d-qm-mm-*.dat", "org-1d-*.dat" ] ))

