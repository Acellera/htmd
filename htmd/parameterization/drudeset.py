
from htmd.parameterization.step import Step, ScriptStep, StepSet

class DrudeSet(StepSet):
    def __init__(self):
        self.add(Step( 0, "input",
            [],
            [ "mol-opt.mol2",
            "mol-opt.pdb" ] ))

        self.add(ScriptStep( 10, "preoptimize",
            [ "0:mol-opt.mol2", "0:mol-opt.pdb" ],
            [] ))

        self.add(ScriptStep( 20, "initial_parameters",
            [ "10:org-mol.rtf"],
            [] ))

        self.add(ScriptStep( 30, "generate_structure",
            [ "20:mol.rtf", "20:mol.prm", "20:mol-opt.pdb" ],
            []
            ))

        self.add(ScriptStep( 110, "esp",
            [ "30:mol.xpsf", "30:mol.inp", "30:mol.crd", "30:mol.prm"],
            []
            ))

        self.add(ScriptStep( 120, "fit_esp",
            [ "110:mol.xpsf", "110:mol-esp.dat", "110:mol-opt.crd", "110:elem-list.txt",
              "30:mol.prm", "30:mol.rtf" ],
            [] ))

        self.add(ScriptStep( 150, "fit_water",
            [ "30:equiv-org.txt", "30:mol-opt.xpsf",
              "120:cg-list.txt", "120:elem-list.txt", "120:mol-opt.crd", "120:mol-esp.xpsf",
              "120:mol.prm", "120:mol.rtf", "120:mol-esp.dat", "120:fit-mol.conf", "120:para-check.dat" ],
            [] ))

        self.add(ScriptStep( 210, "detect_soft_torsion",
            [ "150:mol.prm", "150:new-mol.xpsf", "150:mol-opt.crd" ],
            [] ))


        self.add(ScriptStep( 220, "pes", [], [] ))

        self.add(ScriptStep( 230, "qm_1d_scan", [], [] ))

        self.add(ScriptStep( 240, "1d-fitting", [], [] ))

        self.add(ScriptStep( 250, "refitcharges", [], [] ))

        self.add(ScriptStep( 260, "refittorsions", [], [] ))

        self.add(ScriptStep( 270, "rotamer", [], [] ))

        self.add(ScriptStep( 280, "fit_rotamer", [], [] ))

        self.add(ScriptStep( 999, "results", [], [] ))

