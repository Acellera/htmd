
from parameterize.step import Step, ScriptStep, StepSet

class MatchSet(StepSet):
    def __init__(self):
        super().__init__()

        self.add(ScriptStep( 0, "match", "input",
            []
             ))

        self.add(ScriptStep( 999, "match", "results", [
            "0:mol.rtf", "0:mol.prm", "0:min.xyz",
        ], outputfiles=[ "mol.rtf", "mol.prm", "min.xyz"  ] ))

