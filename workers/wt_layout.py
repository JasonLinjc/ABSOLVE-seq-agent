"""Recognize WT original-library reads with the published fixed-region no-call rule."""
class WTRecognizer:
    def __init__(self, policy, targets):
        self.p = policy
        self.targets = targets
        assert policy['structures'] == {'original': ''} and policy['pbc_length'] == 16

    def recognize(self, sequence):
        constant = sequence[16:39]
        if len(sequence) < 48 or constant.count('N') > 5:
            return None, 'unrecognized_structure'
        if any(a != b and a != 'N' for a,b in zip(constant, self.p['constant'])):
            return None, 'unrecognized_structure'
        tbc = sequence[39:48]
        if tbc not in self.targets:
            return None, 'unknown_tbc'
        return ('original', '', 16, sequence[:16], tbc), None
