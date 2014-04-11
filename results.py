class FDmetResult(object):
    def __init__(self, FragmentResults, ClusterFactor):
        # sum up the total energy and total electron number        
        TotalEnergy = 0.
        TotalElec = 0.
        for (Fragment, r) in FragmentResults:
            Factor = Fragment.GetTotalFactor() * ClusterFactor
            TotalEnergy += Factor * r.EmbEnergy
            TotalElec += Factor * r.EmbElec
        self.TotalElecEnergy = TotalEnergy
        self.TotalEnergy = None
        self.TotalElec = TotalElec
        self.FragmentResults = FragmentResults
        self.FullSystemFock = None
        self.FullSystemVcor = None
