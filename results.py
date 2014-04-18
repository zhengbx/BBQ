class FDmetResult(object):
    def __init__(self, FragmentResults):#, ClusterFactor):
        # sum up the total energy and total electron number        
        TotalEnergy = 0.
        TotalnElec = 0.
        for (Fragment, r) in FragmentResults:
            Factor = Fragment.factor #* ClusterFactor
            if r.EmbEnergy is not None:
                TotalEnergy += Factor * r.EmbEnergy
            TotalnElec += Factor * r.nEmbElec
        self.TotalElecEnergy = TotalEnergy
        self.TotalEnergy = None
        self.TotalnElec = TotalnElec
        self.FragmentResults = FragmentResults
        self.FullSystemFock = None
        self.FullSystemVcor = None
