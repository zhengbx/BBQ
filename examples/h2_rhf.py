#!/usr/bin/env python

coord = [['H', 0.,0.,-1.],
          ['H', 0.,0., 1.]]

keys = {
'HAMILTONIAN': {'hammil': 'qc', },
'GEOMETRY': {'GeometryType' : 'xyz',
             'Coord': coord },
'BASIS': {'OrbBasis': 'STO-3G', },
'DMET': {'conv_threshold': 1e-3,
         'nSitesPerFragment': 1, },
'MFD': {'scf_solver': 'RHF'},
'FORMAT': {'verbose': 0,}
}

if __name__ == '__main__':
    import BBQ
    BBQ.main(keys)
