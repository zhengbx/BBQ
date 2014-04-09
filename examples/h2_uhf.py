#!/usr/bin/env python

import BBQ
import h2_rhf

keys = {
'MFD': {'scf_solver': 'UHF'},
'FORMAT': {'verbose': 5,}
}

keys.update(h2_rhf.keys)

if __name__ == '__main__':
    import BBQ
    BBQ.main(keys)
