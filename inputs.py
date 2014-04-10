#
# File: inputs.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import os, sys

__all__ = ['Input', 'dump_inputfile']
__doc__ = 'Handle the input keywords.'

class Input(object):
# To add more input keywords, create new entry in __doc__ and write the
# keyname, default value, etc. in the *dict format* as follows.  `value` is a
# function defined in _parse_strdict. Its argument `default` will be assigned
# to the default input dict; `allow` can be either a list of allowd values or
# a function to check the given value;  `limits` should be a tuple of
# (lower_bound, upper_bound) for the given value.
    ''' Parse inputs and do very simple check for inputs

    The keys and values for the input-dict are case insensitive.  But they are
    case sensitive when accesse via inputobj.xxx.yyy.
    Available functions:
        get_default_dict
        self_check
        sanity_check
        dump_keys

    Keys and default values
    -----------------------
    'HAMILTONIAN': {
        'Type': value(default='Hubbard',
                      allow=('Hubbard','qc','bcs','hubhol')),
            # type of hamiltonian that is to be solved
        't': None,
            # hopping matrix
        'U': None,
            # Hubbard
        'g': None,
            # eph interaction
        'W': None,
            # phonon energy
        'Occs': None,
            # filling
    },
    'GEOMETRY': {
        'lattice' : value(default='hub1d',
                          allow=('Hub1d','Hub2dSquare','Hub2dTilted','2dHoneyComb')),
            # lattice used in Hubbard lattice calc
        'GeometryType' : value(default='xyz',
                               allow=('xyz','ring','chain','grid',
                                      'chain-h2abs','chain-habs')),
        'Coord' : None,
        'UnitCell': None,
        'ClusterSize': None,
        'LatticeSize': None,
        'Fragments': None,
        'BoundaryCondition': None,
    },
    'BASIS': {
        'OrbBasis': None,
            # basis to be used in ab-initio calc
        'FitBasis': None,
        'BasisLibs': None,
    },
    'DMET': {
        'max_iter': value(default=20, allow=intType, limits=(1,100)),
            # maximum number of DMET self-consistency cycles
        'conv_threshold': value(default=1e-5, limits=(0,0.1)),
        'init_guess': None,
            # initial guess for dmet potential
            # AF = antiferromagnetic,
            # BCS,
            # RAND = random,
            # MAN = user specified,
            # None
        'diis_thr': value(default=1e-2,limits=(0,0.1)),
        'diis_start': value(default=2),
        'diis_dim': value(default=4),
        'frag_group': [],
            # a group of frag-list, every frag-list is a list of atom id which
            # indicates the atoms in the fragment
        'nSitesPerFragment': 1,
            #divides the entire basis into consecutive, equal length fragments
        'basidx_group': [],
            # a group of list, the list contains the id of additional basis
            # functions for the fragment
    },
    'FITTING': {
        'global_fit_dm': value(default=1, allow=(1,2,3,4)),
            # method for global fitting scheme
            # =1, neither FCI nor HF density matrix fixed
            # =2, fix FCI density matrix
            # =3, fix HF density matrix
            # =4, backwards fitting, no fixed DM
        'v_fit_domain': value(default=2, allow=(1,2,3,6,8)),
            # domain of fitting potential
            # =1, whole embedding basis
            # =2, impurity sites
            # =3, diagonal entry of whole embedding
            # =6, diagonal entry of impurity
            # =8, trace of impurity
        'dm_fit_constraint': False,
        'vfit_method': value(default='local',
                             allow=('local','global','oneshot')),
            # detailed fitting method
        'fitpot_damp_fac': value(default=0.5, limits=(0,1)),
        'vfit_init': 0,
            # initial guess for fitting method.  This key only affects the DM
            # fitting procedure, i.e. which fitting potential to start with.
            # It is different from DMET.init_guess.
    },
    'IMPSOLVER': {
        'imp_solver': value(default='FCI', allow=('FCI','CC')),
        'env_pot_for_fci': 0,
            # add the fitting potential to impurity solver or not
    },
    'MFD': {
        'scf_solver': value(default='RHF',
                            allow=('RHF','UHF','BCS')),
        'hf_follow_state': False,
            # mean field solver follows the initial guess or not
        'max_iter': 40,
            # maximum number of HF iterations
        'conv_thrsh': value(1e-6, limits=(0,0.1)),
            # convergence threshold for HF
    },
    'FORMAT': {
        'verbose': value(default=5, allow=range(6)),
            # print level, big value outputs more debug messages
        'sav_vfit': value(default=None, allow=stringType),
            # file to save the fitting potential
    },


    Examples
    --------
    >>> inp = Input({'DMET':{'max_iter':10},'MFD':{'scf_solver':'UHF'}})
    >>> print inp.MFD.max_iter
    40
    >>> inp.MFD.max_iter = 20
    Traceback (most recent call last):
        ...
    SyntaxError: Overwriting  MFD.max_iter  is not allowed
    >>> import copy; inp1 = copy.copy(inp)
    >>> inp1.MFD.max_iter = 20
    >>> print inp1.MFD.max_iter
    20
    >>> inp.self_check()
    True
    >>> inp.sanity_check('IMPSOLVER', 'imp_solver', 'MP2')
    Traceback (most recent call last):
        ...
    ValueError: MP2 is not one of ('FCI', 'CC')
    '''
    def __init__(self, input_dict={}, ofname=None):
        if ofname == None:
            self.output = sys.stdout
        else:
            self.output = open(ofname, 'w')

        # protect the input key.
        self._readonly = True
        self._input_dict = self.get_default_dict()
        self._checker = _parse_strdict(self._inpdict_in_doc(), True)
        self._phony = ['_checker', '_input_dict'] + self.__dict__.keys()

        # merge input keywords with the default keywords
        # input keywords are case insensitive
        all_mods = self._input_dict.keys()
        all_MODS = [i.upper() for i in all_mods]
        MOD2mod = dict(zip(all_MODS, all_mods))
        for mod, subdict in input_dict.items():
            if mod.upper() not in all_MODS:
                raise KeyError('module %s is not found' % mod)
            std_mod = MOD2mod[mod.upper()]
            all_keys = self._input_dict[std_mod].keys()
            all_KEYS = [i.upper() for i in all_keys]
            KEY2key = dict(zip(all_KEYS, all_keys))
            for k,v in subdict.items():
                if k.upper() in all_KEYS:
                    std_key = KEY2key[k.upper()]
                    self._input_dict[std_mod][std_key] = v
                else:
                    raise KeyError('key %s.%s is not found' % (mod, k))

        # phony keywords must be saved before _merge_input2self 
        self._merge_input2self(self._input_dict, self._readonly)

    def _merge_input2self(self, input_dict, read_only):
        # merge input_dict with the class attributes, so that the keyword xxx
        # can be accessed directly by self.xxx
        if read_only:
            # @property are added by function set_prop to RdOnlyClass instead
            # of _InlineClass
            class RdOnlyClass(_InlineClass): pass
            # use closure set_prop to hold key for get_x
            def set_prop(modname, key):
                def get_x(obj):
                    return getattr(obj, '__'+key)
                def set_x(obj, v):
                    raise SyntaxError('Overwriting  %s.%s  is not allowed'
                                      % (modname, key))
                setattr(RdOnlyClass, k, property(get_x, set_x))
            for mod, subdict in input_dict.items():
                self_mod = RdOnlyClass()
                setattr(self, mod, self_mod)
                for k,v in subdict.items():
                    setattr(self_mod, '__'+k, v)
                    set_prop(mod, k)
        else:
# TODO: use closure set_prop, as such _input_dict can be updated in set_x
            for mod, subdict in input_dict.items():
                self_mod = _InlineClass()
                setattr(self, mod, self_mod)
                for k,v in subdict.items():
                    setattr(self_mod, k, v)

    def _inpdict_in_doc(self):
        doc = _find_between(self.__doc__, \
                            '    -----------------------', '    Examples')
        return '{' + doc + '}'

    def _remove_phony(self):
        # Since the input keys have been merged to the class, remove the
        # intrinic member and functions to get input keys
        return filter(lambda s: not (s.startswith('_') or s in self._phony), \
                      self.__dict__.keys())

    # this __copy__ has bug if call twice since some key may diff from _input_dict
    # but if set_prop updates _input_dict, this __copy__ is OK
    #def __copy__(self):
    #    res = Input.__new__(Input)
    #    inpkeys = self._remove_phony()
    #    phonydic = [kv for kv in self.__dict__.items() \
    #                if kv[0] not in inpkeys]
    #    res.__dict__.update(phonydic)
    #    res._readonly = False
    #    res._merge_input2self(self._input_dict, False)
    #    return res
    def __copy__(self):
        res = Input.__new__(Input)
        inpkeys = self._remove_phony()
        for m in dir(self):
            self_mod = getattr(self, m)
            if m in inpkeys and isinstance(self_mod, _InlineClass):
                res_mod = _InlineClass()
                setattr(res, m, res_mod)
                for k in self_mod.__dict__.keys():
                    pk = k.strip('_')
                    v = getattr(self_mod, k)
                    res_mod.__dict__.update(((pk, v),))
            else:
                res.__dict__.update(((m, self_mod),))
        return res
    def __deepcopy__(self, memo):
        import copy
        shallow = self.__copy__().__dict__
        res = Input.__new__(Input)
        memo[id(self)] = res
        for k in self.__dict__.keys():
            setattr(res, k, copy.deepcopy(shallow[k], memo))
        return res

    def get_default_dict(self):
        '''Return the default input dict generated from __doc__'''
        return _parse_strdict(self._inpdict_in_doc())

    def self_check(self):
        '''Check sanity for all inputs'''
        for mod in self._remove_phony():
            if self._checker.has_key(mod):
                self_mod = self.__getattribute__(mod)
                for k in dir(self_mod):
                    if self._checker[mod].has_key(k):
                        if callable(self._checker[mod][k]):
                            self._checker[mod][k](getattr(self_mod, k))
                    else:
                        raise KeyError('key %s.%s is not found' % (mod, k))
            else:
                raise KeyError('module %s is not found' % mod)
        return True

    def sanity_check(self, modname, keyname, val):
        '''Check the given modname.keyname and its value.
        This check is case sensitive for modname, and keyname'''
        if self._checker.has_key(modname):
            if self._checker[modname].has_key(keyname):
                if callable(self._checker[modname][keyname]):
                    return self._checker[modname][keyname](val)
                else:
                    return True
            else:
                raise KeyError('key %s.%s is not found' % (modname, keyname))
        else:
            raise KeyError('module %s is not found' % modname)
        return True

    def dump_keys(self):
        self.output.write('**** modules and keywords ****\n')
        all_mods = self._remove_phony()
        for mod in sorted(all_mods):
            self_mod = getattr(self, mod)
            if isinstance(self_mod, _InlineClass):
                self.output.write('module: %s\n' % mod)
                for k in dir(self_mod):
                    v = getattr(self_mod, k)
                    self.output.write('    key: %s = %s\n' % (k,v))
            else:
                self.output.write('key: %s = %s\n' % (mod, self_mod))

def dump_inputfile(fileobj):
    try:
        filename = os.path.join(os.getcwd(), sys.argv[0])
        contents = open(filename, 'r').read()
        fileobj.write('**** input file is %s ****\n' % filename)
        fileobj.write(contents)
        fileobj.write('******************** input file end ********************\n')
    except:
        pass


class _InlineClass(object):
    def __dir__(self):
        return filter(lambda s: not s.startswith('_'), self.__dict__.keys())

def _find_between(s, start, end):
    ''' sub-strings between the start string and end string
    >>> _find_between('abcdefg', 'b', 'ef')
    'cd'
    '''
    s0 = s.find(start) + len(start)
    s1 = s.find(end)
    return s[s0:s1]

def _member(v, lst):
    if isinstance(v, str):
        return v.upper() in [i.upper() for i in lst \
                             if isinstance(i, str)]
    else:
        return v in lst

def _parse_strdict(docdict, checker=False):
    '''parse the __doc__ of Input class'''
    def stringType(s): return isinstance(s, str)
    def intType(s): return isinstance(s, int)
    if checker:
        # require the function 'value' in the __doc__ to generate the function
        # which can check the sanity for each key
        def value(default=None, allow=None, limits=None, **keywords):
            def check_sanity(v):
                if allow and callable(allow):
                    return allow(v)
                elif allow:
                    if not _member(v, allow):
                        raise ValueError('%s is not one of %s' % (v, allow))
                if limits:
                    if not (limits[0] <= v <= limits[1]):
                        raise ValueError('%s is not in %s' % (v, limits))
                return True
            return check_sanity
    else:
        def value(default=None, **keywords):
            return default
    dic = eval(docdict)
    return dic

if __name__ == '__main__':
    #doc = _find_between(Input.__doc__, \
    #                    '    -----------------------', '    Examples')
    #docdict = '{' + doc + '}'
    #print _parse_strdict(docdict)
    #print _parse_strdict(docdict, True)

    #import copy
    #inp = copy.copy(Input())
    #inp.MFD.scf_solver = 'UHF'
    #inp.self_check()
    #inp.xyz = 1
    #try:
    #    inp.self_check()
    #except KeyError:
    #    pass
    #inp.MFD.conv_thrsh = .9
    #try:
    #    inp.self_check()
    #except ValueError:
    #    pass
    #inp.dump_keys()

    import doctest
    doctest.testmod()
