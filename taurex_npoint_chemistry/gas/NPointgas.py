"""
N-point free chemistry profile for TauREx 3.

Each molecule's volume mixing ratio (VMR) is defined at N discrete
pressure knot points and interpolated in log–log space across the full
atmospheric pressure grid.  Outside the knot range the boundary VMR is
held constant (flat extrapolation).
"""

import math

import numpy as np

from taurex.chemistry import Gas
from taurex.util import molecule_texlabel
from taurex.data.fittable import fitparam
from taurex.util import movingaverage
from taurex.exceptions import InvalidModelException


class InvalidChemistryException(InvalidModelException):
    """
    Exception that is called when atmosphere mix is greater
    than unity
    """
    pass


class NPointGas(Gas):
    """
    Single-molecule free chemistry profile with N pressure–VMR knots.

    Supply positive linear VMRs and pressures in Pa from surface to top.
    Each interior pressure has one corresponding abundance. Interpolation
    and moving-average smoothing operate on log10(VMR).
    """

    def __init__(
        self,
        molecule_name: str = "H2O",
        mix_ratio_surface: float = 1e-4,
        mix_ratio_top: float = 1e-4,
        P_surface: float = None,
        P_top: float = None,
        mix_ratios_points: list = [],
        pressure_points: list = [],
        limit_slope=9999999,
        smoothing_window=10,
    ):

        super().__init__('NPointGas', molecule_name)

        if not hasattr(mix_ratios_points, '__len__'):
            raise Exception('mix_ratios_points is not an iterable')

        if len(mix_ratios_points) != len(pressure_points):
            self.error('Number of mix_ratios_points != number of '
                       'pressure points')
            self.error('len(mix_ratios_points) = %s /= '
                       'len(pressure_points) = %s',
                       len(mix_ratios_points),
                       len(pressure_points))
            raise Exception('Incorrect_number of mix_ratios_points and pressure points')
        
        self.info('N-point gas profile is initialized for {}'.format(molecule_name))
        self.debug('Passed mix_ratios points %s', mix_ratios_points)
        self.debug('Passed pressure points %s', pressure_points)
        self.molecule_name = molecule_name
        self._mr_points = mix_ratios_points
        self._p_points = pressure_points
        self._mix_ratio_surface = mix_ratio_surface
        self._mix_ratio_top = mix_ratio_top
        self._P_surface = P_surface
        self._P_top = P_top
        self._limit_slope = limit_slope
        self._smooth_window = smoothing_window

        self.generate_pressure_fitting_params()
        self.generate_mix_ratio_fitting_params()
        
        # Surface and top volume mixing ratios
        def get_surface_vmr(self):
            return self._mix_ratio_surface

        def set_surface_vmr(self, value):
            self._mix_ratio_surface = value

        def get_top_vmr(self):
            return self._mix_ratio_top

        def set_top_vmr(self, value):
            self._mix_ratio_top = value

        self.add_fittable_param(
            f"{self.molecule_name}_surface",
            rf"${self.molecule_name}_{{\mathrm{{surf}}}}$",
            get_surface_vmr, set_surface_vmr,
            "log", False, [1e-12, 1e-1],
        )

        self.add_fittable_param(
            f"{self.molecule_name}_top",
            rf"${self.molecule_name}_{{\mathrm{{top}}}}$",
            get_top_vmr, set_top_vmr,
            "log", False, [1e-12, 1e-1],
        )
        
        
    # Function to generate the fitting parameters for the pressure points
    def generate_pressure_fitting_params(self):
        """Generates the fitting parameters for the pressure points
        These are given the name ``P_point(number)`` for example, if two extra
        pressure points are defined between the top and surface then the
        fitting parameters generated are ``P_point0`` and ``P_point1``
        """

        bounds = [1e0, 1e5]
        for idx, val in enumerate(self._p_points):
            point_num = idx+1
            param_name = f'{self.molecule_name}_' + 'P_point{}'.format(point_num)
            param_latex = f'{self.molecule_name}_' + '$P_{}$'.format(point_num)

            def read_point(self, idx=idx):
                return self._p_points[idx]

            def write_point(self, value, idx=idx):
                self._p_points[idx] = value

            fget_point = read_point
            fset_point = write_point
            self.debug('FGet_location %s', fget_point)
            default_fit = False
            self.add_fittable_param(param_name, param_latex, fget_point,
                                    fset_point, 'log', default_fit, bounds)
        
            
    # Function to generate the fitting parameters for the mixing ratio points
    def generate_mix_ratio_fitting_params(self):
        """Generates the fitting parameters for the mixing ratio points
        These are given the name ``mi_point(number)`` for example, if two extra
        mixing ratio points are defined between the top and surface then the
        fitting parameters generated are ``mi_point0`` and ``mi_point1``
        """

        bounds = [1e-10, 1.0]
        for idx, val in enumerate(self._mr_points):
            point_num = idx+1
            param_name = f'{self.molecule_name}_' + 'mr_point{}'.format(point_num)
            param_latex = f'{self.molecule_name}_' + '$mr_{}$'.format(point_num)

            def read_point(self, idx=idx):
                return self._mr_points[idx]

            def write_point(self, value, idx=idx):
                self._mr_points[idx] = value

            fget_point = read_point
            fset_point = write_point
            self.debug('FGet_location %s %s', fget_point, fget_point(self))
            default_fit = False
            self.add_fittable_param(param_name, param_latex, fget_point,
                                    fset_point, 'linear', default_fit, bounds)
        
    def check_profile(self, Ppt, MRp):
        
        if(any(Ppt[i] <= Ppt[i + 1] for i in range(len(Ppt)-1))): 
            self.warning('Chemistry profile is not valid - a pressure point is inverted')
            raise InvalidChemistryException

        if(any(abs((MRp[i+1]-MRp[i])/(np.log10(Ppt[i+1])-np.log10(Ppt[i]))) >= self._limit_slope for i in range(len(Ppt)-1))): 
            self.warning('Mixing ratio profile is not valid - profile slope too high')
            raise InvalidChemistryException
        
        
        
    def initialize_profile(self, nlayers=None, temperature_profile=None,
                           pressure_profile=None, altitude_profile=None):
        self.pressure_profile =  pressure_profile
        self.nlayers = nlayers
        
        
    #Creating the profile property to return the mixing ratio profile for the molecule
    @property
    def mixProfile(self):

        MRnodes = [self._mix_ratio_surface, *self._mr_points, self._mix_ratio_top]

        Psurface = self._P_surface
        if Psurface is None or Psurface < 0:
            Psurface = self.pressure_profile[0]

        Ptop = self._P_top
        if Ptop is None or Ptop < 0:
            Ptop = self.pressure_profile[-1]

        Pnodes = [Psurface, *self._p_points, Ptop]

        self.check_profile(Pnodes, MRnodes)

        smooth_window = self._smooth_window

        if all(value == MRnodes[0] for value in MRnodes):
            return np.ones_like(self.pressure_profile)*MRnodes[0]
        
        MR = np.interp((np.log10(self.pressure_profile[::-1])),
                       np.log10(Pnodes[::-1]), np.log10(MRnodes[::-1]))
        # smoothing MR-P profile
        wsize = int(self.nlayers*(smooth_window / 100.0))
        if (wsize % 2 == 0):
            wsize += 1
        MR_smooth = movingaverage(MR, wsize)
        border = int((len(MR) - len(MR_smooth))/2)

        log_vmr = MR.copy()
        if border > 0:
            log_vmr[border:-border] = MR_smooth

        # Return linear VMR in TauREx order, from surface to top.
        return 10.0 ** log_vmr[::-1]

        
        
        
    def write(self, output):
        chemistry = super().write(output)

        chemistry.write_scalar(f'{self.molecule_name}_vmr_surface', self._mix_ratio_surface)
        chemistry.write_scalar(f'{self.molecule_name}_vmr_top', self._mix_ratio_top)
        chemistry.write_array(f'{self.molecule_name}_vmr_points', np.array(self._mr_points))

        P_surface = self._P_surface
        P_top = self._P_top
        if not P_surface:
            P_surface = -1
        if not P_top:
            P_top = -1

        chemistry.write_scalar(f'{self.molecule_name}_P_surface', P_surface)
        chemistry.write_scalar(f'{self.molecule_name}_P_top', P_top)
        chemistry.write_array(f'{self.molecule_name}_pressure_points', np.array(self._p_points))
        chemistry.write_scalar('smoothing_window', self._smooth_window)

        return chemistry

    @classmethod
    def input_keywords(cls):
        """
        Return all input keywords
        """
        return ['chem_npoint', 'NPoint_chemistry', 'NP_chem',  'NPointGas']

