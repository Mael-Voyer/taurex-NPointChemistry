# TauREx N-point chemistry

A TauREx 3 gas plugin for retrieving a molecule's volume mixing ratio (VMR) at several pressures. Other molecules in the same atmosphere can use their own profiles.

## Installation

From the repository root, with your TauREx environment active, run:

```bash
python -m pip install .
taurex --keywords gas
```

## Four-node water profile

Add this example to your TauREx parfile. Pressures are in Pa, ordered from deepest to highest altitude. Supply linear VMRs and one pressure per interior abundance. If you omit the endpoint pressures, the plugin uses the first and last atmospheric layer pressures.

```ini
[Chemistry]
chemistry_type = taurex
fill_gases = H2, He
ratio = 0.17

    [[H2O]]
    gas_type = NPointGas
    P_surface = 1e7
    P_top = 1e-1
    mix_ratio_surface = 1e-8
    mix_ratio_top = 1e-1
    pressure_points = 1e5, 1e3
    mix_ratios_points = 1e-5, 1e-4
    smoothing_window = 10

[Fitting]
H2O_surface:fit = True
H2O_surface:prior = "LogUniform(lin_bounds=(1e-12, 1.0))"
H2O_mr_point1:fit = True
H2O_mr_point1:prior = "LogUniform(lin_bounds=(1e-12, 1.0))"
H2O_mr_point2:fit = True
H2O_mr_point2:prior = "LogUniform(lin_bounds=(1e-12, 1.0))"
H2O_top:fit = True
H2O_top:prior = "LogUniform(lin_bounds=(1e-12, 1.0))"
```

The plugin interpolates in log pressure and log VMR. The optional moving average smooths log VMR across an odd number of layers. The default window is 10% of the grid; set it to zero for pure interpolation. The returned profile contains linear VMRs, including at the edges.

The example fits all four abundances at fixed pressures. Interior pressures are also available as fitting parameters, starting with `H2O_P_point1`. Choose priors that maintain decreasing pressure order. Endpoint pressures are configuration inputs.

## License

BSD 3-Clause. See [LICENSE](LICENSE).
