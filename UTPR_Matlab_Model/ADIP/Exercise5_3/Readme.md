# ATC Exercise Problem 5.3

To run all the Matlab scripts associated with exercise problem 5.3, you will need to have working installations of [_Yalmip_](https://yalmip.github.io/) and [_acados_](https://github.com/acados/acados).
_Yalmip_ just needs to be download and added to your Matlab path.

_Acados_ however requires an installed C compiler inside Matlab and is more complicated.
For Windows, there is an [installation guide](https://docs.acados.org/installation/index.html#windows-for-use-with-matlab) in the _acados_ documentation which you might follow.
Afterwards, you need to download [_CasADi_](https://web.casadi.org/) and install it into the `external` directory of _acados_, as described on the _acados_ [interface page](https://docs.acados.org/matlab_octave_interface/index.html#setup-casadi).
Finally, add the `ACADOS_INSTALL_DIR` environmental variable pointing to the main _acados_ directory.
You should expect this installation to take 20 to 30 minutes if you haven't done this before, but it does not need to be repeated each time you want to use _acados_ afterwards.
