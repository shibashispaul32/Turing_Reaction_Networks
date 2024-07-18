Computational pipepline for manuscript **Widespread biochemical reaction networks enable Turing patterns without imposed feedback** 
Shibashis Paul, Joy Adetunji and Tian Hong (under review)

<br><br><br>
Install Anaconda or Miniconda and run the following code in terminal to create the necessary conda environment named 'turing_network': <br>
conda env create -f environment.yml

Estimated time to install all packages: 30min.

<br><br><br>
<img src='https://github.com/shibashispaul32/Turing_Reaction_Networks/blob/main/comp_details(README).png' width='600'>

Estimated time to run 10 sample simulations for reproducing Fig 2B: 30min.<br>
Estimated time to run numerical continuations for a small randomized parameter set (100): 2min.<br>
Estimated time to run numerical continuations for a large randomized parameter set (10000): 2.5hrs.<br>

<a href="https://visualpde.com/sim/?options=N4IghiBcAuBOCuBTANCAxgIynJq1RADcBLWYgE2IGcRVyBVeggFgDpmBGZgBgDYAmAKy0Q5AGpiCvVgHYAHAE4+QkeQDqaghwWtBg5rzkzhdAIqmC-VrwDM-ORxV1oUbq24f+eJpBBuPNsbGvABUYBhUABQcANRugiEASgCCAHIAIqkAlCJokr4c7nLcCnJhEdFxuklpmTl4mgXu-ArM5VGx8TUZ2bkWvm4CzDY2oeEdVQkpPfUgANYEcwD6AEIAvBwA3FTJG5vLu-4KXILcgiMKMgo2cgI2mwDSazIem2BLycnrcqwcvJwcF6XGwKBSCW5vD7rNw2QTyIyCP4eZhGOT8TYYD5fNaCX7FGw8GzcOT6biAjFQtZuDjgzi8YzXDg3Ypye5zAD2ADNOXsqOs2P8uJsRABbCBNcG8Wxg-iwvh-Ayi4gEfweUH8fj-FH8bSasEcEQAOygzFQAAcCABheBUaDskUiWA+EA7AC0BxC8BiHO5IX4IQA7q6HiEfZy-Z6QoRvasQpjPisIwBHb1c8P+pPB0NpiMBz0xxPxxP+gOO-Iulbu2PRsO5rO1-3wKMx5JhKG5x2NEDunMBmIh2tNmsfNsJkIp2uZge9z1VwvtoMe95fQOO-rd5bzseZpdY4vjnu+jP97NHwOekRUJPYBAoEDQcg33AgAOP3zwEQBi3vkAAXyAA">VisualPDE link for simulating example systems</a>

<br><br>
References to key parameter settings for models and algorithms<br>
| Parameters | Referencing text with rationale | Referencing code | Values/Ranges |
| :---:   | :---: | :---: | :---: |
| Distributions of ODE parameters (rate constants except &sigma;) | Supplementary text B.1.1 and Table S1   | run_ode_cont.py | Log-uniform with ranges listed in Table S1 |
| Control parameter for continuation &sigma;| Main text Methods   | run_ode_cont.py | 0-120 |
| Number of sampled parameter sets for ODE models| Main text Methods   | run_ode_cont.py | 10,000 (100 for testing) |
| Diffusion coefficients| Supplementary text B.1.2 and B.2.2   | ????? | Uniform 1-20 |
| Stepsize for dispersion &Delta;p| Supplementary text B.3 | script_scan_disp_all.py| 0.05 |
