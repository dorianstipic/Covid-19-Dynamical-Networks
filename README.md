# Covid-19-Dynamical-Networks
This repo contains (old) codes and experiments regarding the scientific paper [Effects of quarantine disobedience and mobility restrictions on COVID-19 pandemic waves in dynamical networks](https://www.sciencedirect.com/science/article/pii/S0960077921005543) co-authored with my friend Mislav Bradač and two dear professors Tomislav Lipić and Boris Podobnik. 

There is also an earlier and simpler version of the paper completely in Croatian, of which I am the sole author. [Scenario analysis of the COVID-19 pandemic with graph simulations](https://drive.google.com/file/d/1EzaToeRD-5ksVCNBZd3sZklBeo3JKGSY/view?usp=sharing).

### Model Description
We proposed a plausible explanation of pandemic waves caused by changes in mobility restrictions and quarantine disobedience. The majority of the model parameters were determined using empirical data from medical sources. For simulation purposes the model was left with three (3) free parameters:


| Parameter | Description |
| ------ | ----------- |
| *N*   | Average 'unavoidable' daily contacts. |
| *p<sub>1</sub>* | Probability of mobility. |
| *p<sub>2</sub>* | Probability of quarantine measures disobedience. |



The model was applied to real data of the Croatian Covid-19 pandemic. 

*N* was approximated with the average household size, *p<sub>1</sub>* was extrapolated from the [Stringency index](https://ourworldindata.org/explorers/covid?uniformYAxis=0&hideControls=true&Metric=Stringency+index&Interval=7-day+rolling+average&Relative+to+Population=true&Color+by+test+positivity=false&country=USA%7EITA%7ECAN%7EDEU%7EGBR%7EFRA&hideControls=true) developed by the University of Oxford and *p<sub>2</sub>* was estimated using
Bayesian optimization applied to real mortality data observed in Croatia. 

### Install the Environment and Compile the code
To create the conda environment run the command
```
conda create --name your_env_name python==3.8.1 --file requirements.txt
```
If some packages fail to install (given the older version). It may be needed to install them separately with `pip`.

To compile the C++ code run `make`. Or,
 
```
g++ -O2 -std=c++11 -o model_cluster_trip_v2 model_cluster_trip_v2.cpp
```

To run a simulation,

```
model_cluster_trip_v2 config_cluster_trip_v2_example.json 0 > output.json
```

### Files Description

- `model_cluster_trip_v2.cpp` is where the model is implemented.
- `real_grid_search.py` runs grid search for the given input parameters. See readme_run scripts/real_grid_search_k2.5_mu5_base.sge for exec details.
- `crit_bound_grid.py` finds the 'critical boundary' of the 'catastrophe zone' of healthcare. See readme_run scripts/crit_bound_grid_k2.5_mu5_base.sge for exec details.
- `covid19plots.ipynb` is an old experimental notebook. Is left here for legacy reasons, the notebook is poorly written and should be used as **read-only**.  The code is broken and is too complex (and too bad) to be fixed at this stage.
- `covid19plots cleaner.ipynb` is a recent 'cleaner' version of the nb which can be run and shows some experiments.
- `plot.py` is a vanilla plot.
- `config*.json` are configs that are used to run the C++ code, see how to run a simulation above.
- `R0simul.py` contains some simulations regarding the basic reproduction number R~0~.
- `random_sampling_mortality_r.py`, `random_sampling_plot.py` and `test.py` are some codes related to fat-tail testing. The codes are not particularly important, just a bunch of simulations. The output is `out2.json`. Last section of the notebooks contains some related stuff and analysis.
- `\inputs` folder containing real Croatian Covid-19 data.
- `\log bayes` folder containing inputs and outputs of the Bayesian optimization.
- `\readme_run scripts` folder containing SGE submit scripts.
- `\outputs` folder containing experiment outputs and images (from regular python scripts). Mostly relevant for the [early version](https://drive.google.com/file/d/1EzaToeRD-5ksVCNBZd3sZklBeo3JKGSY/view?usp=sharing).
- `\outputsBilj` folder containing more refined experiments and image (mostly from jupyter notebooks).
- `\tmp` folder containing some vanilla examples.
