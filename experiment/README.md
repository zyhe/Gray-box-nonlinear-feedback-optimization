# Closed-Loop Optimization

## Main contents

`├── data/`: store convergence results  
`├── figure/`: store figures  
`├── params/`  
`│   ├── params_0.mat`: parameter file  
`├── tool_functions/`  
│   `├── dynamics.m`: implement system dynamics  
│   `├── extractProblem.m`: extract problem data  
│   `├── firstResponse.m`: first-order feedback opt.  
│   `├── generateProblem.m`: generate problems  
│   `├── grayBoxResponse.m`: gray-box feedback opt.  
│   `├── stochESResponse.m`: stochastic extremum seeking  
│   `├── zerothResponse.m`: model-free feedback opt.  
`├── closedLoopSim.m`: analyze the closed-loop performance  
`├── gridSearchParams.m`: explore different parameter combinations
`├── plotResult.m`: plot convergence results  

## How to Run

- Call `closedLoopSim.m`. The results will be stored in the folders of `data` and `figure`.  
  If `closedLoopSim(i)` is called, then `params_{i}.mat` is loaded, where `i = 0, 1, ...`. In this case, run `gridSearchParams.m` to generate various files of parameter combinations.
- `plotResult.m` allows adjusting the layout of figures based on the stored data in the folder of `data`.
